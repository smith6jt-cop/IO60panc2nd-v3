/*
 * check_measurements.groovy
 *
 * Diagnostic. Reports per image whether detections already have intensity
 * measurements from a prior run (InstanSeg with "Make measurements" enabled,
 * or stage 04, or any other source).
 *
 * Run via: Run -> Run for project. Output goes to the log and is appended to
 * results/measurement_audit.csv for review across all images.
 *
 * What it checks on the first detection of each image:
 *   - presence of a sentinel intensity column ("Cell: DAPI: Mean" by default)
 *   - presence of a shape column ("Cell: Area um^2" or similar)
 *   - total measurement count
 * The sentinel column is a representative probe; if it is present, the full
 * set from ObjectMeasurements.addIntensityMeasurements is almost certainly
 * present too, since that function writes all combinations atomically.
 */

evaluate(new File(buildFilePath(PROJECT_BASE_DIR, 'scripts', 'lib', 'common.groovy')))

// ==================== PARAMETERS ====================
String SENTINEL_CHANNEL    = 'DAPI'   // any channel guaranteed to be in the panel
String SENTINEL_COMPARTMENT = 'Cell'
String SENTINEL_STAT       = 'Mean'
// ====================================================

def imageName = (getProjectEntry()?.getImageName() ?: 'unknown').replaceAll(',', '_')
def detections = getDetectionObjects()
def n = detections.size()

if (n == 0) {
    println "[audit] ${imageName}: NO DETECTIONS"
    return
}

def sample = detections[0]
def ml = sample.getMeasurementList()

// Pull the full list of measurement names. The API differs slightly across
// versions; try the cleanest path first.
def names = []
try {
    names = ml.getNames() as List
} catch (Throwable t1) {
    try {
        names = ml.getMeasurementNames() as List
    } catch (Throwable t2) {
        // Fallback: probe the measurements map view
        try { names = sample.measurements.keySet() as List }
        catch (Throwable t3) { names = [] }
    }
}

def sentinelKey = "${SENTINEL_COMPARTMENT}: ${SENTINEL_CHANNEL}: ${SENTINEL_STAT}"
double sentinelVal = getMeas(sample, sentinelKey)
boolean hasSentinel = !Double.isNaN(sentinelVal)

// Look for any intensity-style column to handle the case where the sentinel
// channel happens to be missing or named differently.
def intensityPattern = ~'^(Cell|Cytoplasm|Membrane|Nucleus): [^:]+: (Mean|Min|Max|Std\\.?Dev|Median)$'
def intensityLike = names.findAll { intensityPattern.matcher(it).matches() }
def shapeLike = names.findAll {
    def s = it.toLowerCase()
    s.contains('area') || s.contains('circularity') || s.contains('diameter') || s.contains('solidity')
}

println "[audit] ${imageName}"
println "[audit]   detections: ${n}"
println "[audit]   total measurements on sample: ${names.size()}"
println "[audit]   sentinel '${sentinelKey}' present: ${hasSentinel}" +
        (hasSentinel ? " (value=${sentinelVal})" : "")
println "[audit]   intensity-style columns found: ${intensityLike.size()}"
println "[audit]   shape-style columns found: ${shapeLike.size()}"

// Decision rule: if sentinel is present AND there are at least ~50 intensity
// columns, treat as "already measured". 50 is a conservative threshold;
// 59 channels x 4 compartments x 4 stats would be 944, but a partial run
// (e.g. one compartment, one stat) could be far fewer and still be valid.
String verdict
if (hasSentinel && intensityLike.size() >= 50) {
    verdict = 'ALREADY_MEASURED'
} else if (hasSentinel || intensityLike.size() >= 50) {
    verdict = 'PARTIAL'
} else {
    verdict = 'NOT_MEASURED'
}
println "[audit]   verdict: ${verdict}"

// Append to a single project-level CSV for cross-image review
def outDir = buildFilePath(PROJECT_BASE_DIR, 'results'); mkdirs(outDir)
def csv = new File(buildFilePath(outDir, 'measurement_audit.csv'))
boolean isNew = !csv.exists()
csv.withWriterAppend('UTF-8') { w ->
    if (isNew) w.write("image,n_detections,n_measurements,sentinel_present,sentinel_value,n_intensity_cols,n_shape_cols,verdict\n")
    w.write("${imageName},${n},${names.size()},${hasSentinel},${hasSentinel ? sentinelVal : ''},${intensityLike.size()},${shapeLike.size()},${verdict}\n")
}
println "[audit]   appended to ${csv}"