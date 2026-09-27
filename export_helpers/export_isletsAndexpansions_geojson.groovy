/*
 * export_islets_to_geojson.groovy
 *
 * Replicates the manual QuPath workflow:
 *   1. select every annotation classed 'Islet' and 'IsletExpanded'
 *   2. File -> Export objects -> GeoJSON, as a FeatureCollection,
 *      measurements INCLUDED, no compression.
 *
 * Differences from export_class_to_geojson.groovy:
 *   - Handles a LIST of classes in one combined export (not a single class).
 *   - Writes a true FeatureCollection ({"type":"FeatureCollection",...}) via
 *     QuPath's own exporter. GsonTools(...).toJson(list) emits a bare JSON
 *     array instead, which does NOT match the GUI export.
 *   - Includes measurements, uncompressed, compact JSON (matches the example).
 *   - Does NOT delete the annotations afterwards.
 *
 * Produces: <project>/results/<OUTPUT_LABEL>__<imageTag>.geojson
 */

// ==================== PARAMETERS ====================
// Class names EXACTLY as they appear in the project.
// The supplied example file uses 'Islet' (singular) and 'IsletExpanded'.
// Change 'Islet' -> 'Islets' here if your project uses the plural.
def EXPORT_CLASSES = ['Islet', 'IsletExpanded'] as Set
String OUTPUT_LABEL = 'Islets'   // filename prefix
// ====================================================

def server = getCurrentServer()
if (server == null) { println "ERROR: no image open."; return }

// Gather by class directly - no need to pre-select anything in the viewer.
def matching = getAnnotationObjects().findAll {
    it.getROI() != null && (it.getPathClass()?.toString() in EXPORT_CLASSES)
}
if (matching.isEmpty()) {
    println "No annotations found in classes ${EXPORT_CLASSES}."
    return
}

def tally = matching.groupBy { it.getPathClass()?.toString() }
                    .collectEntries { k, v -> [k, v.size()] }
println "[export] Matched ${matching.size()} annotations: ${tally}"

// Build output path under <project>/results/
def baseDir = buildFilePath(PROJECT_BASE_DIR, 'results')
mkdirs(baseDir)
def imageTag = (getProjectEntry()?.getImageName() ?: 'image')
        .replaceAll('[^A-Za-z0-9._-]', '_')
def outPath = buildFilePath(baseDir, "${OUTPUT_LABEL}__${imageTag}.geojson")

// Export exactly like the GUI dialog:
//   "FEATURE_COLLECTION"   -> {"type":"FeatureCollection","features":[...]}
//   (no "EXCLUDE_MEASUREMENTS") -> measurements are included
//   (no "GZIP")                -> uncompressed .geojson
//   (no "PRETTY_JSON")         -> compact JSON, matching your example file
long t0 = System.currentTimeMillis()
exportObjectsToGeoJson(matching, outPath, "FEATURE_COLLECTION")

def sizeMB = new File(outPath).length() / (1024.0 * 1024.0)
println String.format("[export] Wrote %s (%.1f MB) in %.1fs",
        outPath, sizeMB, (System.currentTimeMillis() - t0) / 1000.0)
println "[export] Annotations left untouched in the hierarchy. Done."