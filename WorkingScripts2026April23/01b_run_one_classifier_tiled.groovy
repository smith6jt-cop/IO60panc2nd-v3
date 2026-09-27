/*
 * 01b_run_one_classifier_tiled_safe.groovy
 *
 * Safe tile-based classifier runner. Compared to the earlier version that
 * caused data loss, this script:
 *
 *   1. NEVER modifies the PathClass, ROI, name, or hierarchy membership of
 *      any existing annotation. All preserved-class annotations are
 *      untouchable.
 *
 *   2. Uses a unique throwaway PathClass '__TEMP_TILE_RUN__' for the
 *      per-tile region object the classifier is pointed at. This class is
 *      only used by this script.
 *
 *   3. Snapshots the per-class annotation count BEFORE the run. After
 *      EACH tile, re-counts. If any preserved class decreased, aborts
 *      immediately with a clear diff. Prints the final counts at end.
 *
 *   4. Uses a TEST_MODE=true default. In test mode, only ONE tile runs,
 *      and the output stays in the hierarchy for visual inspection —
 *      no GeoJSON written, nothing removed. Verify outputs look right and
 *      preserved counts are unchanged before rerunning with TEST_MODE=false.
 *
 *   5. In full run mode (TEST_MODE=false), writes per-tile objects
 *      incrementally to a JSON Lines scratch file, then on completion
 *      consolidates into ONE final GeoJSON. Crash-safe and produces a
 *      single file for downstream work. Resumable: tiles whose entries
 *      are already in the scratch file are skipped.
 *
 * Before running: run 00_backup_and_export.groovy first to back up
 * data.qpdata and export preserved-class annotations to GeoJSON.
 */

evaluate(new File(buildFilePath(PROJECT_BASE_DIR, 'scripts', 'lib', 'common.groovy')))

import qupath.lib.io.GsonTools
import qupath.lib.objects.PathObjects
import qupath.lib.objects.classes.PathClass
import qupath.lib.roi.GeometryTools
import qupath.lib.roi.ROIs
import org.locationtech.jts.geom.Polygon
import org.locationtech.jts.simplify.TopologyPreservingSimplifier

// ==================== PARAMETERS ====================
String MARKER          = 'Podoplanin'        // 'Caveolin' / 'CD34' / 'B3Tubb' / 'Podoplanin'
double TILE_SIZE_UM    = 4000.0
double MIN_AREA_UM2    = 10.0
double SIMPLIFY_TOL_PX = 1.0
boolean TEST_MODE      = false            // true = run only ONE tile then stop
int TEST_TILE_ROW      = 1               // which tile to run in TEST_MODE
int TEST_TILE_COL      = 1
// ====================================================

// Marker -> (classifier filename, output class)
def MARKER_TABLE = [
    'Caveolin'  : [classifier:'Cavdil1200',                rawClass:'Vessel_Cav_raw'],
    'CD34'      : [classifier:'CD34dil1430',             rawClass:'Vessel_CD34_raw'],
    'B3Tubb'    : [classifier:'B3Tub_Ero0.5_Thr600',     rawClass:CLASS_NERVE],
    'Podoplanin': [classifier:'PDPNmorphClos0-5_Thr300', rawClass:CLASS_LYMPHATIC]
]
def entry = MARKER_TABLE[MARKER]
if (entry == null) {
    println "ERROR: unknown MARKER '${MARKER}'. Must be one of: ${MARKER_TABLE.keySet()}"; return
}
String CLASSIFIER_NAME = entry.classifier
String RAW_CLASS       = entry.rawClass

// Throwaway class used ONLY by this script for the per-tile region object.
// Nothing else in the hierarchy should share this name.
String TEMP_CLASS = '__TEMP_TILE_RUN__'

// Classes whose count MUST NOT decrease during this run. If any of these
// lose annotations, the script aborts.
def PRESERVED_CLASSES = [
    CLASS_TISSUE, CLASS_ISLET, CLASS_ISLET_EXPANDED,
    CLASS_NERVE, CLASS_LYMPHATIC,
    'Vessel_Cav_CD34', 'Vessel_CD34_only', 'Vessel_Cav_only',
    'Vessel_Cav_raw', 'Vessel_CD34_raw', 'Blood_vessel'
] as Set

// The output class we'll be adding to. Not preserved — will grow during run.
// (If MARKER is B3Tubb/Podoplanin, RAW_CLASS is CLASS_NERVE/CLASS_LYMPHATIC
// which ARE in PRESERVED_CLASSES. Remove RAW_CLASS from preserved set so
// the count is allowed to grow.)
PRESERVED_CLASSES.remove(RAW_CLASS)

// ==================== SETUP ====================
def server = getCurrentServer()
if (server == null) { println "ERROR: no image open."; return }
def cal = server.getPixelCalibration()
double sx = (cal?.getPixelWidthMicrons())  ?: (cal?.getAveragedPixelSizeMicrons()) ?: 1.0
double sy = (cal?.getPixelHeightMicrons()) ?: (cal?.getAveragedPixelSizeMicrons()) ?: 1.0
double areaScale = sx * sy

def classifierPath = buildFilePath(PROJECT_BASE_DIR, 'classifiers', 'pixel_classifiers', "${CLASSIFIER_NAME}.json")
if (!new File(classifierPath).exists()) {
    println "ERROR: classifier '${CLASSIFIER_NAME}.json' not found at ${classifierPath}"; return
}

def tissueAnns = getAnnotationObjects().findAll {
    it.getROI() && it.getPathClass()?.toString() == CLASS_TISSUE
}
if (tissueAnns.isEmpty()) { println "ERROR: no '${CLASS_TISSUE}' annotations."; return }

// Read Tissue geometry ONCE into a local variable. After this point the
// script does not depend on Tissue hierarchy membership.
def tissueGeoms = tissueAnns.collect { it.getROI().getGeometry() }
def tissueUnion = (tissueGeoms.size() == 1) ? tissueGeoms[0] : safeUnionPolys(tissueGeoms)
def plane = tissueAnns[0].getROI().getImagePlane()
def env = tissueUnion.getEnvelopeInternal()

// ==================== PRE-RUN HIERARCHY SNAPSHOT ====================
def snapshotCounts = { ->
    def m = [:]
    getAnnotationObjects().each { ann ->
        def c = ann.getPathClass()?.toString() ?: '(unclassified)'
        m[c] = (m[c] ?: 0) + 1
    }
    m
}

def preCounts = snapshotCounts()
println "[safe] Pre-run annotation counts:"
preCounts.each { cls, n -> println "[safe]   ${cls}: ${n}" }
println "[safe] Preserved classes (must not decrease): ${PRESERVED_CLASSES}"

def verifyPreserved = { String context ->
    def now = snapshotCounts()
    def violations = []
    PRESERVED_CLASSES.each { cls ->
        int before = preCounts[cls] ?: 0
        int after = now[cls] ?: 0
        if (after < before) {
            violations << "${cls}: ${before} -> ${after} (LOST ${before - after})"
        }
    }
    if (!violations.isEmpty()) {
        println "[safe] ABORT: preserved class count decreased at ${context}"
        violations.each { println "[safe]   ${it}" }
        return false
    }
    return true
}

// ==================== TILING SETUP ====================
double tileSizePxX = TILE_SIZE_UM / sx
double tileSizePxY = TILE_SIZE_UM / sy
int nCols = (int) Math.ceil((env.getMaxX() - env.getMinX()) / tileSizePxX)
int nRows = (int) Math.ceil((env.getMaxY() - env.getMinY()) / tileSizePxY)
println "[safe] classifier=${CLASSIFIER_NAME} -> ${RAW_CLASS}"
println "[safe] Tile size ${TILE_SIZE_UM}um = ${tileSizePxX as int}x${tileSizePxY as int} px; grid ${nCols}x${nRows} = ${nCols * nRows} tiles"
if (TEST_MODE) println "[safe] TEST_MODE=true: will run only tile r${TEST_TILE_ROW}_c${TEST_TILE_COL} then stop"

// Output paths: single final GeoJSON + JSON Lines scratch file for crash safety.
def baseDir = buildFilePath(PROJECT_BASE_DIR, 'results')
mkdirs(baseDir)
def imageTag = (getProjectEntry()?.getImageName() ?: server?.getShortServerName() ?: 'img')
        .replaceAll('[^A-Za-z0-9._-]', '_')
def safeClass = RAW_CLASS.replaceAll('[^A-Za-z0-9._-]', '_')
def finalPath   = buildFilePath(baseDir, "${safeClass}__${imageTag}.geojson")
def scratchPath = buildFilePath(baseDir, "${safeClass}__${imageTag}.scratch.jsonl")
println "[safe] Final output: ${finalPath}"
println "[safe] Scratch (crash-safe incremental): ${scratchPath}"

// Build set of (row, col) pairs already present in scratch (for resume)
def completedTiles = [] as Set
def scratchFile = new File(scratchPath)
if (!TEST_MODE && scratchFile.exists()) {
    scratchFile.eachLine('UTF-8') { line ->
        def m = (line =~ /^\{"row":\s*(\d+),\s*"col":\s*(\d+)/)
        if (m.find()) completedTiles << "${m.group(1)}_${m.group(2)}"
    }
    println "[safe] Resume: ${completedTiles.size()} tiles already in scratch"
}

def rawPC  = PathClass.fromString(RAW_CLASS)
def tempPC = PathClass.fromString(TEMP_CLASS)

// Safety: fail fast if any existing annotation has our throwaway class
def existingTemp = getAnnotationObjects().findAll { it.getPathClass()?.toString() == TEMP_CLASS }
if (!existingTemp.isEmpty()) {
    println "[safe] ABORT: ${existingTemp.size()} annotations already have class '${TEMP_CLASS}'. Remove them before running."
    return
}

// ==================== TILE LOOP ====================
long tStart = System.currentTimeMillis()
int tilesProcessed = 0, tilesSkipped = 0, tilesEmpty = 0, tilesFailed = 0
int totalObjects = 0

def runOneTile = { int row, int col ->
    String tileKey = "${row}_${col}"
    if (completedTiles.contains(tileKey)) {
        tilesSkipped++
        return 'skip'
    }

    double tx = env.getMinX() + col * tileSizePxX
    double ty = env.getMinY() + row * tileSizePxY
    double tw = Math.min(tileSizePxX, env.getMaxX() - tx)
    double th = Math.min(tileSizePxY, env.getMaxY() - ty)
    if (tw <= 0 || th <= 0) return 'skip'

    def tileROI = ROIs.createRectangleROI(tx, ty, tw, th, plane)
    def tileGeom = tileROI.getGeometry()

    def overlap = null
    try { overlap = tileGeom.intersection(tissueUnion) } catch (Throwable t) { overlap = null }
    if (overlap == null || overlap.isEmpty() || overlap.getArea() < 1.0) {
        // Record as empty in scratch so we don't revisit on resume
        if (!TEST_MODE) {
            scratchFile << /{"row":${row},"col":${col},"empty":true}/ + "\n"
        }
        tilesEmpty++
        return 'empty'
    }

    // Temp tile annotation with TEMP_CLASS.
    def overlapROI = GeometryTools.geometryToROI(overlap, plane)
    def tempAnn = PathObjects.createAnnotationObject(overlapROI, tempPC)
    tempAnn.setName("__temp_tile_${row}_${col}")
    addObject(tempAnn)
    fireHierarchyUpdate()

    long t0 = System.currentTimeMillis()
    int tileObjectCount = 0
    boolean tileFailed = false
    String tileError = null
    try {
        removeSelectedObjects()
        getCurrentHierarchy().getSelectionModel().setSelectedObject(tempAnn)

        def before = new HashSet(getAnnotationObjects())
        def newObjs = null
        int attempts = 0
        while (attempts < 2 && newObjs == null) {
            attempts++
            try {
                createAnnotationsFromPixelClassifier(CLASSIFIER_NAME, MIN_AREA_UM2, 0.0)
                newObjs = getAnnotationObjects().findAll {
                    !before.contains(it) && !it.is(tempAnn)
                }
            } catch (Throwable e) {
                tileError = e.getMessage()?.take(200)
                println "[safe]   tile r${row}_c${col} attempt ${attempts} FAILED: ${tileError}"
                def partial = getAnnotationObjects().findAll {
                    !before.contains(it) && !it.is(tempAnn)
                }
                if (!partial.isEmpty()) removeObjects(partial, true)
                if (attempts < 2) { System.gc(); Thread.sleep(500) }
                else tileFailed = true
            }
        }

        if (tileFailed) {
            if (!TEST_MODE) {
                scratchFile << /{"row":${row},"col":${col},"failed":true}/ + "\n"
            }
        } else if (newObjs == null || newObjs.isEmpty()) {
            if (!TEST_MODE) {
                scratchFile << /{"row":${row},"col":${col},"empty":true}/ + "\n"
            }
            tilesEmpty++
        } else {
            newObjs.each { it.setPathClass(rawPC) }

            // Split multipolygons
            def splitAdd = []
            def splitRemove = []
            newObjs.each { ann ->
                def roi = ann.getROI(); if (roi == null) return
                def g = roi.getGeometry(); if (g == null || g.isEmpty()) return
                int nParts = g.getNumGeometries()
                if (nParts > 1) {
                    splitRemove << ann
                    for (int i = 0; i < nParts; i++) {
                        def part = g.getGeometryN(i)
                        if (!(part instanceof Polygon) || part.isEmpty()) continue
                        double areaUm2 = part.getArea() * areaScale
                        if (areaUm2 < MIN_AREA_UM2) continue
                        def partROI = GeometryTools.geometryToROI(part, roi.getImagePlane())
                        splitAdd << PathObjects.createAnnotationObject(partROI, rawPC)
                    }
                }
            }
            if (!splitRemove.isEmpty()) {
                removeObjects(splitRemove, true)
                addObjects(splitAdd)
                newObjs = (newObjs - splitRemove) + splitAdd
            }

            // Simplify
            newObjs.each { ann ->
                def roi = ann.getROI(); if (roi == null) return
                def g = roi.getGeometry(); if (g == null || g.isEmpty()) return
                def s = TopologyPreservingSimplifier.simplify(g, SIMPLIFY_TOL_PX)
                if (s != null && !s.isEmpty()) {
                    ann.setROI(GeometryTools.geometryToROI(s, roi.getImagePlane()))
                }
            }

            tileObjectCount = newObjs.size()
            totalObjects += tileObjectCount

            if (TEST_MODE) {
                // Leave newObjs IN the hierarchy for visual inspection.
                // Do not remove, do not serialize.
            } else {
                // Serialize this tile's objects as one JSON Lines record.
                // Format: {"row":R,"col":C,"features":[<geojson feature>,...]}
                def featuresJson = GsonTools.getInstance(false).toJson(newObjs)
                scratchFile << /{"row":${row},"col":${col},"features":/ + featuresJson + "}\n"
                // Remove from hierarchy — they are persisted in scratch
                removeObjects(newObjs, true)
            }
        }
    } finally {
        // ALWAYS remove the temp tile annotation, no matter what.
        removeObjectAndDescendants(tempAnn)
        fireHierarchyUpdate()
        System.gc()
    }

    double sec = (System.currentTimeMillis() - t0) / 1000.0
    if (tileFailed) { tilesFailed++; return [status:'fail', sec:sec, n:0] }
    else { tilesProcessed++; return [status:'ok', sec:sec, n:tileObjectCount] }
}

def processOrder
if (TEST_MODE) {
    processOrder = [[TEST_TILE_ROW, TEST_TILE_COL]]
} else {
    processOrder = []
    for (int r = 0; r < nRows; r++) for (int c = 0; c < nCols; c++) processOrder << [r, c]
}

def aborted = false
for (def rc : processOrder) {
    int row = rc[0], col = rc[1]
    def result = runOneTile(row, col)
    if (result == 'skip' || result == 'empty') continue

    // Verify preserved classes after each tile
    if (!verifyPreserved("after tile r${row}_c${col}")) {
        aborted = true
        break
    }

    if (result instanceof Map) {
        double elapsed = (System.currentTimeMillis() - tStart) / 1000.0
        int done = tilesProcessed + tilesSkipped + tilesEmpty + tilesFailed
        int remaining = processOrder.size() - done
        double etaMin = (done > 0 && remaining > 0) ? (elapsed / done * remaining / 60.0) : 0.0
        def tag = (result.status == 'fail') ? "FAILED" : "${result.n} objects"
        println String.format("[safe] tile %d/%d (r%d c%d): %s, %.1fs  | done=%d remain=%d failed=%d ETA=%.1fmin",
                done, processOrder.size(), row, col, tag, result.sec, done, remaining, tilesFailed, etaMin)
    }
}

// ==================== CONSOLIDATE SCRATCH TO FINAL GEOJSON ====================
// In full-run mode, convert the JSON Lines scratch file into one FeatureCollection.
if (!TEST_MODE && !aborted) {
    println "[safe] Consolidating scratch file to final GeoJSON..."
    long tCons = System.currentTimeMillis()
    def featuresArrays = []
    int featureCount = 0
    scratchFile.eachLine('UTF-8') { line ->
        // Extract the "features":[...] array portion of each line, if present.
        // Lines with "empty":true or "failed":true have no features.
        def fmIdx = line.indexOf('"features":')
        if (fmIdx < 0) return
        // The remainder from "features": to the last '}' is the array + closing brace.
        def fragment = line.substring(fmIdx + '"features":'.length())
        // Strip trailing '}' that closes the outer tile record.
        int lastBrace = fragment.lastIndexOf('}')
        if (lastBrace < 0) return
        def arr = fragment.substring(0, lastBrace).trim()
        if (arr.startsWith('[') && arr.endsWith(']') && arr.length() > 2) {
            featuresArrays << arr.substring(1, arr.length() - 1)  // drop [ and ]
            // Rough feature count: count top-level '{' entries. Good enough for log.
        }
    }
    // Build the final FeatureCollection by joining all feature arrays with commas.
    def joined = featuresArrays.findAll { it }.join(',')
    new File(finalPath).withWriter('UTF-8') { w ->
        w.write('{"type":"FeatureCollection","features":[')
        w.write(joined)
        w.write(']}')
    }
    double mb = new File(finalPath).length() / (1024.0 * 1024.0)
    double sec = (System.currentTimeMillis() - tCons) / 1000.0
    println String.format("[safe] Wrote final GeoJSON: %s (%.1f MB) in %.1fs",
            finalPath, mb, sec)
    println "[safe] Scratch file retained at ${scratchPath} (delete if no longer needed)."
}

// ==================== POST-RUN REPORT ====================
def postCounts = snapshotCounts()
println "[safe] Post-run annotation counts:"
postCounts.each { cls, n -> println "[safe]   ${cls}: ${n}" }

println "[safe] Diff vs pre-run:"
def allClasses = (preCounts.keySet() + postCounts.keySet()) as Set
allClasses.sort().each { cls ->
    int before = preCounts[cls] ?: 0
    int after = postCounts[cls] ?: 0
    if (before != after) {
        println "[safe]   ${cls}: ${before} -> ${after}  (${after - before >= 0 ? '+' : ''}${after - before})"
    }
}

if (aborted) {
    println "[safe] ABORTED due to preserved-class violation. Inspect hierarchy."
} else if (TEST_MODE) {
    println "[safe] TEST_MODE complete."
    println "[safe] Objects are IN the hierarchy for inspection (class: ${RAW_CLASS})."
    println "[safe] Once you've verified, delete them (or use resetDetectionClassifications on annotations),"
    println "[safe] set TEST_MODE=false, and rerun to process the full grid."
} else {
    println String.format("[safe] Done. processed=%d skipped=%d empty=%d failed=%d totalObjects=%d in %.1fmin",
            tilesProcessed, tilesSkipped, tilesEmpty, tilesFailed, totalObjects,
            (System.currentTimeMillis() - tStart) / 60000.0)
    println "[safe] Final GeoJSON: ${finalPath}"
}