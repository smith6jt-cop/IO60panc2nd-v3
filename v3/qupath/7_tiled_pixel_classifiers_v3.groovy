/*
 * 7_tiled_pixel_classifiers_v3.groovy   (v3 copy of scripts/7_tiled_pixel_classifiers.groovy; the original is unchanged)
 *
 * v3 changes:
 *   - the tiled area is the union of every Tissue AND every pLN annotation (REGION_CLASSES), so a lymph
 *     node drawn disjoint from the Tissue polygon (6533, HDL112, HDL115, 6551) is classified too; pLN is
 *     also added to the preserved classes the per-tile safety check counts. (A pLN inside Tissue adds
 *     nothing: the union is the same area.)
 *   - 2026-10-01, for run_v3.sh step '7': args[0] names the output root. Files go to
 *     <root>/<Class>/<Class>__<image tag>.geojson -- one file per image and class, the layout the Senior
 *     typing lane reads (context.load_objects) -- and the resume files to <root>/scratch/. Without
 *     args the files go to results/ exactly as before.
 *   - a tile that FAILED is no longer counted as done on resume (it was: the resume reader matched every
 *     row/col line), so a re-run retries it, and the outcome goes to V3_RESULT: the step passes only when
 *     every marker wrote its file and no tile failed. Nothing is saved to the data file.
 * The classifiers, the tile size, the minimum area and the simplification are the original's.
 *
 * Original header:
 * 07_tiled_pixel_classifiers.groovy
 *
 * Workflow step 9. Run the three structural pixel classifiers (B3TUBB for
 * Nerve, CD34 for Blood_vessel, Podoplanin for Lymphatic) tile by tile across
 * each whole-slide image. Outputs are written as per-marker GeoJSON under
 * results/ and are NOT retained in the QuPath hierarchy - the CD34 object
 * coverage in particular is too large and complex to keep loaded in the
 * project without overwhelming system resources.
 *
 * Caveolin was removed from this stage per current workflow. If Caveolin
 * processing is reinstated later, add an entry back to MARKER_TABLE.
 *
 * GeoJSON class labels inside each file:
 *   Nerve      (from the B3TUBB classifier)
 *   Blood_vessel (from the CD34 classifier)
 *   Lymphatic  (from the Podoplanin classifier)
 *
 * Safety model (unchanged from the original 01b_run_one_classifier_tiled):
 *   - a throwaway '__TEMP_TILE_RUN__' class is used for the per-tile ROI
 *   - preserved-class counts are snapshotted before the run and verified
 *     after every tile; a decrease aborts
 *   - per-marker scratch JSONL file permits crash-safe resume
 *
 * Run via: Run -> Run for project.
 */

evaluate(new File(buildFilePath(PROJECT_BASE_DIR, 'scripts', 'lib', 'common.groovy')))

import qupath.lib.io.GsonTools
import qupath.lib.io.PathIO
import qupath.lib.objects.PathObjects
import qupath.lib.objects.classes.PathClass
import qupath.lib.roi.GeometryTools
import qupath.lib.roi.ROIs
import org.locationtech.jts.geom.Polygon
import org.locationtech.jts.simplify.TopologyPreservingSimplifier

// ==================== PARAMETERS ====================
List<String> REGION_CLASSES = [CLASS_TISSUE, 'pLN']   // v3: tiles cover the union of these annotations
double TILE_SIZE_UM    = 2000.0
double MIN_AREA_UM2    = 10.0
double SIMPLIFY_TOL_PX = 1.0
boolean TEST_MODE      = false
int TEST_TILE_ROW      = 1
int TEST_TILE_COL      = 1

// Order matters only for logging; each marker runs to completion before the next.
def MARKERS = ['B3Tubb', 'CD34', 'Podoplanin']

// Per-marker config: classifier name, output PathClass, whether to retain in hierarchy.
// keepInHierarchy is always false in the current workflow; see header.
def MARKER_TABLE = [
    'CD34'      : [classifier:'CD34dil1430',             rawClass:CLASS_BLOOD_VESSEL, keepInHierarchy:false],
    'B3Tubb'    : [classifier:'B3Tub_Ero0.5_Thr600',     rawClass:CLASS_NERVE,        keepInHierarchy:false],
    'Podoplanin': [classifier:'PDPNmorphClos0-5_Thr300', rawClass:CLASS_LYMPHATIC,    keepInHierarchy:false]
]
// ====================================================

String TEMP_CLASS = '__TEMP_TILE_RUN__'

def server = getCurrentServer()
if (server == null) { println "[07] ERROR: no image open."; return }
def cal = server.getPixelCalibration()
double sx = (cal?.getPixelWidthMicrons())  ?: (cal?.getAveragedPixelSizeMicrons()) ?: 1.0
double sy = (cal?.getPixelHeightMicrons()) ?: (cal?.getAveragedPixelSizeMicrons()) ?: 1.0
double areaScale = sx * sy

def tissueAnns = getAnnotationObjects().findAll {
    it.getROI() && it.getPathClass()?.toString() in REGION_CLASSES
}
if (tissueAnns.isEmpty()) { println "[07] ERROR: no ${REGION_CLASSES} annotations."; return }
println "[07] tile area = union of ${tissueAnns.countBy { it.getPathClass().toString() }}"

def tissueGeoms = tissueAnns.collect { it.getROI().getGeometry() }
def tissueUnion = (tissueGeoms.size() == 1) ? tissueGeoms[0] : safeUnionPolys(tissueGeoms)
def plane = tissueAnns[0].getROI().getImagePlane()
def env = tissueUnion.getEnvelopeInternal()

double tileSizePxX = TILE_SIZE_UM / sx
double tileSizePxY = TILE_SIZE_UM / sy
int nCols = (int) Math.ceil((env.getMaxX() - env.getMinX()) / tileSizePxX)
int nRows = (int) Math.ceil((env.getMaxY() - env.getMinY()) / tileSizePxY)

// v3: args[0] = output root (run_v3.sh passes {OUT}/structures); the tag then follows the driver's convention
boolean V3_LAYOUT = binding.hasVariable('args') && args && args.length > 0 && args[0]?.trim()
def baseDir = V3_LAYOUT ? args[0].trim() : buildFilePath(PROJECT_BASE_DIR, 'results')
mkdirs(baseDir)
def imageTag = (getProjectEntry()?.getImageName() ?: server?.getShortServerName() ?: 'img')
        .replaceAll(V3_LAYOUT ? '[^A-Za-z0-9._-]+' : '[^A-Za-z0-9._-]', '_')
def markerResults = [:]     // marker -> what happened, for V3_RESULT

def tempPC = PathClass.fromString(TEMP_CLASS)

// ==================== PER-MARKER RUNNER ====================
def runMarker = { String marker ->
    def entry = MARKER_TABLE[marker]
    markerResults[marker] = [ok: false, message: 'did not finish']
    if (entry == null) { println "[07] Skipping unknown marker '${marker}'"; markerResults[marker].message = 'unknown marker'; return }
    String CLASSIFIER_NAME = entry.classifier
    String RAW_CLASS       = entry.rawClass
    boolean KEEP_IN_HIER   = entry.keepInHierarchy

    println "\n[07] ===== ${marker} -> ${CLASSIFIER_NAME} -> class '${RAW_CLASS}' (keep=${KEEP_IN_HIER}) ====="

    def classifierPath = buildFilePath(PROJECT_BASE_DIR, 'classifiers', 'pixel_classifiers', "${CLASSIFIER_NAME}.json")
    if (!new File(classifierPath).exists()) {
        println "[07] ERROR: classifier '${CLASSIFIER_NAME}.json' not found at ${classifierPath}"
        markerResults[marker].message = "classifier ${CLASSIFIER_NAME}.json not found".toString()
        return
    }

    def safeClass = RAW_CLASS.replaceAll('[^A-Za-z0-9._-]', '_')
    def finalDir   = V3_LAYOUT ? buildFilePath(baseDir, safeClass) : baseDir
    def scratchDir = V3_LAYOUT ? buildFilePath(baseDir, 'scratch') : baseDir
    mkdirs(finalDir); mkdirs(scratchDir)
    def finalPath   = buildFilePath(finalDir, "${safeClass}__${imageTag}.geojson")
    def scratchPath = buildFilePath(scratchDir, "${safeClass}__${imageTag}.scratch.jsonl")
    println "[07]   final:  ${finalPath}"
    println "[07]   scratch:${scratchPath}"

    // Preserved-class snapshot. RAW_CLASS may coincide with a preserved class
    // (Blood_vessel, Nerve, Lymphatic); exclude it so the count is allowed
    // to grow during this marker's run.
    def PRESERVED_CLASSES = [
        CLASS_TISSUE, 'pLN', CLASS_ISLET, CLASS_ISLET_EXPANDED,
        CLASS_NERVE, CLASS_LYMPHATIC, CLASS_BLOOD_VESSEL
    ] as Set
    PRESERVED_CLASSES.remove(RAW_CLASS)

    def snapshotCounts = { ->
        def m = [:]
        getAnnotationObjects().each { ann ->
            def c = ann.getPathClass()?.toString() ?: '(unclassified)'
            m[c] = (m[c] ?: 0) + 1
        }
        m
    }
    def preCounts = snapshotCounts()
    def verifyPreserved = { String context ->
        def now = snapshotCounts()
        def violations = []
        PRESERVED_CLASSES.each { cls ->
            int before = preCounts[cls] ?: 0
            int after = now[cls] ?: 0
            if (after < before) violations << "${cls}: ${before} -> ${after}"
        }
        if (!violations.isEmpty()) {
            println "[07]   ABORT at ${context}:"
            violations.each { println "[07]     ${it}" }
            return false
        }
        return true
    }

    // Resume support
    def completedTiles = [] as Set
    def scratchFile = new File(scratchPath)
    if (!TEST_MODE && scratchFile.exists()) {
        scratchFile.eachLine('UTF-8') { line ->
            if (line.contains('"failed":true')) return      // v3: a failed tile is retried, not counted as done
            def m = (line =~ /^\{"row":\s*(\d+),\s*"col":\s*(\d+)/)
            if (m.find()) completedTiles << "${m.group(1)}_${m.group(2)}"
        }
        println "[07]   Resume: ${completedTiles.size()} tiles already in scratch"
    }

    // Guard against orphaned temp annotations
    def existingTemp = getAnnotationObjects().findAll { it.getPathClass()?.toString() == TEMP_CLASS }
    if (!existingTemp.isEmpty()) {
        println "[07]   Cleaning ${existingTemp.size()} orphaned ${TEMP_CLASS} annotations"
        removeObjects(existingTemp, true)
    }

    def rawPC = PathClass.fromString(RAW_CLASS)

    long tStart = System.currentTimeMillis()
    int tilesProcessed = 0, tilesSkipped = 0, tilesEmpty = 0, tilesFailed = 0
    int totalObjects = 0

    def runOneTile = { int row, int col ->
        String tileKey = "${row}_${col}"
        if (completedTiles.contains(tileKey)) { tilesSkipped++; return 'skip' }

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
            if (!TEST_MODE) scratchFile << /{"row":${row},"col":${col},"empty":true}/ + "\n"
            tilesEmpty++; return 'empty'
        }

        def overlapROI = GeometryTools.geometryToROI(overlap, plane)
        def tempAnn = PathObjects.createAnnotationObject(overlapROI, tempPC)
        tempAnn.setName("__temp_tile_${row}_${col}")
        addObject(tempAnn)
        fireHierarchyUpdate()

        long t0 = System.currentTimeMillis()
        int tileObjectCount = 0
        boolean tileFailed = false
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
                    def partial = getAnnotationObjects().findAll {
                        !before.contains(it) && !it.is(tempAnn)
                    }
                    if (!partial.isEmpty()) removeObjects(partial, true)
                    if (attempts < 2) { System.gc(); Thread.sleep(500) } else tileFailed = true
                }
            }

            if (tileFailed) {
                if (!TEST_MODE) scratchFile << /{"row":${row},"col":${col},"failed":true}/ + "\n"
            } else if (newObjs == null || newObjs.isEmpty()) {
                if (!TEST_MODE) scratchFile << /{"row":${row},"col":${col},"empty":true}/ + "\n"
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
                    // Leave newObjs in hierarchy for inspection.
                } else {
                    // Always serialize to scratch, then remove from hierarchy.
                    // The keepInHierarchy flag guards the re-import step below,
                    // not the tile-loop removal, so that scratch is always the
                    // single source of truth for resume.
                    def featuresJson = GsonTools.getInstance(false).toJson(newObjs)
                    scratchFile << /{"row":${row},"col":${col},"features":/ + featuresJson + "}\n"
                    removeObjects(newObjs, true)
                }
            }
        } finally {
            removeObjectAndDescendants(tempAnn)
            fireHierarchyUpdate()
//            System.gc()
        }

        double sec = (System.currentTimeMillis() - t0) / 1000.0
        if (tileFailed) { tilesFailed++; return [status:'fail', sec:sec, n:0] }
        else { tilesProcessed++; return [status:'ok', sec:sec, n:tileObjectCount] }
    }

    def processOrder = []
    if (TEST_MODE) processOrder = [[TEST_TILE_ROW, TEST_TILE_COL]]
    else for (int r = 0; r < nRows; r++) for (int c = 0; c < nCols; c++) processOrder << [r, c]

    boolean aborted = false
    for (def rc : processOrder) {
        int row = rc[0], col = rc[1]
        def result = runOneTile(row, col)
        if (result == 'skip' || result == 'empty') continue
        if (!verifyPreserved("after tile r${row}_c${col}")) { aborted = true; break }
        if (result instanceof Map) {
            double elapsed = (System.currentTimeMillis() - tStart) / 1000.0
            int done = tilesProcessed + tilesSkipped + tilesEmpty + tilesFailed
            int remaining = processOrder.size() - done
            double etaMin = (done > 0 && remaining > 0) ? (elapsed / done * remaining / 60.0) : 0.0
            def tag = (result.status == 'fail') ? "FAILED" : "${result.n} objects"
            println String.format("[07]   tile %d/%d (r%d c%d): %s, %.1fs  ETA=%.1fmin",
                    done, processOrder.size(), row, col, tag, result.sec, etaMin)
        }
    }

    // ==================== CONSOLIDATE ====================
    if (!TEST_MODE && !aborted) {
        println "[07]   Consolidating scratch to final GeoJSON..."
        def featuresArrays = []
        scratchFile.eachLine('UTF-8') { line ->
            def fmIdx = line.indexOf('"features":'); if (fmIdx < 0) return
            def fragment = line.substring(fmIdx + '"features":'.length())
            int lastBrace = fragment.lastIndexOf('}'); if (lastBrace < 0) return
            def arr = fragment.substring(0, lastBrace).trim()
            if (arr.startsWith('[') && arr.endsWith(']') && arr.length() > 2) {
                featuresArrays << arr.substring(1, arr.length() - 1)
            }
        }
        def joined = featuresArrays.findAll { it }.join(',')
        new File(finalPath).withWriter('UTF-8') { w ->
            w.write('{"type":"FeatureCollection","features":[')
            w.write(joined)
            w.write(']}')
        }
        double mb = new File(finalPath).length() / (1024.0 * 1024.0)
        println String.format("[07]   Wrote %s (%.1f MB)", finalPath, mb)
        markerResults[marker] = [ok: tilesFailed == 0, output: finalPath, classifier: CLASSIFIER_NAME, raw_class: RAW_CLASS,
                                 tiles_processed: tilesProcessed, tiles_resumed: tilesSkipped, tiles_empty: tilesEmpty,
                                 tiles_failed: tilesFailed, objects_this_run: totalObjects, megabytes: mb,
                                 message: (tilesFailed == 0 ? 'complete' : "${tilesFailed} tile(s) failed".toString())]

        // For CD34 only: re-import as Blood_vessel annotations
        if (KEEP_IN_HIER) {
            try {
                def imported = PathIO.readObjectsFromGeoJson(new File(finalPath).toPath())
                def anns = imported.findAll { it.isAnnotation() || it.hasROI() }
                anns.each { it.setPathClass(rawPC) }
                addObjects(anns)
                fireHierarchyUpdate()
                println "[07]   Re-imported ${anns.size()} '${RAW_CLASS}' annotations into hierarchy"
            } catch (Throwable t) {
                println "[07]   WARN: re-import failed (${t.class.simpleName}: ${t.message})"
                println "[07]         Check QuPath version; PathIO.readObjectsFromGeoJson requires >= 0.4."
            }
        }
    }

    if (aborted) markerResults[marker].message = 'aborted: a preserved annotation class lost objects'
    println String.format("[07]   %s done. processed=%d skipped=%d empty=%d failed=%d objects=%d in %.1fmin",
            marker, tilesProcessed, tilesSkipped, tilesEmpty, tilesFailed, totalObjects,
            (System.currentTimeMillis() - tStart) / 60000.0)
}

// ==================== RUN ALL MARKERS ====================
long runStart = System.currentTimeMillis()
MARKERS.each { marker -> runMarker(marker) }
println String.format("[07] All markers done in %.1f min",
        (System.currentTimeMillis() - runStart) / 60000.0)

// v3: the verdict for the run_v3 footer. Nothing was saved: every object went to GeoJSON and the temp tiles were removed.
boolean allOk = !TEST_MODE && MARKERS.every { markerResults[it]?.ok == true }
def leftovers = getAnnotationObjects().count { it.getPathClass()?.toString() == TEMP_CLASS }
V3_RESULT = [ok: allOk && leftovers == 0, changed: false,
             message: MARKERS.collect { "${it}: ${markerResults[it]?.message}" }.join('; ').toString(),
             details: [markers: markerResults, regions: tissueAnns.countBy { it.getPathClass().toString() },
                       tile_um: TILE_SIZE_UM, min_area_um2: MIN_AREA_UM2, simplify_px: SIMPLIFY_TOL_PX,
                       temp_annotations_left: leftovers, minutes: (System.currentTimeMillis() - runStart) / 60000.0,
                       outputs: MARKERS.collect { markerResults[it]?.output }.findAll { it }]]
println "[07] ${V3_RESULT.ok ? 'OK' : 'FAILED'}: ${V3_RESULT.message}"
