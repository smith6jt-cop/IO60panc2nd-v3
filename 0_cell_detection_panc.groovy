/*
 * 0_cell_detection.groovy
 *
 * Workflow step 1: InstanSeg cell detection over the Tissue annotation ONLY.
 *
 * WHY THIS SCRIPT EXISTS
 * ----------------------
 * Detection used to be a manual GUI step — the only unscripted stage in the pipeline — and it
 * silently corrupted 8 of 20 donors. The QuPath workflow history shows the affected runs used:
 *
 *     selectAnnotations()            <-- selects Tissue AND every nested Islet_N child
 *     InstanSeg...detectObjects()
 *
 * InstanSeg's detectObjects() runs inference once per SELECTED parent object and does not
 * deduplicate across parents (PruneObjectOutputHandler dedupes only across tiles WITHIN one
 * parent). So every nucleus inside an islet was segmented twice — once under Tissue, once under
 * its islet. The two passes tile from different ROI origins, so the copies are near-identical but
 * not equal: centroids ~1 px apart, areas within a fraction of a percent, different vertex counts,
 * two distinct UUIDs. 8_hierarchy_parenting.groovy then re-parents both to the islet, which is why
 * 100% of the duplicates landed in islet cores.
 *
 * Measured damage: 155,669 duplicate pairs across donors 6450, 6505, 6523, 6534, 6566, 6591, 6623.
 * Roughly half to two-thirds of every islet cell in those donors was a duplicate. Nothing
 * downstream detected it — the UUIDs are distinct, so every uniqueness check passed.
 *
 * THE RULES THIS SCRIPT ENFORCES
 * ------------------------------
 *   1. NEVER selectAnnotations(). Select the Tissue annotation explicitly, by classification.
 *   2. Protect annotations below Tissue. QuPath 0.7.0's PixelProcessor (used by InstanSeg) calls
 *      clearChildObjects() on every parent it processes, which would detach Islet / IsletExpanded /
 *      vessel annotations. They are parked at the root during detection and restored afterwards
 *      (names and "Islet: UID" untouched). A Tissue nested inside another Tissue is refused.
 *   3. removeDetections() first, so a re-run replaces rather than accumulates.
 *   4. Verify the InstanSeg result and that no two detections share a centroid; throw (and do not
 *      save) if either check fails.
 *
 * Normally run BEFORE 1_panc_objects.groovy. Re-detecting an image that already has islets is
 * supported (rule 2); afterwards re-run 2 -> 3 -> 8 (4 is independent of detections).
 *
 * Headless:
 *   /home/smith6jt/QuPath/bin/QuPath script \
 *       --project /home/smith6jt/IO60panc2nd/project.qpproj \
 *       --image "6539_Scan1.er.qptiff - resolution #1" \
 *       /home/smith6jt/IO60panc2nd/scripts/0_cell_detection_panc.groovy
 */

import qupath.lib.images.servers.ColorTransforms

// ==================== PARAMETERS ====================
// Matches every run recorded in the project workflow history (QuPath 0.7.0).
String MODEL_PATH   = "/home/smith6jt/QuPath/v0.7/instanseg/downloaded/fluorescence_nuclei_and_cells-0.1.1"
String DEVICE       = "gpu0"      // "cpu" if no CUDA device
int    TILE_DIMS    = 512
int    PAD          = 32
int    N_THREADS    = 12
boolean MAKE_MEASUREMENTS = false // measurements come from 2_fast_cell_measurements.groovy
String TISSUE_CLASS = "Tissue"
double PIXEL_SIZE_UM = 0.5078     // same value 1_panc_objects.groovy sets; InstanSeg derives its downsample from it
List<String> INPUT_CHANNELS = ["DAPI", "CD68", "CD44", "Vimentin", "CD99", "E-cadherin", "M2Gal3", "HLA-A",
                               "Pan-Cytokeratin", "CD79a", "CD3e", "Beta-actin", "Iba1"]
double DUP_RADIUS_UM = 1.0        // post-run duplicate check; matches phenocycler.cell_qc
double DUP_AREA_TOL  = 0.15
// ====================================================

def server = getCurrentServer()
if (server == null) { println "ERROR: no image open."; return }
// Calibration must be in place BEFORE detection: InstanSeg computes its downsample from it, and the
// duplicate-check radius below converts um -> px with it.
setImageType('FLUORESCENCE')
setPixelSizeMicrons(PIXEL_SIZE_UM, PIXEL_SIZE_UM)
server = getCurrentServer()
def presentNames = server.getMetadata().getChannels().collect { it.getName() }
def missing = INPUT_CHANNELS.findAll { !(it in presentNames) }
if (!missing.isEmpty()) {
    throw new IllegalStateException("[detect] Channel(s) not found: ${missing}. Normalize channel names before detection.")
}
def imageName = getProjectEntry()?.getImageName() ?: server.getMetadata()?.getName()
println "[detect] image: ${imageName}"

// ---- 1. select ONLY the tissue annotation ---------------------------------
resetSelection()
selectObjectsByClassification(TISSUE_CLASS)
def selected = getSelectedObjects()
if (selected.isEmpty()) {
    println "ERROR: no '${TISSUE_CLASS}' annotation found. Create it before running detection."
    return
}

// ---- 2. nested Tissue is refused; other annotations below Tissue are parked --
// selectObjectsByClassification() already guarantees the selection is Tissue-only. The hazard is
// below the Tissue: PixelProcessor.processTiled() -> parent.clearChildObjects() removes EVERY child
// of each processed parent, annotations included. A Tissue nested in another Tissue would also be
// processed twice, so that is refused outright.
def selectedSet = Collections.newSetFromMap(new IdentityHashMap())
selectedSet.addAll(selected)
def hasSelectedAncestor = { obj ->
    def p = obj.getParent()
    while (p != null) { if (selectedSet.contains(p)) return true; p = p.getParent() }
    return false
}
def nestedTissue = selected.findAll { hasSelectedAncestor(it) }
if (!nestedTissue.isEmpty()) {
    println "ERROR: ${nestedTissue.size()} ${TISSUE_CLASS} annotation(s) are nested inside another ${TISSUE_CLASS}."
    println "       They would be segmented twice. Refusing to run."
    return
}
println "[detect] selected ${selected.size()} '${TISSUE_CLASS}' annotation(s) — nothing else."

// ---- 3. clear existing detections so a re-run replaces, never accumulates --
int before = getDetectionObjects().size()
if (before > 0) {
    println "[detect] clearing ${before} existing detections"
    removeDetections()
    resetSelection()
    selectObjectsByClassification(TISSUE_CLASS)
}

// ---- park non-detection children of Tissue at the root ----------------------
// (after removeDetections, so only annotations etc. remain). Their subtrees move with them.
def root = getCurrentHierarchy().getRootObject()
def parked = new IdentityHashMap()   // child -> original Tissue parent
selected.each { t -> new ArrayList(t.getChildObjects()).each { c -> parked[c] = t } }
if (!parked.isEmpty()) {
    root.addChildObjects(new ArrayList(parked.keySet()))
    fireHierarchyUpdate()
    println "[detect] parked ${parked.size()} annotation(s) at root during detection: " +
            parked.keySet().countBy { it.getPathClass()?.toString() ?: 'Unclassified' }
}
def restoreParked = {
    parked.groupBy { it.value }.each { t, entries -> t.addChildObjects(new ArrayList(entries.keySet())) }
    fireHierarchyUpdate()
    if (!parked.isEmpty()) println "[detect] restored ${parked.size()} annotation(s) under ${TISSUE_CLASS}"
}

// ---- 4. run InstanSeg ----------------------------
//def presentNames = server.getMetadata().getChannels()
//    .collect { it?.getName() }.findAll { it && it.trim() }
//def extractors = presentNames.collect { ColorTransforms.createChannelExtractor(it) }
//println "[detect] ${presentNames.size()} input channels, model ${MODEL_PATH.split('/')[-1]}"
//
long t0 = System.currentTimeMillis()
def results
try {
results = qupath.ext.instanseg.core.InstanSeg.builder()
    .modelPath(MODEL_PATH)
    .device(DEVICE)
    .inputChannels(INPUT_CHANNELS.collect { ColorTransforms.createChannelExtractor(it) })
    .tileDims(TILE_DIMS)
    .interTilePadding(PAD)
    .nThreads(N_THREADS)
    .makeMeasurements(MAKE_MEASUREMENTS)
    .randomColors(false)
    .outputType("Default")
    .build()
    .detectObjects()
} finally {
    restoreParked()   // also runs if detection throws, so the in-memory hierarchy is never left altered
}
fireHierarchyUpdate()
// detectObjects() catches its own exceptions and per-tile failures; it never throws. Check the
// returned record so a failed or partial run is not saved over the cleared detections.
println "[detect] InstanSeg: ${results}"
if (results.wasInterrupted() || results.nTilesFailed() > 0 || results.nObjectsDetected() == 0) {
    throw new IllegalStateException("[detect] InstanSeg run incomplete (interrupted=${results.wasInterrupted()}, " +
            "failed tiles=${results.nTilesFailed()}/${results.nTilesProcessed()}, objects=${results.nObjectsDetected()}). NOT saving.")
}

def cells = getDetectionObjects()
println String.format("[detect] %,d detections in %.0fs", cells.size(),
        (System.currentTimeMillis() - t0) / 1000.0)

// ---- 5. verify no duplicates ----------------------------------------------
// Grid-bucket the centroids so this is O(n) rather than O(n^2) on ~2M cells.
def cal = server.getPixelCalibration()
double pxSize = cal.hasPixelSizeMicrons() ? cal.getAveragedPixelSizeMicrons() : 1.0
double rPx = DUP_RADIUS_UM / pxSize
// Plain map (not withDefault): withDefault's get() inserts missing keys as a side effect, which
// caused a ConcurrentModificationException when neighbor lookups mutated the map mid-iteration.
// Pair ordering uses a unique per-object index (identityHashCode is not unique and can hide a pair).
def cellList = new ArrayList(cells)
def buckets = [:]
for (int i = 0; i < cellList.size(); i++) {
    def roi = cellList[i].getROI()
    if (roi == null) continue
    def key = [(int) (roi.getCentroidX() / rPx), (int) (roi.getCentroidY() / rPx)]
    buckets.computeIfAbsent(key) { [] } << i
}
int dupPairs = 0
for (entry in new ArrayList(buckets.entrySet())) {
    def (bx, by) = entry.key
    def near = []
    for (dx in -1..1) for (dy in -1..1) near.addAll(buckets.getOrDefault([bx + dx, by + dy], []))
    for (int ia in entry.value) {
        def ra = cellList[ia].getROI()
        for (int ib in near) {
            if (ib <= ia) continue
            def rb = cellList[ib].getROI()
            double d = Math.hypot(ra.getCentroidX() - rb.getCentroidX(),
                                  ra.getCentroidY() - rb.getCentroidY())
            if (d > rPx) continue
            double aa = ra.getArea(), ab = rb.getArea()
            if (Math.abs(aa - ab) / Math.max(Math.max(aa, ab), 1e-9) < DUP_AREA_TOL) dupPairs++
        }
    }
}
if (dupPairs > 0) {
    println "ERROR: ${dupPairs} duplicate detection pair(s) found after detection."
    println "       Two objects share a centroid within ${DUP_RADIUS_UM} um with matching area."
    throw new IllegalStateException("[detect] ${dupPairs} duplicate pairs - detections NOT saved.")
} else {
    println "[detect] duplicate check: PASS (0 pairs within ${DUP_RADIUS_UM} um)"
}
getProjectEntry().saveImageData(getCurrentImageData())
println "[detect] Done."
