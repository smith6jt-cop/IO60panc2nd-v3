// Group classified cell detections of one class into annotations at a target resolution.
// Tiled parallel buffer-union: partition cells by centroid into tiles, union per tile in
// parallel, then global union, then negative buffer.

import qupath.lib.objects.PathObjects
import qupath.lib.roi.GeometryTools
import org.locationtech.jts.geom.util.GeometryFixer
import org.locationtech.jts.operation.overlayng.OverlayNGRobust
import java.util.concurrent.ForkJoinPool
import java.util.concurrent.Callable
import java.util.concurrent.atomic.AtomicLong

// --- Parameters ---
String targetClassName = "Other: Gal3Norm"
double resolutionUm    = 4.0     // collaborator's minimum resolution
double minAreaUm2      = 0.0      // drop groups smaller than this (0 keeps all)
boolean closeOnly      = true     // true: morphological closing; false: leave dilated
int tileSizePx         = 5000     // tile edge in pixels for parallel partitioning
int nThreads           = 52

def imageData = getCurrentImageData()
def cal = imageData.getServer().getPixelCalibration()
if (!cal.hasPixelSizeMicrons()) {
    println "ERROR: image has no pixel size in microns"
    return
}
double pxSize     = cal.getAveragedPixelSizeMicrons()
double bufferPx   = (resolutionUm / 2.0) / pxSize
double minAreaPx2 = minAreaUm2 / (pxSize * pxSize)

def targetClass = getPathClass(targetClassName)
def cells = getDetectionObjects().findAll { it.getPathClass() == targetClass }
println "Cells of class ${targetClassName}: ${cells.size()}"
if (cells.isEmpty()) return

def plane = cells[0].getROI().getImagePlane()

// --- Partition cells by centroid into tiles ---
def tiles = [:].withDefault { [] }
cells.each { c ->
    def roi = c.getROI()
    int ix = (int) Math.floor(roi.getCentroidX() / tileSizePx)
    int iy = (int) Math.floor(roi.getCentroidY() / tileSizePx)
    tiles["${ix}_${iy}"] << c
}
println "Non-empty tiles: ${tiles.size()}"

// --- Per-tile buffered union in parallel ---
def nOk = new AtomicLong(0)
def nSkipped = new AtomicLong(0)
def pool = new ForkJoinPool(nThreads)
def perTileResults = []

try {
    def tasks = tiles.values().collect { tileCells ->
        { ->
            def buffered = new ArrayList(tileCells.size())
            int ok = 0, skipped = 0
            tileCells.each { cell ->
                try {
                    def g = cell.getROI().getGeometry()
                    if (!g.isValid()) g = GeometryFixer.fix(g)
                    def b = g.buffer(bufferPx)
                    if (b != null && !b.isEmpty()) { buffered << b; ok++ }
                    else skipped++
                } catch (Throwable t) {
                    skipped++
                }
            }
            nOk.addAndGet(ok); nSkipped.addAndGet(skipped)
            if (buffered.isEmpty()) return null
            try {
                return OverlayNGRobust.union(buffered)
            } catch (Throwable t) {
                println "  tile union failed: ${t.message}"
                return null
            }
        } as Callable
    }
    perTileResults = pool.invokeAll(tasks)
                         .collect { it.get() }
                         .findAll { it != null && !it.isEmpty() }
} finally {
    pool.shutdown()
}

println "Buffered ok: ${nOk.get()}, skipped: ${nSkipped.get()}"
println "Per-tile unions produced: ${perTileResults.size()}"
if (perTileResults.isEmpty()) return

// --- Global union of per-tile results ---
def unioned
try {
    unioned = OverlayNGRobust.union(perTileResults)
} catch (Throwable t) {
    println "ERROR during global union: ${t.message}"
    return
}
if (unioned == null || unioned.isEmpty()) return

// --- Negative buffer (close) ---
def result = unioned
if (closeOnly) {
    try {
        def fixed = unioned.isValid() ? unioned : GeometryFixer.fix(unioned)
        result = fixed.buffer(-bufferPx)
    } catch (Throwable t) {
        println "ERROR during un-buffer: ${t.message}; leaving dilated"
    }
}
if (result == null || result.isEmpty()) return

// --- Build annotations ---
def newAnnotations = []
for (int i = 0; i < result.getNumGeometries(); i++) {
    def g = result.getGeometryN(i)
    if (g.getArea() < minAreaPx2) continue
    def roi = GeometryTools.geometryToROI(g, plane)
    def ann = PathObjects.createAnnotationObject(roi, targetClass)
    ann.setLocked(true)
    newAnnotations << ann
}

imageData.getHierarchy().addObjects(newAnnotations)
fireHierarchyUpdate()
println "Created ${newAnnotations.size()} group annotations"