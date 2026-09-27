/*
 * 2_fast_cell_measurements.groovy
 *
 * Shape + intensity measurements on all detections (post Instanseg cell segmentation).
 *
 * Threading: 32 threads by default. The bottleneck for 59-channel QPTIFF is
 * tile read/decompression, which is CPU-bound. On the 64-core Threadripper,
 * 32 leaves headroom against the system-instability zone documented in the
 * QuPath JVM tuning notes (>=64 with full heap usage was the failure mode).
 * Going below ~24 leaves performance on the table.
 *
 * Optional: PREWARM_TILES = true triggers a coarse parallel pre-pass over
 * the Tissue bounds to decode tiles into the QuPath cache, so the per-cell
 * readRegion calls hit cache rather than disk on the second pass. This is
 * worth ~10-20% on cold caches; toggle off if memory pressure is a concern.
 */

import qupath.lib.analysis.features.ObjectMeasurements
import qupath.lib.analysis.features.ObjectMeasurements.Measurements
import qupath.lib.analysis.features.ObjectMeasurements.Compartments
import qupath.lib.regions.RegionRequest
import java.util.concurrent.ForkJoinPool
import java.util.concurrent.atomic.AtomicInteger

// ==================== TUNING ====================
int N_THREADS        = 18
boolean PREWARM_TILES = true
int PREWARM_TILE_PX   = 2048    // size of pre-warm tiles in pixels
// ================================================

def server = getCurrentServer()
def downsample = 1.0

def compartments = [
    Compartments.NUCLEUS,
    Compartments.CYTOPLASM,
    Compartments.MEMBRANE,
    Compartments.CELL
]

def measurements = [
    Measurements.MEAN,
    Measurements.MIN,
    Measurements.MAX,
    Measurements.STD_DEV
]

def detections = getDetectionObjects()
def total = detections.size()
selectDetections()
addShapeMeasurements("AREA", "LENGTH", "CIRCULARITY", "SOLIDITY", "MAX_DIAMETER", "MIN_DIAMETER", "NUCLEUS_CELL_RATIO")

println("Sorting ${total} detections spatially...")
detections.sort { it.getROI().getCentroidY() * 100000 + it.getROI().getCentroidX() }

// ----- Optional tile pre-warm -----
if (PREWARM_TILES) {
    println("Pre-warming tile cache...")
    long t0 = System.currentTimeMillis()
    def tissues = getAnnotationObjects().findAll { it.getPathClass()?.toString() == 'Tissue' }
    if (tissues.isEmpty()) {
        println("  No Tissue annotations - skipping pre-warm")
    } else {
        // Compute bounding box across all Tissue annotations
        double minX = Double.POSITIVE_INFINITY, minY = Double.POSITIVE_INFINITY
        double maxX = Double.NEGATIVE_INFINITY, maxY = Double.NEGATIVE_INFINITY
        tissues.each { ann ->
            def env = ann.getROI().getGeometry().getEnvelopeInternal()
            if (env.getMinX() < minX) minX = env.getMinX()
            if (env.getMinY() < minY) minY = env.getMinY()
            if (env.getMaxX() > maxX) maxX = env.getMaxX()
            if (env.getMaxY() > maxY) maxY = env.getMaxY()
        }
        def plane = tissues[0].getROI().getImagePlane()
        def tileRequests = []
        for (int y = (int) minY; y < maxY; y += PREWARM_TILE_PX) {
            for (int x = (int) minX; x < maxX; x += PREWARM_TILE_PX) {
                int w = (int) Math.min(PREWARM_TILE_PX, maxX - x)
                int h = (int) Math.min(PREWARM_TILE_PX, maxY - y)
                tileRequests << RegionRequest.createInstance(server.getPath(), 1.0,
                        x, y, w, h, plane.getZ(), plane.getT())
            }
        }
        println("  ${tileRequests.size()} pre-warm tiles")
        def warmCounter = new AtomicInteger(0)
        def warmPool = new ForkJoinPool(N_THREADS)
        warmPool.submit {
            tileRequests.parallelStream().forEach { req ->
                try { server.readRegion(req) } catch (Throwable t) { /* skip on error */ }
                int n = warmCounter.incrementAndGet()
                if (n % 200 == 0) println("    ${n}/${tileRequests.size()} tiles")
            }
        }.get()
        warmPool.shutdown()
        println String.format("  Pre-warm done in %.1f min", (System.currentTimeMillis() - t0) / 60000.0)
    }
}

// ----- Main per-cell measurement loop -----
println("Measuring ${total} detections with ${N_THREADS} threads...")

def counter = new AtomicInteger(0)
def startTime = System.currentTimeMillis()
def lastReport = new AtomicInteger(0)
def reportInterval = 50000

def pool = new ForkJoinPool(N_THREADS)
pool.submit {
    detections.parallelStream().forEach { detection ->
        ObjectMeasurements.addIntensityMeasurements(
            server, detection, downsample, measurements, compartments
        )
        def count = counter.incrementAndGet()
        def prev = lastReport.get()
        if (count - prev >= reportInterval && lastReport.compareAndSet(prev, count)) {
            def elapsed = (System.currentTimeMillis() - startTime) / 1000.0
            def rate = count / elapsed
            def remaining = (total - count) / rate
            def pct = String.format("%.1f", count * 100.0 / total)
            println("  ${count}/${total} (${pct}%) - ${String.format('%.0f', rate)}/sec - ~${String.format('%.1f', remaining / 60)} min remaining")
        }
    }
}.get()
pool.shutdown()

def totalTime = (System.currentTimeMillis() - startTime) / 1000.0
println("Done. ${total} detections in ${String.format('%.1f', totalTime / 60)} min (${String.format('%.0f', total / totalTime)}/sec)")

fireHierarchyUpdate()