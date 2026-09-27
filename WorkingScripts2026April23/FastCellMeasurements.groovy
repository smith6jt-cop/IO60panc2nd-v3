import qupath.lib.analysis.features.ObjectMeasurements
import qupath.lib.analysis.features.ObjectMeasurements.Measurements
import qupath.lib.analysis.features.ObjectMeasurements.Compartments
import java.util.concurrent.ForkJoinPool
import java.util.concurrent.atomic.AtomicInteger

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
selectDetections();
addShapeMeasurements("AREA", "LENGTH", "CIRCULARITY", "SOLIDITY", "MAX_DIAMETER", "MIN_DIAMETER", "NUCLEUS_CELL_RATIO")

println("Sorting ${total} detections spatially...")

// Sort by tile-aligned blocks to maximize cache hits
detections.sort { it.getROI().getCentroidY() * 100000 + it.getROI().getCentroidX() }

println("Measuring ${total} detections with 12 threads...")

def counter = new AtomicInteger(0)
def startTime = System.currentTimeMillis()
def lastReport = new AtomicInteger(0)
def reportInterval = 50000

def pool = new ForkJoinPool(12)
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
            println("  ${count}/${total} (${pct}%) — ${String.format('%.0f', rate)}/sec — ~${String.format('%.1f', remaining / 60)} min remaining")
        }
    }
}.get()
pool.shutdown()

def totalTime = (System.currentTimeMillis() - startTime) / 1000.0
println("Done. ${total} detections in ${String.format('%.1f', totalTime / 60)} min (${String.format('%.0f', total / totalTime)}/sec)")

fireHierarchyUpdate()