/*
 * 05_cell_measurements.groovy
 *
 * Stage 05: add shape and intensity measurements to all detection objects.
 *
 * Two modes, selectable via MEASURE_CHANNELS:
 *   null (default) -> measure all 59 channels. Heavy but matches prior runs.
 *   List<String>    -> measure only the listed channels (via a TransformedServer
 *                      per channel). Significantly faster when downstream only
 *                      needs a subset (e.g. the object-classifier channels).
 *
 * Output: per-detection measurements of the form
 *   "<compartment>: <channel>: <stat>"   e.g. "Cell: CD34: Mean"
 *
 * Stage 06 (channel normalization) adds "(norm)" variants of the Mean columns
 * for the OBJECT_CLF_CHANNELS defined in common.groovy.
 */

evaluate(new File(buildFilePath(PROJECT_BASE_DIR, 'scripts', 'lib', 'common.groovy')))

import qupath.lib.analysis.features.ObjectMeasurements
import qupath.lib.analysis.features.ObjectMeasurements.Measurements
import qupath.lib.analysis.features.ObjectMeasurements.Compartments
import qupath.lib.images.servers.TransformedServerBuilder
import java.util.concurrent.ForkJoinPool
import java.util.concurrent.atomic.AtomicInteger

// ==================== PARAMETERS ====================
double DOWNSAMPLE    = 1.0
int    N_THREADS     = 12
// Set to null for all channels, or to a list of channel names to restrict.
// Example (fast): MEASURE_CHANNELS = OBJECT_CLF_CHANNELS
List<String> MEASURE_CHANNELS = null
// ====================================================

def server = getCurrentServer()
def detections = getDetectionObjects()
int total = detections.size()
if (total == 0) { println "[05] No detections. Exiting."; return }

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

selectDetections()
addShapeMeasurements("AREA", "LENGTH", "CIRCULARITY", "SOLIDITY",
                     "MAX_DIAMETER", "MIN_DIAMETER", "NUCLEUS_CELL_RATIO")

println "[05] Sorting ${total} detections spatially..."
detections.sort { it.getROI().getCentroidY() * 100000 + it.getROI().getCentroidX() }

// ---------- Choose servers to iterate over ----------
def serversToMeasure = []
if (MEASURE_CHANNELS == null || MEASURE_CHANNELS.isEmpty()) {
    serversToMeasure << server
    println "[05] Measuring all ${server.getMetadata().getChannels().size()} channels"
} else {
    def chanNames = server.getMetadata().getChannels().collect { it.getName() }
    MEASURE_CHANNELS.each { name ->
        int ci = chanNames.indexOf(name)
        if (ci < 0) { println "[05] WARN: channel '${name}' not found - skipping"; return }
        def single = new TransformedServerBuilder(server).extractChannels(ci).build()
        serversToMeasure << single
    }
    println "[05] Measuring ${serversToMeasure.size()} selected channels: ${MEASURE_CHANNELS}"
}

// ---------- Run measurements ----------
long startTime = System.currentTimeMillis()
int reportInterval = 50000

serversToMeasure.eachWithIndex { s, si ->
    def counter = new AtomicInteger(0)
    def lastReport = new AtomicInteger(0)
    println "[05] Server ${si+1}/${serversToMeasure.size()}: ${s.getMetadata().getChannels().collect{ it.getName() }}"

    def pool = new ForkJoinPool(N_THREADS)
    pool.submit {
        detections.parallelStream().forEach { det ->
            ObjectMeasurements.addIntensityMeasurements(
                s, det, DOWNSAMPLE, measurements, compartments
            )
            int count = counter.incrementAndGet()
            int prev = lastReport.get()
            if (count - prev >= reportInterval && lastReport.compareAndSet(prev, count)) {
                double elapsed = (System.currentTimeMillis() - startTime) / 1000.0
                double rate = count / elapsed
                double remaining = (total - count) / rate
                println String.format("[05]   %d/%d (%.1f%%) - %.0f/s - ~%.1f min left",
                    count, total, count * 100.0 / total, rate, remaining / 60)
            }
        }
    }.get()
    pool.shutdown()
}

double totalTime = (System.currentTimeMillis() - startTime) / 1000.0
println String.format("[05] Done. ${total} detections in %.1f min", totalTime / 60.0)
fireHierarchyUpdate()
