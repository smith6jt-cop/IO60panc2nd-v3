/*
 * 5_normalize_channels.groovy
 * Compute per-image anchors (p50, p99) for the object-classifier
 * channels (M2Gal3, CD68, Iba1) from pixels inside the reference mask, then
 * add normalized copies of the per-detection intensity measurements.  Normalize channels for pixel classifiers
 * creating the objects: CD34: Blood_vessel, Podoplanin: Lymphatic, B3TUBB: Nerve.
 *
 * Method: robust rescaling.
 *   norm = (raw - p50) / max(p99 - p50, eps)
 *
 *   p50 = 50th percentile of channel intensity across reference-mask pixels
 *   p99 = 99th percentile of channel intensity across reference-mask pixels
 *   Reference mask = full Tissue union (no exclusions).
 *
 *
 * Output:
 *   - per-detection measurements "<comp>: <chan>: Mean (norm)" for
 *     comp in {Cell, Cytoplasm, Membrane, Nucleus} and
 *     chan in OBJECT_CLF_CHANNELS.
 *   - results/norm_stats.csv with per-image anchors for audit.
 *
 * ---------------------------------------------------------------------------
 * THREADING (this version)
 *   Both phases were single-threaded in the original; this rewrite saturates
 *   the box within each image, so "Run for project" (sequential across the 22
 *   images, parallel within each) fully uses the hardware -- same execution
 *   model as stages 1-4.
 *
 *   Phase A (percentile anchors): the in-mask pixel coordinates are enumerated
 *     ONCE and reused, then the 6 classifier channels are sampled/sorted in
 *     parallel (was 7 full mask-grid walks + sequential channels).
 *   Phase B (write norm columns): the per-detection write loop runs over a
 *     ForkJoinPool. Each detection writes only its own measurement list, so
 *     there is no cross-object contention (identical safety guarantee to
 *     stage 2's parallel addIntensityMeasurements pass).
 *
 *   Numerical output is unchanged: same downsample, same reference mask, same
 *   percentile-index formula, same set of writes -- only the execution order
 *   differs, so norm_stats.csv and the (norm) columns match the serial run.
 *
 *   N_THREADS defaults to 32 to match stage 2's documented JVM-tuning note
 *   (>=64 threads with full heap usage was the instability failure mode on the
 *   64-core Threadripper). Phase B here is low-heap (no tile decode, just
 *   measurement-map ops), so 32 already saturates the physical cores with
 *   headroom; it can be raised toward 64 if desired.
 */

evaluate(new File(buildFilePath(PROJECT_BASE_DIR, 'scripts', 'lib', 'common.groovy')))

import qupath.lib.regions.RegionRequest
import qupath.lib.images.servers.TransformedServerBuilder
import qupath.lib.roi.GeometryTools

import java.awt.image.BufferedImage
import java.awt.Color
import java.awt.geom.AffineTransform

import java.util.concurrent.ForkJoinPool
import java.util.concurrent.ConcurrentHashMap
import java.util.concurrent.atomic.AtomicInteger

// ==================== PARAMETERS ====================
double SAMPLE_DOWNSAMPLE = 16.0   // pixel-sampling downsample for percentile estimation
double P_LOW            = 0.50
double P_HIGH           = 0.99
def COMPARTMENTS        = ['Cell', 'Cytoplasm', 'Membrane', 'Nucleus']
def STATS_TO_NORMALIZE  = ['Mean']       // extend with 'Max' etc. if desired
int N_THREADS           = 32             // see THREADING note in header
int REPORT_INTERVAL     = 200000         // detections between progress lines (phase B)
// ====================================================

def server = getCurrentServer()
if (server == null) { println "[05] ERROR: no image open."; return }

def chanNames = server.getMetadata().getChannels().collect { it.getName() }
def chanIdx = OBJECT_CLF_CHANNELS.collectEntries { [(it): chanNames.indexOf(it)] }
def missing = chanIdx.findAll { k, v -> v < 0 }.keySet()
if (!missing.isEmpty()) {
    println "[05] ERROR: classifier channel(s) missing: ${missing}"
    println "[05]        Check stage 00 alias map and rerun."
    return
}
println "[05] Target channels and indices: ${chanIdx}"

// ---------- Build reference mask: full Tissue union ----------
def tissueAnns = getAnnotationObjects().findAll {
    it.getPathClass()?.toString() == CLASS_TISSUE && it.getROI() != null
}
if (tissueAnns.isEmpty()) { println "[05] ERROR: no ${CLASS_TISSUE} annotations."; return }

def tissueGeoms = tissueAnns.collect { it.getROI().getGeometry() }
def refGeom = (tissueGeoms.size() == 1) ? tissueGeoms[0] : safeUnionPolys(tissueGeoms)
if (refGeom == null || refGeom.isEmpty()) {
    println "[05] ERROR: tissue union is empty."; return
}
println "[05] Reference mask = ${CLASS_TISSUE} (full; no exclusions)"
def refArea = refGeom.getArea()
println String.format("[05] Reference mask area: %.1f px^2 (%.1f mm^2)",
        refArea,
        refArea * (server.getPixelCalibration().getPixelWidthMicrons() ?: 1.0)
              * (server.getPixelCalibration().getPixelHeightMicrons() ?: 1.0) / 1_000_000.0)

// ---------- Rasterize reference mask at SAMPLE_DOWNSAMPLE ----------
int maskW = (int) Math.ceil(server.getWidth()  / SAMPLE_DOWNSAMPLE)
int maskH = (int) Math.ceil(server.getHeight() / SAMPLE_DOWNSAMPLE)
def maskImg = new BufferedImage(maskW, maskH, BufferedImage.TYPE_BYTE_GRAY)
def g2 = maskImg.createGraphics()
g2.setColor(Color.BLACK); g2.fillRect(0, 0, maskW, maskH)
g2.setColor(Color.WHITE)
g2.setTransform(AffineTransform.getScaleInstance(1.0 / SAMPLE_DOWNSAMPLE, 1.0 / SAMPLE_DOWNSAMPLE))
def shape = GeometryTools.geometryToShape(refGeom)
g2.fill(shape)
g2.dispose()
def maskRaster = maskImg.getRaster()

// ---------- Enumerate in-mask pixel coordinates ONCE (reused by all channels) ----------
// Pass 1: count. Pass 2: fill flat coordinate arrays. Replaces the original
// per-channel full-grid rescan (7 mask walks -> 2), and each channel now only
// visits in-mask pixels.
int[] mbuf = new int[1]
long maskPxCountL = 0
for (int y = 0; y < maskH; y++) {
    for (int x = 0; x < maskW; x++) {
        maskRaster.getPixel(x, y, mbuf)
        if (mbuf[0] > 0) maskPxCountL++
    }
}
if (maskPxCountL > Integer.MAX_VALUE) {
    println "[05] ERROR: mask pixel count exceeds Integer.MAX_VALUE (${maskPxCountL}); raise SAMPLE_DOWNSAMPLE."
    return
}
int maskPxCount = (int) maskPxCountL
int[] maskXs = new int[maskPxCount]
int[] maskYs = new int[maskPxCount]
int mi = 0
for (int y = 0; y < maskH; y++) {
    for (int x = 0; x < maskW; x++) {
        maskRaster.getPixel(x, y, mbuf)
        if (mbuf[0] > 0) { maskXs[mi] = x; maskYs[mi] = y; mi++ }
    }
}
println "[05] Mask raster: ${maskW}x${maskH} at downsample=${SAMPLE_DOWNSAMPLE}, ${maskPxCount} in-mask pixels"

// ---------- Compute p_low and p_high per channel (parallel across channels) ----------
def plane = tissueAnns[0].getROI().getImagePlane()
def stats = new ConcurrentHashMap()

int nChan = OBJECT_CLF_CHANNELS.size()
def phaseAPool = new ForkJoinPool(Math.max(1, Math.min(N_THREADS, nChan)))
phaseAPool.submit {
    OBJECT_CLF_CHANNELS.parallelStream().forEach { name ->
        int ci = chanIdx[name]
        def single = new TransformedServerBuilder(server).extractChannels(ci).build()
        def req = RegionRequest.createInstance(
            single.getPath(), SAMPLE_DOWNSAMPLE,
            0, 0, single.getWidth(), single.getHeight(), plane.getZ(), plane.getT()
        )
        def img = single.readRegion(req)
        def raster = img.getRaster()
        int rw = raster.getWidth()
        int rh = raster.getHeight()

        // Sample the precomputed in-mask coordinates; skip any that fall outside
        // this channel raster (matches the original min(maskW,rasterW) clamp).
        float[] values = new float[maskPxCount]
        int k = 0
        for (int i = 0; i < maskPxCount; i++) {
            int x = maskXs[i]
            int y = maskYs[i]
            if (x >= rw || y >= rh) continue
            values[k++] = raster.getSampleFloat(x, y, 0)
        }
        int n = k
        if (n == 0) {
            println "[05] WARN: no pixels sampled for ${name}"
            stats[name] = [p_low: 0.0d, p_high: 1.0d, n: 0]
            return
        }
        // Sort the populated prefix only.
        float[] sorted = (n == values.length) ? values : java.util.Arrays.copyOf(values, n)
        java.util.Arrays.sort(sorted)
        double pLow  = sorted[Math.min(n - 1, (int)(n * P_LOW))]
        double pHigh = sorted[Math.min(n - 1, (int)(n * P_HIGH))]
        stats[name] = [p_low: pLow, p_high: pHigh, n: n]
        println String.format("[05]   %-8s n=%d  p%.0f=%.3f  p%.0f=%.3f  span=%.3f",
                name, n, P_LOW*100, pLow, P_HIGH*100, pHigh, pHigh - pLow)
    }
}.get()
phaseAPool.shutdown()

// ---------- Precompute write specs (src -> dst, pLow, span) ----------
// Independent of detections; built once so the hot loop does no string work.
def writeSpecs = []
OBJECT_CLF_CHANNELS.each { name ->
    def s = stats[name]
    if (s == null) return
    double pLow = s.p_low
    double span = Math.max(s.p_high - pLow, 1e-6)
    COMPARTMENTS.each { comp ->
        STATS_TO_NORMALIZE.each { stat ->
            writeSpecs << [
                src : "${comp}: ${name}: ${stat}".toString(),
                dst : "${comp}: ${name}: ${stat} (norm)".toString(),
                pLow: pLow,
                span: span
            ]
        }
    }
}

// ---------- Write normalized measurements onto detections (parallel) ----------
def detections = getDetectionObjects()
if (detections.isEmpty()) {
    println "[05] WARN: no detections - nothing to normalize."
} else if (writeSpecs.isEmpty()) {
    println "[05] WARN: no valid channel stats - nothing to write."
} else {
    // Serial write-verification probe (once), before the parallel pass. Catches
    // API-version mismatches in putMeas that would otherwise silently produce
    // zero norm columns. Scans a few detections until a non-NaN source is found.
    boolean verified = false
    outer:
    for (def det : detections) {
        for (def spec : writeSpecs) {
            double v = getMeas(det, spec.src)
            if (Double.isNaN(v)) continue
            double expected = (v - spec.pLow) / spec.span
            putMeas(det, spec.dst, expected)
            double readBack = getMeas(det, spec.dst)
            if (Double.isNaN(readBack)) {
                throw new RuntimeException(
                    "[05] putMeas failed: wrote '${spec.dst}' but read-back returned NaN. " +
                    "Check common.groovy putMeas against this QuPath version's MeasurementList API.")
            }
            println "[05] Write verification: '${spec.dst}' = ${readBack} (expected ${expected})"
            verified = true
            break outer
        }
    }
    if (!verified) {
        println "[05] WARN: every source measurement was NaN across all detections - no norm columns written."
    } else {
        int total = detections.size()
        def nUpdated = new AtomicInteger(0)
        def counter  = new AtomicInteger(0)
        def lastReport = new AtomicInteger(0)
        long startTime = System.currentTimeMillis()
        println "[05] Writing norm columns for ${total} detections with ${N_THREADS} threads (${writeSpecs.size()} columns each)..."

        def pool = new ForkJoinPool(N_THREADS)
        pool.submit {
            detections.parallelStream().forEach { det ->
                for (def spec : writeSpecs) {
                    double v = getMeas(det, spec.src)
                    if (Double.isNaN(v)) continue
                    putMeas(det, spec.dst, (v - spec.pLow) / spec.span)
                    nUpdated.incrementAndGet()
                }
                int count = counter.incrementAndGet()
                int prev = lastReport.get()
                if (count - prev >= REPORT_INTERVAL && lastReport.compareAndSet(prev, count)) {
                    double elapsed = (System.currentTimeMillis() - startTime) / 1000.0
                    double rate = count / Math.max(elapsed, 1e-6)
                    double remaining = (total - count) / Math.max(rate, 1e-6)
                    println String.format("[05]   %d/%d (%.1f%%) - %.0f det/s - ~%.1f min remaining",
                            count, total, count * 100.0 / total, rate, remaining / 60.0)
                }
            }
        }.get()
        pool.shutdown()

        println "[05] Wrote ${nUpdated.get()} normalized measurement values across ${total} detections"
    }
}

// ---------- Persist anchors ----------
def outDir = buildFilePath(PROJECT_BASE_DIR, 'results'); mkdirs(outDir)
def csvFile = new File(buildFilePath(outDir, 'norm_stats.csv'))
boolean newFile = !csvFile.exists()
csvFile.withWriterAppend('UTF-8') { w ->
    if (newFile) w.write("image,channel,p_low,p_high,span,n_pixels,p_low_quantile,p_high_quantile\n")
    def imgName = (getProjectEntry()?.getImageName() ?: 'unknown').replaceAll(',', '_')
    OBJECT_CLF_CHANNELS.each { name ->
        def s = stats[name]; if (s == null) return
        w.write("${imgName},${name},${s.p_low},${s.p_high},${s.p_high - s.p_low},${s.n},${P_LOW},${P_HIGH}\n")
    }
}
println "[05] Wrote anchors to ${csvFile}"
fireHierarchyUpdate()