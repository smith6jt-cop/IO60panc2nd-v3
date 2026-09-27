/*
 * 05_normalize_channels.groovy
 *
 * Workflow step 6.5 (inserted between FastCellMeasurements and the composite
 * classifier). Compute per-image anchors (p50, p99) for the object-classifier
 * channels (M2Gal3, CD34, Iba1) from pixels inside the reference mask, then
 * add normalized copies of the per-detection intensity measurements.
 *
 * Method: robust rescaling.
 *   norm = (raw - p50) / max(p99 - p50, eps)
 *
 *   p50 = 50th percentile of channel intensity across reference-mask pixels
 *   p99 = 99th percentile of channel intensity across reference-mask pixels
 *   Reference mask = full Tissue union (no exclusions).
 *
 * The reference mask is deliberately the full Tissue area: Islet and
 * Blood_vessel regions are NOT excluded. Rationale: CD34 vessel coverage
 * is too large to reliably represent in the hierarchy in this pipeline
 * (it is written only as GeoJSON in stage 02), and the user has specified
 * a single Tissue-wide baseline. See NORMALIZATION.md for the implications.
 *
 * Output:
 *   - per-detection measurements "<comp>: <chan>: Mean (norm)" for
 *     comp in {Cell, Cytoplasm, Membrane, Nucleus} and
 *     chan in OBJECT_CLF_CHANNELS.
 *   - results/norm_stats.csv with per-image anchors for audit.
 *
 * Run via: Run -> Run for project.
 */

evaluate(new File(buildFilePath(PROJECT_BASE_DIR, 'scripts', 'lib', 'common.groovy')))

import qupath.lib.regions.RegionRequest
import qupath.lib.images.servers.TransformedServerBuilder
import qupath.lib.roi.GeometryTools

import java.awt.image.BufferedImage
import java.awt.Color
import java.awt.geom.AffineTransform

import it.unimi.dsi.fastutil.floats.FloatArrayList

// ==================== PARAMETERS ====================
double SAMPLE_DOWNSAMPLE = 16.0   // pixel-sampling downsample for percentile estimation
double P_LOW            = 0.50
double P_HIGH           = 0.99
def COMPARTMENTS        = ['Cell', 'Cytoplasm', 'Membrane', 'Nucleus']
def STATS_TO_NORMALIZE  = ['Mean']       // extend with 'Max' etc. if desired
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

long maskPxCount = 0
int[] mbuf = new int[1]
for (int y = 0; y < maskH; y++) {
    for (int x = 0; x < maskW; x++) {
        maskRaster.getPixel(x, y, mbuf)
        if (mbuf[0] > 0) maskPxCount++
    }
}
println "[05] Mask raster: ${maskW}x${maskH} at downsample=${SAMPLE_DOWNSAMPLE}, ${maskPxCount} in-mask pixels"

// ---------- Compute p_low and p_high per channel ----------
def plane = tissueAnns[0].getROI().getImagePlane()
def stats = [:]

OBJECT_CLF_CHANNELS.each { name ->
    int ci = chanIdx[name]
    def single = new TransformedServerBuilder(server).extractChannels(ci).build()
    def req = RegionRequest.createInstance(
        single.getPath(), SAMPLE_DOWNSAMPLE,
        0, 0, single.getWidth(), single.getHeight(), plane.getZ(), plane.getT()
    )
    def img = single.readRegion(req)
    def raster = img.getRaster()
    int iw = Math.min(maskW, raster.getWidth())
    int ih = Math.min(maskH, raster.getHeight())

    def values = new FloatArrayList()
    for (int y = 0; y < ih; y++) {
        for (int x = 0; x < iw; x++) {
            maskRaster.getPixel(x, y, mbuf)
            if (mbuf[0] == 0) continue
            values.add(raster.getSampleFloat(x, y, 0))
        }
    }
    int n = values.size()
    if (n == 0) {
        println "[05] WARN: no pixels sampled for ${name}"
        stats[name] = [p_low: 0.0d, p_high: 1.0d, n: 0]
        return
    }
    values.sort()
    double pLow  = values.getFloat(Math.min(n - 1, (int)(n * P_LOW)))
    double pHigh = values.getFloat(Math.min(n - 1, (int)(n * P_HIGH)))
    stats[name] = [p_low: pLow, p_high: pHigh, n: n]
    println String.format("[05]   %-8s n=%d  p%.0f=%.3f  p%.0f=%.3f  span=%.3f",
            name, n, P_LOW*100, pLow, P_HIGH*100, pHigh, pHigh - pLow)
}

// ---------- Write normalized measurements onto detections ----------
def detections = getDetectionObjects()
if (detections.isEmpty()) {
    println "[05] WARN: no detections - nothing to normalize."
} else {
    int nUpdated = 0
    OBJECT_CLF_CHANNELS.each { name ->
        def s = stats[name]
        if (s == null) return
        double pLow = s.p_low
        double span = Math.max(s.p_high - pLow, 1e-6)
        COMPARTMENTS.each { comp ->
            STATS_TO_NORMALIZE.each { stat ->
                String src = "${comp}: ${name}: ${stat}"
                String dst = "${comp}: ${name}: ${stat} (norm)"
                detections.each { det ->
                    double v = getMeas(det, src)
                    if (!Double.isNaN(v)) {
                        putMeas(det, dst, (v - pLow) / span)
                        nUpdated++
                    }
                }
            }
        }
    }
    println "[05] Wrote ${nUpdated} normalized measurement values across ${detections.size()} detections"
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
