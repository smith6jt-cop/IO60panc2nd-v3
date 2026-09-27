import qupath.lib.objects.PathObjects
import qupath.lib.roi.GeometryTools
import org.locationtech.jts.geom.Geometry
import org.locationtech.jts.geom.Polygon
import org.locationtech.jts.geom.LinearRing
import org.locationtech.jts.geom.GeometryFactory
import org.locationtech.jts.geom.util.GeometryFixer
import org.locationtech.jts.operation.union.UnaryUnionOp
import org.locationtech.jts.operation.overlayng.OverlayNGRobust
import java.util.concurrent.ForkJoinPool
import java.util.concurrent.Callable
import java.util.concurrent.atomic.AtomicLong
import java.util.concurrent.atomic.AtomicInteger

// --- Parameters ---
String targetClassName = "Other: Gal3Norm"
double bridgeUm          = 3.0       // dilation per cell; cells within 2*bridgeUm edge-to-edge merge
double neighborRadiusUm  = 50.0      // density pre-filter search radius
int    minNeighbors      = 10        // required Gal3+ neighbors within neighborRadiusUm; 0 disables filter
double minAreaUm2        = 2500.0    // drop output features smaller than this (one 50-µm pixel = 2500 µm²)
boolean closeOnly        = false     // dilation-only merge (no negative buffer)
int    tileSizePx        = 1500
int    nThreads          = 32
int    mergeBatchSize    = 8

def imageData = getCurrentImageData()
def cal = imageData.getServer().getPixelCalibration()
if (!cal.hasPixelSizeMicrons()) { println "ERROR: image has no pixel size"; return }
double pxSize     = cal.getAveragedPixelSizeMicrons()
double bufferPx   = bridgeUm / pxSize
double minAreaPx2 = minAreaUm2 / (pxSize * pxSize)

def targetClass = getPathClass(targetClassName)
def allCells = getDetectionObjects().findAll { it.getPathClass() == targetClass }
println "Cells of class ${targetClassName}: ${allCells.size()}"
if (allCells.isEmpty()) return

// --- Density pre-filter: drop cells with too few Gal3+ neighbors within neighborRadiusUm ---
def cells
if (minNeighbors > 0) {
    long tFilter = System.currentTimeMillis()
    double rPx  = neighborRadiusUm / pxSize
    double rPx2 = rPx * rPx
    int gridPx  = (int) Math.ceil(rPx)

    int n = allCells.size()
    double[] cx = new double[n]
    double[] cy = new double[n]
    def grid = new HashMap<String, List<Integer>>()
    for (int i = 0; i < n; i++) {
        def roi = allCells[i].getROI()
        cx[i] = roi.getCentroidX()
        cy[i] = roi.getCentroidY()
        String key = ((int) Math.floor(cx[i] / gridPx)) + "_" + ((int) Math.floor(cy[i] / gridPx))
        def bucket = grid.get(key)
        if (bucket == null) { bucket = new ArrayList<Integer>(); grid.put(key, bucket) }
        bucket.add(i)
    }

    int[] counts = new int[n]
    def pool0 = new ForkJoinPool(nThreads)
    def tasks = (0..<n).collect { int i ->
        { ->
            int gx = (int) Math.floor(cx[i] / gridPx)
            int gy = (int) Math.floor(cy[i] / gridPx)
            int c = 0
            for (int dx = -1; dx <= 1; dx++) {
                for (int dy = -1; dy <= 1; dy++) {
                    def bucket = grid.get((gx+dx) + "_" + (gy+dy))
                    if (bucket == null) continue
                    for (int j : bucket) {
                        if (j == i) continue
                        double ex = cx[j] - cx[i]
                        double ey = cy[j] - cy[i]
                        if (ex*ex + ey*ey <= rPx2) c++
                    }
                }
            }
            counts[i] = c
            return null
        } as Callable
    }
    pool0.invokeAll(tasks).each { it.get() }
    pool0.shutdown()

    def sortedCounts = counts.toList().toSorted()
    def pct = { double p -> sortedCounts[Math.min(sortedCounts.size()-1, (int)(p * sortedCounts.size()))] }
    println "Neighbor-count distribution (p10/p50/p90/p99/max): " +
            "${pct(0.1)} / ${pct(0.5)} / ${pct(0.9)} / ${pct(0.99)} / ${sortedCounts[-1]}"

    cells = new ArrayList()
    for (int i = 0; i < n; i++) if (counts[i] >= minNeighbors) cells << allCells[i]
    println "Density filter kept ${cells.size()}/${n} cells (>= ${minNeighbors} neighbors within ${neighborRadiusUm} µm) " +
            "in ${(System.currentTimeMillis()-tFilter)/1000}s"
} else {
    cells = allCells
}
if (cells.isEmpty()) return
def plane = cells[0].getROI().getImagePlane()

// --- Partition by centroid into tiles ---
def tiles = [:].withDefault {