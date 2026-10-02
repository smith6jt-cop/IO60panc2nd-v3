/*
 * 3_signed_distance_v3.groovy   (v3 copy of scripts/3_signed_distance.groovy; the original is unchanged)
 *
 * Workflow step 4. Signed distance from each detection centroid to the boundary of each
 * non-ignored annotation class (IsletExpanded excluded). Negative inside, positive outside.
 *
 * v3 changes (nothing else differs from the original):
 *   - these columns are written on EVERY image, so all images share one schema:
 *         Islet, EndocrineCluster, EndocrineSingle, Hormone Fragment, Artifact, pLN, Tissue
 *     A class with no annotation on the image gets NaN on every detection (not a missing column).
 *     Any other non-ignored class present (e.g. Nerve, Blood_vessel) is still computed, as before.
 *   - no image / no detections throws instead of returning (under run_v3 an early return would also be a
 *     failure, but without a message); no Islet annotations is a warning, not a stop -- the other classes
 *     are still computed and the Islet column is NaN.
 *
 * Output column (unchanged from the previous version of this script):
 *   "Signed distance to annotation with <class> <unit>"   e.g. "... with Islet µm"
 * NOTE: QuPath's own command (DistanceTools.detectionToAnnotationDistances, splitClassNames=false)
 * names it "Signed distance to annotation <class> <unit>" (no "with"), and the unit string from
 * PixelCalibration is "µm", not "um". Anything downstream that keys on the column name must match.
 *
 * Semantics match DistanceTools.centroidToBoundsDistance2D (QuPath 0.7.0):
 *   - source point = centroid of the NUCLEUS ROI when the detection is a cell with a nucleus
 *   - distance to the combined boundary of all annotations of that class on the same plane
 *   - EXTERIOR -> +d, INTERIOR -> -d, BOUNDARY -> 0
 *
 * Speed: QuPath's DistanceToPoint is a brute-force scan over every target vertex per detection.
 * The previous version made that tolerable by simplifying the targets by 2 um, which shifted
 * boundaries by up to ~2 um and could flip the sign for cells near islet edges (so the sign no
 * longer agreed with 8_hierarchy_parenting, which uses the full-precision geometry). This version
 * uses JTS IndexedFacetDistance on the exact geometry instead: identical results to the
 * brute-force method, O(log n) per query, no simplification.
 *
 * Because classes are filtered explicitly, running after 4_islet_expander no longer adds an
 * IsletExpanded column.
 */

import qupath.lib.objects.PathObjectTools
import qupath.lib.objects.classes.PathClassTools
import qupath.lib.roi.GeometryTools
import org.locationtech.jts.algorithm.locate.IndexedPointInAreaLocator
import org.locationtech.jts.geom.Coordinate
import org.locationtech.jts.geom.Location
import org.locationtech.jts.geom.Polygonal
import org.locationtech.jts.geom.util.AffineTransformation
import org.locationtech.jts.geom.util.GeometryCombiner
import org.locationtech.jts.operation.distance.IndexedFacetDistance

evaluate(new File(buildFilePath(PROJECT_BASE_DIR, 'scripts', 'lib', 'common.groovy')))

// ==================== PARAMETERS ====================
// Columns every image gets (NaN where the class is absent), in addition to any other class present.
// 'Hormone Fragment' (step 1b, 2026-10-02): an image without one still gets the column, so the context tables agree
List<String> ENSURED_CLASSES = ['Islet', 'EndocrineCluster', 'EndocrineSingle', 'Hormone Fragment', 'Artifact', 'pLN', 'Tissue']
// ====================================================

def imageData = getCurrentImageData()
if (imageData == null) throw new IllegalStateException("[02v3] No image open.")

def detections = getDetectionObjects()
if (detections.isEmpty()) throw new IllegalStateException("[02v3] No detections.")

def islets = getAnnotationObjects().findAll {
    it.getPathClass()?.toString() == CLASS_ISLET && it.getROI() != null
}
if (islets.isEmpty())
    println "[02v3] WARN: no ${CLASS_ISLET} annotations (run 1_panc_objects_v3.groovy first); its column will be NaN."

def annotationClasses = getAnnotationObjects()
    .collect { it.getPathClass() }
    .findAll { it != null && it.isValid() && !PathClassTools.isIgnoredClass(it) && it.toString() != CLASS_ISLET_EXPANDED }
    .toSet()
println "[02v3] Computing signed distances. Classes: ${annotationClasses}"
println "[02v3] Detections: ${detections.size()}"

def cal = imageData.getServer().getPixelCalibration()
if (cal.getPixelWidthUnit() != cal.getPixelHeightUnit())
    throw new IllegalArgumentException("Pixel width/height units differ: ${cal.getPixelWidthUnit()} / ${cal.getPixelHeightUnit()}")
double pixelWidth  = cal.getPixelWidth().doubleValue()
double pixelHeight = cal.getPixelHeight().doubleValue()
String unit = cal.getPixelWidthUnit()
def transform = (pixelWidth == 1 && pixelHeight == 1) ? null : AffineTransformation.scaleInstance(pixelWidth, pixelHeight)
def gf = GeometryTools.getDefaultFactory()

// Detections grouped by plane, using the same centroid ROI QuPath's DistanceTools uses
def planeKey = { roi -> [roi.getZ(), roi.getT()] }   // List key: a GString key fails map lookup in Groovy
def detsByPlane = detections.groupBy { planeKey(PathObjectTools.getROI(it, true)) }

long t0 = System.currentTimeMillis()
def computed = new HashSet<String>()
annotationClasses.each { pc ->
    String measurementName = "Signed distance to annotation with ${pc} ${unit}"
    def targets = getAnnotationObjects().findAll { it.getPathClass() == pc && it.getROI() != null && it.getROI().isArea() }
    def targetsByPlane = targets.groupBy { planeKey(it.getROI()) }
    targetsByPlane.each { key, planeTargets ->
        def dets = detsByPlane[key]
        if (!dets) return
        def polys = []
        planeTargets.each { ann ->
            def g = ann.getROI().getGeometry()
            if (transform != null) g = transform.transform(g)
            for (int i = 0; i < g.getNumGeometries(); i++) {
                def gi = g.getGeometryN(i)
                if (gi instanceof Polygonal) polys << gi
            }
        }
        if (polys.isEmpty()) return
        def combined = polys.size() == 1 ? polys[0] : GeometryCombiner.combine(polys)
        long tClass = System.currentTimeMillis()
        def precisionModel = combined.getPrecisionModel()   // QuPath snaps query points to this grid; do the same
        def facetDistance = new IndexedFacetDistance(combined)
        def locator = new IndexedPointInAreaLocator(combined)
        locator.locate(new Coordinate(0, 0))   // build index before parallel use (JTS issue #571, as DistanceTools does)

        dets.parallelStream().forEach { p ->
            def roi = PathObjectTools.getROI(p, true)
            def c = new Coordinate(roi.getCentroidX() * pixelWidth, roi.getCentroidY() * pixelHeight)
            precisionModel.makePrecise(c)
            int loc = locator.locate(c)
            double d
            if (loc == Location.BOUNDARY) d = 0.0
            else {
                d = facetDistance.distance(gf.createPoint(c))
                if (loc == Location.INTERIOR) d = -d
            }
            def ml = p.getMeasurementList()
            ml.put(measurementName, d)
            ml.close()
        }
        computed.add(pc.toString())
        println String.format("[02v3]   %s (plane z=%d,t=%d): %d annotation(s), %d vertices, %d detections in %.1f s",
                pc.toString(), key[0], key[1], planeTargets.size(), combined.getNumPoints(), dets.size(),
                (System.currentTimeMillis() - tClass) / 1000.0)
    }
}

// v3: one schema for every image -- absent ensured classes become NaN columns
def filled = ENSURED_CLASSES.findAll { !computed.contains(it) }
filled.each { String cls ->
    String measurementName = "Signed distance to annotation with ${cls} ${unit}"
    detections.parallelStream().forEach { p ->
        def ml = p.getMeasurementList()
        ml.put(measurementName, Double.NaN)
        ml.close()
    }
}
if (!filled.isEmpty()) println "[02v3] No annotations of ${filled}: wrote NaN columns so every image has the same schema"
getCurrentHierarchy().fireObjectMeasurementsChangedEvent(this, detections)
double sec = (System.currentTimeMillis() - t0) / 1000.0

// Verify the ensured columns exist, and that the Islet column was written when islets exist.
def expectedKey = "Signed distance to annotation with ${CLASS_ISLET} ${unit}"
long missingCount = detections.count { Double.isNaN(it.getMeasurementList().get(expectedKey)) }
if (islets.isEmpty())
    println "[02v3] '${expectedKey}' is NaN on all ${detections.size()} detections (no islets)"
else if (missingCount > 0)
    println "[02v3] WARN: '${expectedKey}' missing on ${missingCount} detection(s)."
else
    println "[02v3] Verified: '${expectedKey}' present on all ${detections.size()} detections"
def absentCols = ENSURED_CLASSES.findAll { cls ->
    !detections[0].getMeasurementList().containsKey("Signed distance to annotation with ${cls} ${unit}".toString())
}
if (!absentCols.isEmpty()) throw new IllegalStateException("[02v3] ensured columns missing: ${absentCols}")

println String.format("[02v3] Done in %.1f s", sec)
