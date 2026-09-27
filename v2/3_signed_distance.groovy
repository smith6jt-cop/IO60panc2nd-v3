/*
 * 3_signed_distance.groovy
 *
 * Workflow step 4. Scripted equivalent of QuPath UI:
 *   Analyze -> Spatial analysis -> Signed Distance to Annotations 2D
 *
 * Adds "Signed distance to annotation with <class> um" measurements to every
 * detection, where <class> is each non-ignored annotation classification
 * present in the hierarchy. Negative inside, positive outside.
 *
 * MUST run BEFORE 4_islet_expander.groovy creates IsletExpanded
 * annotations. The classifications present at this point are Tissue and
 * Islet, so the relevant output column for downstream classifier work is
 *   "Signed distance to annotation with Islet um"
 *
 * The "Signed distance to annotation with Tissue um" column will also be
 * written; it is harmless and may be useful.
 *
 * The QuPath API does not provide a built-in way to restrict this command
 * to a single class. If filtering is needed, the alternatives are:
 *   - run this stage now (only Tissue and Islet exist), then ignore the
 *     Tissue distance column downstream; OR
 *   - temporarily set Tissue's PathClass to a class with the *Ignore*
 *     marker before running, then restore.
 *
 * Implementation: replicates DistanceTools.detectionToAnnotationDistancesSigned
 * per-class, but calls the lower-level DistanceTools.centroidToBoundsDistance2D
 * against disposable, simplified copies of each class's annotation geometry
 * (see note below) instead of the raw hierarchy objects. splitClassNames is
 * effectively false, since pancreas classes are simple names with no
 * colon-separated subcomponents.
 */

evaluate(new File(buildFilePath(PROJECT_BASE_DIR, 'scripts', 'lib', 'common.groovy')))

import qupath.lib.analysis.DistanceTools
import qupath.lib.objects.PathObjects
import qupath.lib.objects.classes.PathClassTools
import qupath.lib.roi.GeometryTools
import org.locationtech.jts.simplify.TopologyPreservingSimplifier

// DistanceTools' boundary-distance step (JTS DistanceToPoint) is brute
// force over every vertex of the target geometry, with no spatial index -
// only the inside/outside test is indexed. Islet/Tissue annotations come
// straight out of createAnnotationsFromPixelClassifier (step 01) and are
// never simplified before this stage, so they carry near pixel-perfect
// vertex counts. That is what made this take 34 min on 2.87M detections:
// runtime is O(detections x target_vertices), not O(detections) as the
// "faster than the GUI" framing assumed. We simplify disposable copies of
// the annotation geometries (never touching the real hierarchy objects,
// so step 04's expansion still sees full-precision boundaries) purely for
// this distance query.
SIMPLIFY_TOL_UM = 2.0

def imageData = getCurrentImageData()
if (imageData == null) { println "[02] No image open."; return }

def detections = getDetectionObjects()
if (detections.isEmpty()) {
    println "[02] No detections."
    return
}

def islets = getAnnotationObjects().findAll {
    it.getPathClass()?.toString() == CLASS_ISLET && it.getROI() != null
}
if (islets.isEmpty()) {
    println "[02] No ${CLASS_ISLET} annotations. Run 00_panc_objects.groovy first."
    return
}

def alreadyExpanded = getAnnotationObjects().any {
    it.getPathClass()?.toString() == CLASS_ISLET_EXPANDED
}
if (alreadyExpanded) {
    println "[02] WARN: ${CLASS_ISLET_EXPANDED} annotations exist."
    println "[02]       Distances will also be computed to that class."
    println "[02]       Run this BEFORE step 5 to avoid that."
}

def annotationClasses = getAnnotationObjects()
    .collect { it.getPathClass() }
    .findAll { it != null && it.isValid() && !PathClassTools.isIgnoredClass(it) }
    .toSet()
println "[02] Computing signed distances. Classes present: ${annotationClasses}"
println "[02] Detections: ${detections.size()}"
System.out.flush()

def cal = imageData.getServer().getPixelCalibration()
double pixelWidth  = cal.getPixelWidth().doubleValue()
double pixelHeight = cal.getPixelHeight().doubleValue()
double pxSize      = cal.getAveragedPixelSize().doubleValue()
double simplifyTolPx = SIMPLIFY_TOL_UM / pxSize
String unit = cal.getPixelWidthUnit()
println "[02] Simplify tolerance: ${SIMPLIFY_TOL_UM} um (${String.format('%.3f', simplifyTolPx)} px)"
System.out.flush()

long t0 = System.currentTimeMillis()
annotationClasses.each { pc ->
    def targets = getAnnotationObjects().findAll { it.getPathClass() == pc && it.getROI() != null }
    if (targets.isEmpty()) return
    long totalVerticesBefore = targets.sum { it.getROI().getGeometry().getNumPoints() }
    println "[02]   ${pc}: ${targets.size()} annotation(s), ${totalVerticesBefore} vertices - simplifying..."
    System.out.flush()
    long tClass = System.currentTimeMillis()
    def simplifiedTargets = targets.collect { ann ->
        def roi = ann.getROI()
        def raw = roi.getGeometry()
        def geom
        try {
            // isValid() is a cheap single pass; skip the buffer(0) topology
            // repair (expensive on large pixel-classifier polygons) unless
            // the geometry actually needs it.
            def fixed = raw.isValid() ? raw : raw.buffer(0)
            geom = TopologyPreservingSimplifier.simplify(fixed, simplifyTolPx)
        } catch (Throwable t) {
            geom = raw
        }
        def simpRoi = GeometryTools.geometryToROI(geom, roi.getImagePlane())
        return PathObjects.createAnnotationObject(simpRoi, pc)
    }
    long totalVerticesAfter  = simplifiedTargets.sum { it.getROI().getGeometry().getNumPoints() }
    double simplifySec = (System.currentTimeMillis() - tClass) / 1000.0
    println String.format("[02]   ${pc}: vertices %d -> %d (%.1f s), computing distances...", totalVerticesBefore, totalVerticesAfter, simplifySec)
    System.out.flush()
    String measurementName = "Signed distance to annotation with ${pc} ${unit}"
    long tDist = System.currentTimeMillis()
    DistanceTools.centroidToBoundsDistance2D(detections, simplifiedTargets, pixelWidth, pixelHeight, measurementName, true)
    double distSec = (System.currentTimeMillis() - tDist) / 1000.0
    println String.format("[02]   ${pc}: distances done (%.1f s)", distSec)
    System.out.flush()
}
fireHierarchyUpdate()
double sec = (System.currentTimeMillis() - t0) / 1000.0

// Verify the expected Islet column was written.
def expectedKey = "Signed distance to annotation with ${CLASS_ISLET} ${unit}"
def sample = detections[0]
double check = getMeas(sample, expectedKey)
if (Double.isNaN(check)) {
    println "[02] WARN: expected measurement '${expectedKey}' not found on first detection."
    println "[02]       Check the measurement table for the actual column name."
} else {
    println "[02] Verified: '${expectedKey}' = ${check} on first detection"
}

println String.format("[02] Done in %.1f s", sec)