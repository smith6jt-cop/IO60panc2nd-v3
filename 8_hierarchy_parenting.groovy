/*
 * 8_hierarchy_parenting.groovy
 *
 * Parent detection objects to their containing annotation
 *           using an STRtree spatial index and PreparedGeometry.contains on
 *           the detection centroid. Priority order ensures large structural
 *           annotations claim detections before thin bands do:
 *
 *               Islet -> IsletExpanded -> other -> Blood_vessel -> Tissue
 */

evaluate(new File(buildFilePath(PROJECT_BASE_DIR, 'scripts', 'lib', 'common.groovy')))

import org.locationtech.jts.index.strtree.STRtree
import org.locationtech.jts.geom.Coordinate
import org.locationtech.jts.geom.prep.PreparedGeometryFactory
import qupath.lib.objects.PathObjectTools
import qupath.lib.roi.GeometryTools

// true: test the NUCLEUS centroid (cells), the same point 3_signed_distance uses, so a cell parented
// to an Islet always has a negative "Signed distance to annotation with Islet". false: whole-cell
// centroid (what QuPath's own hierarchy resolution uses).
boolean USE_NUCLEUS_CENTROID = true

def hierarchy = getCurrentHierarchy()
def allAnnotations = getAnnotationObjects()
def allDetections  = getDetectionObjects()
println "[08] Annotations: ${allAnnotations.size()}, Detections: ${allDetections.size()}"

if (allDetections.isEmpty()) { println "[08] No detections - nothing to parent."; return }

// Index the test point itself (not the cell envelope): cheaper, and exact for the containment test
def gf = GeometryTools.getDefaultFactory()
def testPoint = new IdentityHashMap()
def tree = new STRtree()
for (det in allDetections) {
    def roi = PathObjectTools.getROI(det, USE_NUCLEUS_CENTROID)
    def pt = gf.createPoint(new Coordinate(roi.getCentroidX(), roi.getCentroidY()))
    testPoint[det] = pt
    tree.insert(pt.getEnvelopeInternal(), det)
}
tree.build()
println "[08] Spatial index built"

// Priority order (first annotation to contain a centroid claims the detection):
//   Islet -> IsletExpanded -> other classes (Nerve, Lymphatic, TissueSample, ...) -> Blood_vessel -> Tissue
// Tissue must be LAST: it contains everything else, so any class processed after it only receives
// detections outside every Tissue annotation. The previous order put Blood_vessel after Tissue, and
// left Tissue vs Nerve/Lymphatic/TissueSample to getAnnotationObjects() order, which comes from the
// hierarchy's spatial-index query and is not defined.
def rank = { ann ->
    switch (ann.getPathClass()?.toString()) {
        case CLASS_ISLET:          return 0
        case CLASS_ISLET_EXPANDED: return 1
        case CLASS_BLOOD_VESSEL:   return 3
        case CLASS_TISSUE:         return 4
        default:                   return 2
    }
}
def allSorted = allAnnotations.findAll { it.getPathClass() != null && it.getROI() != null }.sort { rank(it) }
def nByRank = allSorted.countBy { rank(it) }
println "[08] Processing: ${(nByRank[0] ?: 0) + (nByRank[1] ?: 0)} structural, ${nByRank[2] ?: 0} other, " +
        "${nByRank[3] ?: 0} ${CLASS_BLOOD_VESSEL}, ${nByRank[4] ?: 0} ${CLASS_TISSUE} (last)"

def pgFactory = new PreparedGeometryFactory()
def assigned  = Collections.newSetFromMap(new java.util.concurrent.ConcurrentHashMap())
int count = 0

// ----- Parenting rule: centroid containment -----
for (ann in allSorted) {
    def geom = ann.getROI().getGeometry()
    def prep = pgFactory.create(geom)
    def envelope = geom.getEnvelopeInternal()
    def candidates = tree.query(envelope)

    def children = []
    for (cand in candidates) {
        if (assigned.contains(cand)) continue
        if (prep.contains(testPoint[cand])) {
            children.add(cand)
            assigned.add(cand)
        }
    }
    if (!children.isEmpty()) ann.addChildObjects(children)

    count++
    if (count % 1000 == 0) println "[08] Processed ${count}/${allSorted.size()} annotations"
}

fireHierarchyUpdate()
println "[08] Done. Assigned ${assigned.size()}/${allDetections.size()} detections"
