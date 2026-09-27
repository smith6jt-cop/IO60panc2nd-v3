/*
 * 08_hierarchy_parenting.groovy
 *
 * Workflow step 10. Parent detection objects to their containing annotation
 *           using an STRtree spatial index and PreparedGeometry.contains on
 *           the detection centroid. Priority order ensures large structural
 *           annotations claim detections before thin bands do:
 *
 *               Islet  ->  IsletExpanded  ->  Tissue  ->  Blood_vessel
 *
 *           Note: under the current pipeline, Blood_vessel annotations are
 *           NOT in the hierarchy (stage 07 writes CD34 output only to
 *           GeoJSON). The Blood_vessel bucket below is therefore empty in
 *           practice and is retained only for the case where a user
 *           manually imports a vessel subset from the GeoJSON.
 *
 * Optional tightening: because IsletExpanded bands are thin (~20 um) and a
 * cell body can cross a band boundary with its centroid on the outside, there
 * is a commented section that uses area-fraction-in-annotation (>= 0.5) as
 * the parenting rule instead of centroid containment. Enable it if the
 * centroid rule misplaces peri-islet cells.
 */

evaluate(new File(buildFilePath(PROJECT_BASE_DIR, 'scripts', 'lib', 'common.groovy')))

import org.locationtech.jts.index.strtree.STRtree
import org.locationtech.jts.geom.prep.PreparedGeometryFactory

def hierarchy = getCurrentHierarchy()
def allAnnotations = getAnnotationObjects()
def allDetections  = getDetectionObjects()
println "[08] Annotations: ${allAnnotations.size()}, Detections: ${allDetections.size()}"

if (allDetections.isEmpty()) { println "[08] No detections - nothing to parent."; return }

// Build STRtree over detection envelopes
def tree = new STRtree()
for (det in allDetections) {
    def e = det.getROI().getGeometry().getEnvelopeInternal()
    tree.insert(e, det)
}
tree.build()
println "[08] Spatial index built"

// Priority classes: cores first, then expansions, then tissue, then vessels.
def priorityHigh = [CLASS_ISLET, CLASS_ISLET_EXPANDED]
def priorityLow  = [CLASS_BLOOD_VESSEL]

def structuralAnns = allAnnotations.findAll { it.getPathClass()?.toString() in priorityHigh }
def vesselAnns     = allAnnotations.findAll { it.getPathClass()?.toString() in priorityLow }
def otherAnns      = allAnnotations.findAll {
    def cls = it.getPathClass()?.toString()
    cls != null && !(cls in priorityHigh) && !(cls in priorityLow)
}

def allSorted = structuralAnns + otherAnns + vesselAnns
println "[08] Processing: ${structuralAnns.size()} structural, ${otherAnns.size()} other, ${vesselAnns.size()} ${CLASS_BLOOD_VESSEL}"

def pgFactory = new PreparedGeometryFactory()
def assigned  = Collections.newSetFromMap(new java.util.concurrent.ConcurrentHashMap())
int count = 0

// ----- Parenting rule: centroid containment (default) -----
for (ann in allSorted) {
    def geom = ann.getROI().getGeometry()
    def prep = pgFactory.create(geom)
    def envelope = geom.getEnvelopeInternal()
    def candidates = tree.query(envelope)

    def children = []
    for (cand in candidates) {
        if (assigned.contains(cand)) continue
        def centroid = cand.getROI().getGeometry().getCentroid()
        if (prep.contains(centroid)) {
            children.add(cand)
            assigned.add(cand)
        }
    }
    if (!children.isEmpty()) ann.addChildObjects(children)

    count++
    if (count % 1000 == 0) println "[08] Processed ${count}/${allSorted.size()} annotations"
}

// ----- Optional: area-fraction rule (uncomment to use instead of above) -----
//
// for (ann in allSorted) {
//     def geom = ann.getROI().getGeometry()
//     def prep = pgFactory.create(geom)
//     def envelope = geom.getEnvelopeInternal()
//     def candidates = tree.query(envelope)
//     def children = []
//     for (cand in candidates) {
//         if (assigned.contains(cand)) continue
//         def cg = cand.getROI().getGeometry()
//         double cellArea = cg.getArea()
//         if (cellArea <= 0) continue
//         if (!prep.intersects(cg)) continue
//         double inter = geom.intersection(cg).getArea()
//         if (inter / cellArea >= 0.5) {
//             children.add(cand)
//             assigned.add(cand)
//         }
//     }
//     if (!children.isEmpty()) ann.addChildObjects(children)
// }

fireHierarchyUpdate()
println "[08] Done. Assigned ${assigned.size()}/${allDetections.size()} detections"
