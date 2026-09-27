/*
 * 8_hierarchy_parenting.groovy
 *
 * Parent detection objects to their containing annotation
 *           using an STRtree spatial index and PreparedGeometry.contains on
 *           the detection centroid. Priority order ensures large structural
 *           annotations claim detections before thin bands do:
 *
 *               Islet  ->  IsletExpanded  ->  Tissue
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

// ----- Parenting rule: centroid containment -----
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

fireHierarchyUpdate()
println "[08] Done. Assigned ${assigned.size()}/${allDetections.size()} detections"
