import org.locationtech.jts.index.strtree.STRtree
import org.locationtech.jts.geom.prep.PreparedGeometryFactory

// Separate annotations from detections
def hierarchy = getCurrentHierarchy()
def allAnnotations = getAnnotationObjects()
def allDetections = getDetectionObjects()

println "Annotations: ${allAnnotations.size()}, Detections: ${allDetections.size()}"

// Build spatial index (R-tree) over all detections using their centroids/envelopes
def tree = new STRtree()
for (det in allDetections) {
    def env = det.getROI().getGeometry().getEnvelopeInternal()
    tree.insert(env, det)
}
tree.build()
println "Spatial index built"

// Priority order: large structural annotations first, SmallVessel last
def priorityClasses = ['Islet', 'IsletExpanded']
def structuralAnns = allAnnotations.findAll { it.getPathClass()?.toString() in priorityClasses }
def smallVesselAnns = allAnnotations.findAll { it.getPathClass()?.toString() == 'Blood_Vessel' }
def otherAnns = allAnnotations.findAll { 
    def cls = it.getPathClass()?.toString()
    cls != 'Blood_Vessel' && !(cls in priorityClasses)
}

def allSorted = structuralAnns + otherAnns + smallVesselAnns
println "Processing ${structuralAnns.size()} structural, ${smallVesselAnns.size()} Blood_Vessel, ${otherAnns.size()} other annotations"

def pgFactory = new PreparedGeometryFactory()
def assigned = Collections.newSetFromMap(new java.util.concurrent.ConcurrentHashMap())
int count = 0

for (ann in allSorted) {
    def geom = ann.getROI().getGeometry()
    def prepGeom = pgFactory.create(geom)
    def envelope = geom.getEnvelopeInternal()
    
    // Query R-tree for candidates within bounding box
    def candidates = tree.query(envelope)
    
    def children = []
    for (cand in candidates) {
        if (assigned.contains(cand))
            continue
        // Fast containment test using PreparedGeometry (cached spatial index internally)
        def centroid = cand.getROI().getGeometry().getCentroid()
        if (prepGeom.contains(centroid)) {
            children.add(cand)
            assigned.add(cand)
        }
    }
    
    if (!children.isEmpty()) {
        ann.addChildObjects(children)
    }
    
    count++
    if (count % 1000 == 0)
        println "Processed ${count}/${allSorted.size()} annotations"
}

// Fire single hierarchy update at the end
fireHierarchyUpdate()
println "Done. Assigned ${assigned.size()}/${allDetections.size()} detections"