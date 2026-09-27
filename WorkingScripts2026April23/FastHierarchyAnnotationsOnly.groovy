import org.locationtech.jts.index.strtree.STRtree
import org.locationtech.jts.geom.prep.PreparedGeometryFactory

def hierarchy = getCurrentHierarchy()
def allAnnotations = getAnnotationObjects()
println "Total annotations: ${allAnnotations.size()}"

// Define parent classes (large structural annotations)
def parentClasses = ['Follicle', 'PALS', 'LargeVessel', 'Trabeculae'] as Set

def parents = allAnnotations.findAll { parentClasses.contains(it.getPathClass()?.toString()) }
def children = allAnnotations.findAll { !parentClasses.contains(it.getPathClass()?.toString()) }
println "Parents: ${parents.size()}, Children: ${children.size()}"

// Build R-tree over child annotations
def tree = new STRtree()
for (child in children) {
    def env = child.getROI().getGeometry().getEnvelopeInternal()
    tree.insert(env, child)
}
tree.build()

def pgFactory = new PreparedGeometryFactory()
def assigned = new HashSet()
int count = 0

for (parent in parents) {
    def geom = parent.getROI().getGeometry()
    def prepGeom = pgFactory.create(geom)
    def candidates = tree.query(geom.getEnvelopeInternal())
    
    def matched = []
    for (cand in candidates) {
        if (assigned.contains(cand))
            continue
        def centroid = cand.getROI().getGeometry().getCentroid()
        if (prepGeom.contains(centroid)) {
            matched.add(cand)
            assigned.add(cand)
        }
    }
    
    if (!matched.isEmpty())
        parent.addChildObjects(matched)
    
    count++
    if (count % 500 == 0)
        println "Processed ${count}/${parents.size()} parents"
}

fireHierarchyUpdate()
println "Done. Nested ${assigned.size()}/${children.size()} child annotations into ${parents.size()} parents"