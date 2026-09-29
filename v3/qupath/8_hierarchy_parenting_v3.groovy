/*
 * 8_hierarchy_parenting_v3.groovy   (v3 copy of scripts/8_hierarchy_parenting.groovy; the original is unchanged)
 *
 * Parent detection objects to their containing annotation
 *           using an STRtree spatial index and PreparedGeometry.contains on
 *           the detection centroid. v3 uses an EXPLICIT priority order; the first annotation that
 *           contains a detection's test point claims it:
 *
 *               Islet -> EndocrineCluster -> EndocrineSingle -> Artifact -> pLN -> Tissue   (no IsletExpanded, 2026-09-29)
 *
 *   - endocrine objects first, so a cluster or single lying inside an islet's band keeps its own cells (bands
 *     exist for Islets only; 4_islet_expander_v3.groovy); an endocrine object drawn over an Artifact keeps its
 *     cells too, and carries 'Artifact overlap um^2' (script 4) for exclusion downstream;
 *   - Artifact (folds, debris; drawn at gate 1) and pLN BEFORE the 20 um bands and Tissue, so a band can never
 *     claim cells of a fold or of the lymph node (script 4 also clips bands to Tissue minus Artifact and pLN),
 *     and their cells can be excluded downstream by Parent;
 *   - Tissue LAST: it contains everything else.
 * Every other annotation class (e.g. Nerve, Lymphatic, Blood_vessel from script 7, stale Other: Gal3Norm)
 * is NOT used for parenting; their annotation counts, and how many cells they would have covered, are
 * reported. Detections inside none of the listed classes keep their current parent (reported too).
 */

import org.locationtech.jts.index.strtree.STRtree
import org.locationtech.jts.geom.Coordinate
import org.locationtech.jts.geom.prep.PreparedGeometryFactory
import qupath.lib.objects.PathObjectTools
import qupath.lib.roi.GeometryTools

evaluate(new File(buildFilePath(PROJECT_BASE_DIR, 'scripts', 'lib', 'common.groovy')))

// ==================== PARAMETERS ====================
// no IsletExpanded: peri-islet regions dropped (user, 2026-09-29); distances (script 3) carry that relation
List<String> PARENT_ORDER = ['Islet', 'EndocrineCluster', 'EndocrineSingle', 'Artifact', 'pLN', 'Tissue']
// true: test the NUCLEUS centroid (cells), the same point 3_signed_distance uses, so a cell parented
// to an Islet always has a negative "Signed distance to annotation with Islet". false: whole-cell
// centroid (what QuPath's own hierarchy resolution uses).
boolean USE_NUCLEUS_CENTROID = true
// ====================================================

def hierarchy = getCurrentHierarchy()
def allAnnotations = getAnnotationObjects()
def allDetections  = getDetectionObjects()
println "[08v3] Annotations: ${allAnnotations.size()}, Detections: ${allDetections.size()}"

// No early `return` (run_v3 appends its provenance footer); nothing to parent is a failure.
if (allDetections.isEmpty()) throw new IllegalStateException("[08v3] No detections - nothing to parent.")

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
println "[08v3] Spatial index built"

def clsOf = { ann -> ann.getPathClass()?.toString() }
def rank = { ann -> PARENT_ORDER.indexOf(clsOf(ann)) }   // -1 = not used for parenting
def usable = allAnnotations.findAll { it.getPathClass() != null && it.getROI() != null && it.getROI().isArea() }
def allSorted = usable.findAll { rank(it) >= 0 }.sort { a, b -> rank(a) <=> rank(b) }   // stable within a class
def sortedSet = Collections.newSetFromMap(new IdentityHashMap())
sortedSet.addAll(allSorted)
def unused = allAnnotations.findAll { !sortedSet.contains(it) }
def nByClass = new LinkedHashMap()
PARENT_ORDER.each { c -> nByClass[c] = allSorted.count { clsOf(it) == c } }
println "[08v3] Parenting order and annotation counts: ${nByClass}"
def unusedByClass = new TreeMap(unused.countBy { clsOf(it) ?: 'Unclassified' })
if (!unused.isEmpty()) println "[08v3] Not used for parenting (by design): ${unusedByClass}"

def pgFactory = new PreparedGeometryFactory()
def assigned  = Collections.newSetFromMap(new java.util.concurrent.ConcurrentHashMap())
def claimedByClass = new LinkedHashMap()
PARENT_ORDER.each { claimedByClass[it] = 0 }
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
    if (!children.isEmpty()) {
        ann.addChildObjects(children)
        claimedByClass[clsOf(ann)] += children.size()
    }

    count++
    if (count % 1000 == 0) println "[08v3] Processed ${count}/${allSorted.size()} annotations"
}

// Report what the unused classes would have covered (no parenting change)
def coveredByUnused = new TreeMap()
unused.findAll { it.getROI() != null && it.getROI().isArea() }.groupBy { clsOf(it) ?: 'Unclassified' }.each { c, anns ->
    long n = 0
    anns.each { ann ->
        def prep = pgFactory.create(ann.getROI().getGeometry())
        n += tree.query(ann.getROI().getGeometry().getEnvelopeInternal()).count { prep.contains(testPoint[it]) }
    }
    coveredByUnused[c] = n
}

fireHierarchyUpdate()
println "[08v3] Detections claimed per class: ${claimedByClass}"
if (!coveredByUnused.isEmpty()) println "[08v3] Detections inside classes not used for parenting (kept their assigned parent): ${coveredByUnused}"
println "[08v3] Done. Assigned ${assigned.size()}/${allDetections.size()} detections; " +
        "${allDetections.size() - assigned.size()} lie in none of ${PARENT_ORDER} and keep their current parent"
