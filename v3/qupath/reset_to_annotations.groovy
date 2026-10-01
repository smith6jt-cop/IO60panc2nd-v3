/*
 * reset_to_annotations.groovy -- start an image's data file over before a new cells pass.
 *
 * The Tissue (and pLN) annotations were redrawn on 2026-09-30, so the cells and the pixel islet
 * objects of the 2026-09-28/29 run no longer describe the image. This step keeps the annotations
 * that are the maintainer's own (KEEP classes, with their object IDs) and everything the data file
 * holds about the image itself (image type, display settings), and removes the rest:
 *
 *   strip-all   every detection and every annotation outside KEEP        (Tissue or pLN changed)
 *   keep-cells  every annotation outside KEEP; detections stay            (outlines unchanged, so
 *               the existing cells are still the cells of this Tissue)
 *
 * It refuses to save unless exactly one Tissue annotation remains (or, with a 4th argument "1+", at
 * least one: an IO60pLN slide has one Tissue outline per node section) and, in keep-cells mode, the
 * detection count is unchanged. Success is the sentinel JSON, not QuPath's exit code.
 *
 *   QuPath script -p <private project copy> -i "<image>" \
 *       -a <strip-all|keep-cells> -a <sentinel.json> [-a "Tissue|pLN" [-a 1+]] reset_to_annotations.groovy
 */

import qupath.lib.io.GsonTools   // groovy-json is not bundled with QuPath 0.7

def a = (binding.hasVariable('args') && args != null) ? (args as List).collect { it.toString() } : []
if (a.size() < 2)
    throw new IllegalArgumentException('usage: reset_to_annotations.groovy <strip-all|keep-cells> <sentinel.json> [KEEP classes, |-separated]')
String mode = a[0]
def sentinelFile = new File(a[1])
Set<String> KEEP = ((a.size() > 2 && a[2]) ? a[2].split('\\|').collect { it.trim() } : ['Tissue', 'pLN']) as Set
if (!(mode in ['strip-all', 'keep-cells']))
    throw new IllegalArgumentException("unknown mode '${mode}'")
String tissueRule = (a.size() > 3 && a[3]) ? a[3].trim() : '1'
if (!(tissueRule in ['1', '1+']))
    throw new IllegalArgumentException("unknown Tissue rule '${tissueRule}' (1 or 1+)")

def imageData = getCurrentImageData()
def hierarchy = imageData.getHierarchy()
def entry = getProjectEntry()
def clsOf = { o -> o.getPathClass()?.toString() ?: 'Unclassified' }
def countBy = { objs -> new TreeMap(objs.countBy { clsOf(it) }) }

def annBefore = hierarchy.getAnnotationObjects()
int detBefore = hierarchy.getDetectionObjects().size()
def kept = annBefore.findAll { clsOf(it) in KEEP }
def drop = annBefore.findAll { !(clsOf(it) in KEEP) }
def keptIds = kept.collect { it.getID().toString() }.sort()

def result = [image: entry?.getImageName(), entry_id: entry?.getID(), mode: mode, keep: KEEP.sort(),
              before: [annotations: countBy(annBefore), detections: detBefore],
              kept_ids: keptIds, ok: false, saved: false]

int nTissue = kept.count { clsOf(it) == 'Tissue' }
result.tissue_rule = tissueRule
if (tissueRule == '1' ? nTissue != 1 : nTissue < 1) {
    result.message = "expected ${tissueRule == '1' ? 'exactly one' : 'at least one'} Tissue annotation, found ${nTissue}".toString()
} else {
    if (mode == 'strip-all') removeDetections()
    hierarchy.removeObjects(drop, true)      // keepChildren: detections under a removed object move up
    fireHierarchyUpdate()

    def annAfter = hierarchy.getAnnotationObjects()
    int detAfter = hierarchy.getDetectionObjects().size()
    result.after = [annotations: countBy(annAfter), detections: detAfter]
    boolean idsSame = annAfter.collect { it.getID().toString() }.sort() == keptIds
    boolean detOk = (mode == 'strip-all') ? detAfter == 0 : detAfter == detBefore
    result.ids_same = idsSame
    result.detections_ok = detOk
    if (idsSame && detOk) {
        entry.saveImageData(imageData)
        result.saved = true
        result.ok = true
    } else {
        result.message = 'post-check failed; nothing saved'
    }
}
result.finished_at = java.time.OffsetDateTime.now().toString()
sentinelFile.getParentFile()?.mkdirs()
sentinelFile.text = GsonTools.getInstance(true).toJson(result)
println "[reset] ${result.image}: mode ${mode}, ok=${result.ok}, before ${result.before}, after ${result.after}"
