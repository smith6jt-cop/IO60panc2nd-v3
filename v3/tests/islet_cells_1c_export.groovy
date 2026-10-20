/*
 * islet_cells_1c_export.groovy -- test helper for 1c_islet_objects_to_cells_v3.groovy (READ-ONLY: never saves).
 *
 * Exports what check_islet_cells_1c.py needs to re-check a 1c run independently (shapely, not JTS):
 *   <OUT>.regions.geojson   Tissue / pLN / Artifact, with QuPath's own GeoJSON exporter: the input of the anchors'
 *                           annotation hash (islet_anchors.annotation_hash recomputes it from this file)
 *   <OUT>.islets.geojson    Islet / EndocrineCluster / EndocrineSingle with their measurements
 *   <OUT>.cells.geojson     cells whose bounding box lies within NEAR_PX of an islet-class object's box
 *                           (CELLS=near, default), every cell (CELLS=all) or none; no measurements
 *   <OUT>.cell_ids.txt      every cell's ID, sorted (1c must not change the cells)
 *   <OUT>.info.json         annotation counts per class, cell count, parents of the islet-class objects,
 *                           detections directly under the root
 *
 *   QuPath script -p <project> -i <image> -a OUT=/path/prefix [-a CELLS=near|all|none] [-a NEAR_PX=40] \
 *       scripts/v3/tests/islet_cells_1c_export.groovy
 */
import qupath.lib.io.PathIO
import org.locationtech.jts.geom.Envelope
import org.locationtech.jts.index.strtree.STRtree

def A = [OUT: '', CELLS: 'near', NEAR_PX: '40']
((binding.hasVariable('args') && args != null) ? (args as List) : []).each { String a ->
    int eq = a.indexOf('=')
    if (eq <= 0 || !A.containsKey(a.substring(0, eq))) throw new IllegalArgumentException("bad argument '${a}'; known ${A.keySet()}")
    A[a.substring(0, eq)] = a.substring(eq + 1)
}
if (!A.OUT) throw new IllegalArgumentException('OUT=<path prefix> is required')
if (!(A.CELLS in ['near', 'all', 'none'])) throw new IllegalArgumentException("CELLS must be near|all|none")
new File(A.OUT).getAbsoluteFile().getParentFile().mkdirs()
def ENDO = ['Islet', 'EndocrineCluster', 'EndocrineSingle']
def REGIONS = ['Tissue', 'pLN', 'Artifact']
def cls = { o -> o?.getPathClass()?.toString() }
def hier = getCurrentHierarchy()
def anns = getAnnotationObjects()
def write = { String path, Collection objs, boolean measurements ->
    new File(path).withOutputStream { os ->
        if (measurements) PathIO.exportObjectsAsGeoJSON(os, objs, PathIO.GeoJsonExportOptions.FEATURE_COLLECTION)
        else PathIO.exportObjectsAsGeoJSON(os, objs, PathIO.GeoJsonExportOptions.FEATURE_COLLECTION,
                PathIO.GeoJsonExportOptions.EXCLUDE_MEASUREMENTS)
    }
}
write("${A.OUT}.regions.geojson", anns.findAll { cls(it) in REGIONS }, true)
def endo = anns.findAll { cls(it) in ENDO }
write("${A.OUT}.islets.geojson", endo, true)
def cells = new ArrayList(hier.getCellObjects())
def picked = []
if (A.CELLS == 'all') picked = cells
else if (A.CELLS == 'near' && !endo.isEmpty()) {
    double near = A.NEAR_PX as double
    def t = new STRtree()
    endo.each { o -> def e = new Envelope(o.getROI().getGeometry().getEnvelopeInternal()); e.expandBy(near); t.insert(e, o) }
    t.build()
    picked = cells.findAll { c -> def r = c.getROI()
        !t.query(new Envelope(r.getBoundsX(), r.getBoundsX() + r.getBoundsWidth(), r.getBoundsY(), r.getBoundsY() + r.getBoundsHeight())).isEmpty() }
}
write("${A.OUT}.cells.geojson", picked, false)
new File("${A.OUT}.cell_ids.txt").text = cells.collect { it.getID().toString() }.sort().join('\n') + '\n'
def parentOf = { o -> def p = o.getParent(); p == null ? 'NONE' : (p.isRootObject() ? 'ROOT' : (cls(p) ?: 'Unclassified')) }
def info = [annotations_by_class: new TreeMap(anns.countBy { cls(it) ?: 'Unclassified' }), n_cells: cells.size(),
            n_cells_exported: picked.size(), islet_parents: new TreeMap(endo.countBy { parentOf(it) }),
            detections_under_root: hier.getRootObject().getChildObjects().count { it.isDetection() }]
new File("${A.OUT}.info.json").text = qupath.lib.io.GsonTools.getInstance(true).toJson(info)
println "[export1c] ${A.OUT}: ${info}"
