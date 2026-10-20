/*
 * islet_cells_1c_synthetic.groovy -- test fixture for 1c_islet_objects_to_cells_v3.groovy: known geometric cases.
 *
 * DESTRUCTIVE on the image it runs on: removes every detection and every Islet/EndocrineCluster/EndocrineSingle,
 * then draws synthetic cells (squares on InstanSeg's 0.985 px label grid: 0.5 um model pixels on 0.5078 um image
 * pixels) and synthetic pixel objects. It refuses to run unless PRIVATE_PROJECT=<project directory> names the
 * current project and that project is not under /home/smith6jt/IO60panc2nd. Run with `QuPath script -s` (saves).
 * Writes EXPECT=<path> (JSON) with the answer check_islet_cells_1c.py --expect compares for each scenario box.
 *
 * Scenarios (lp = label pixels of 0.985 px; cells are 15 x 15 lp, ~56 um^2):
 *   A halo_single   3x3 touching cells, a disc covering the centre cell and <20% of each neighbour
 *                   (mean caliper 11 um: script 1 calls it a cluster)            -> 1 EndocrineSingle
 *   B halo_only     a small disc on the junction of 4 cells (13% of each)        -> dropped, no object
 *   C merge         2 pixel objects (route 1 + bridge, route 2), kept cells 1 lp apart -> 1 object of 4 cells,
 *                   n_pixel_objects 2, route 3, bridge 1
 *   D split         1 pixel object over 2 cells 2 lp (1.97 px) apart            -> 2 singles
 *   E corner        2 cells touching at one corner                              -> 1 object of 2 cells
 *   F diagonal      2 cells 1 lp apart diagonally (1.39 px)                      -> 2 singles at TOUCH_GAP_PX 1
 *   G hole          ring-shaped pixel object over 8 of 3x3 cells                -> 1 object of 8, centre cell in its hole
 *   H swallow       16-cell ring + a separate inner cell 18 lp from it (2 objects) -> 1 Islet of 17, 2 pixel objects
 *   I cover         a pixel object over 60% of one cell and 40% of another      -> 1 single, min_cell_cover 0.6
 *   J islet         10x10 touching cells (75 um)                                -> 1 Islet of 100
 *   K border        1 cell 6 px (3 um) from the Tissue edge                     -> 1 single, Tissue border flag 1
 *   L fragments     1 cell whose outline has 3 parts: body, a fragment 1 lp away, one 7 lp away
 *                   -> 1 single holding the whole cell: body + near fragment joined, far fragment a 2nd part
 *   M protrusion    2 member cells 1 lp apart; a NON-member cell pokes a 1 lp wide arm into their gap
 *                   -> 1 object of 2 cells that takes none of the non-member (the gap fill is trimmed)
 */
import qupath.lib.objects.PathObjects
import qupath.lib.objects.classes.PathClass
import qupath.lib.roi.GeometryTools
import org.locationtech.jts.geom.Coordinate
import org.locationtech.jts.geom.Geometry

def A = [PRIVATE_PROJECT: '', EXPECT: '']
((binding.hasVariable('args') && args != null) ? (args as List) : []).each { String a ->
    int eq = a.indexOf('=')
    if (eq <= 0 || !A.containsKey(a.substring(0, eq))) throw new IllegalArgumentException("bad argument '${a}'; known ${A.keySet()}")
    A[a.substring(0, eq)] = a.substring(eq + 1)
}
def projDir = getProject()?.getPath()?.getParent()?.toFile()?.getCanonicalPath()
if (!A.PRIVATE_PROJECT || projDir == null || new File(A.PRIVATE_PROJECT).getCanonicalPath() != projDir ||
        projDir.startsWith('/home/smith6jt/IO60panc2nd'))
    throw new IllegalStateException("refusing: PRIVATE_PROJECT must name this project's directory (${projDir}), outside the real project")
if (!A.EXPECT) throw new IllegalArgumentException('EXPECT=<json path> is required')

def hier = getCurrentHierarchy()
def gf = GeometryTools.getDefaultFactory()
def plane = qupath.lib.regions.ImagePlane.getDefaultPlane()
double s = 0.5d / 0.5078d            // one InstanSeg label pixel in image pixels
double um = 0.5078d
def ENDO = ['Islet', 'EndocrineCluster', 'EndocrineSingle']
removeObjects(new ArrayList(hier.getDetectionObjects()), false)
removeObjects(getAnnotationObjects().findAll { it.getPathClass()?.toString() in ENDO }, false)

def box = { double ox, double oy, double x0, double y0, double x1, double y1 ->   // lp -> image px
    gf.toGeometry(new org.locationtech.jts.geom.Envelope(ox + x0 * s, ox + x1 * s, oy + y0 * s, oy + y1 * s))
}
def cells = []
def cellGeom = { Geometry g, Geometry n ->
    cells << PathObjects.createCellObject(GeometryTools.geometryToROI(g, plane), GeometryTools.geometryToROI(n, plane))
}
def cell = { double ox, double oy, double x0, double y0, double w = 15, double h = 15 ->
    def g = box(ox, oy, x0, y0, x0 + w, y0 + h)
    def n = box(ox, oy, x0 + w / 4, y0 + h / 4, x0 + 3 * w / 4, y0 + 3 * h / 4)
    cells << PathObjects.createCellObject(GeometryTools.geometryToROI(g, plane), GeometryTools.geometryToROI(n, plane))
}
def pix = []
def pixObj = { Geometry g, String c, int route, int bridge ->
    def o = PathObjects.createAnnotationObject(GeometryTools.geometryToROI(g, plane), PathClass.fromString(c))
    o.getMeasurementList().put('Endocrine route', route)
    o.getMeasurementList().put('Bridge flag', bridge)
    o.getMeasurementList().close()
    pix << o
}
def circle = { double ox, double oy, double cx, double cy, double r -> gf.createPoint(new Coordinate(ox + cx * s, oy + cy * s)).buffer(r * s, 32) }
def expect = []
def scen = { String name, double ox, double oy, double wLp, double hLp, Map e ->
    expect << ([name: name, box: [ox - 10, oy - 10, ox + wLp * s + 10, oy + hLp * s + 10]] + e)
}

// A halo_single
double ox = 200, oy = 200
for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) cell(ox, oy, 15 * i, 15 * j)
pixObj(circle(ox, oy, 22.5, 22.5, 11), 'EndocrineCluster', 2, 0)
scen('A_halo_single', ox, oy, 45, 45, [n_objects: 1, classes: ['EndocrineSingle'], n_cells: [1], n_pixel_objects: [1], route: [2],
                                        bridge: [0], border: [0], hole_cells: 0, dropped: 0])
// B halo_only
ox = 400; oy = 200
for (int i = 0; i < 2; i++) for (int j = 0; j < 2; j++) cell(ox, oy, 15 * i, 15 * j)
pixObj(circle(ox, oy, 15, 15, 6), 'EndocrineSingle', 2, 0)
scen('B_halo_only', ox, oy, 30, 30, [n_objects: 0, hole_cells: 0, dropped: 1])
// C merge (kept cells 1 lp apart; routes 1 + 2; bridge on one)
ox = 600; oy = 200
[0, 15, 31, 46].each { x -> cell(ox, oy, x, 0) }
pixObj(box(ox, oy, 0, 0, 30, 15), 'EndocrineCluster', 1, 1)
pixObj(box(ox, oy, 31, 0, 61, 15), 'EndocrineCluster', 2, 0)
scen('C_merge', ox, oy, 61, 15, [n_objects: 1, classes: ['EndocrineCluster'], n_cells: [4], n_pixel_objects: [2], route: [3],
                                  bridge: [1], border: [0], hole_cells: 0, dropped: 0])
// D split (2 lp apart)
ox = 800; oy = 200
cell(ox, oy, 0, 0); cell(ox, oy, 17, 0)
pixObj(box(ox, oy, 0, 0, 32, 15), 'EndocrineCluster', 3, 0)
scen('D_split', ox, oy, 32, 15, [n_objects: 2, classes: ['EndocrineSingle', 'EndocrineSingle'], n_cells: [1, 1], n_pixel_objects: [1, 1],
                                  route: [3, 3], hole_cells: 0, dropped: 0, split: 1])
// E corner touch
ox = 1000; oy = 200
cell(ox, oy, 0, 0); cell(ox, oy, 15, 15)
pixObj(box(ox, oy, 0, 0, 30, 30), 'EndocrineCluster', 2, 0)
scen('E_corner', ox, oy, 30, 30, [n_objects: 1, classes: ['EndocrineCluster'], n_cells: [2], n_pixel_objects: [1], hole_cells: 0, dropped: 0])
// F diagonal 1 lp gap (1.39 px)
ox = 1200; oy = 200
cell(ox, oy, 0, 0); cell(ox, oy, 16, 16)
pixObj(box(ox, oy, 0, 0, 31, 31), 'EndocrineCluster', 2, 0)
scen('F_diagonal', ox, oy, 31, 31, [n_objects: 2, classes: ['EndocrineSingle', 'EndocrineSingle'], n_cells: [1, 1], hole_cells: 0, dropped: 0, split: 1])
// G ring-shaped pixel object: the centre cell ends up in the new object's filled hole
ox = 200; oy = 500
for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) cell(ox, oy, 15 * i, 15 * j)
pixObj(box(ox, oy, 0, 0, 45, 45).difference(box(ox, oy, 15, 15, 30, 30)), 'EndocrineCluster', 2, 0)
scen('G_hole', ox, oy, 45, 45, [n_objects: 1, classes: ['EndocrineCluster'], n_cells: [8], n_pixel_objects: [1], hole_cells: 1, dropped: 0])
// H ring of 16 cells + a separate inner cell (its own pixel object) inside the ring's hole
ox = 500; oy = 500
for (int i = 0; i < 5; i++) for (int j = 0; j < 5; j++) if (!(i in 1..3 && j in 1..3)) cell(ox, oy, 15 * i, 15 * j)
cell(ox, oy, 33, 33, 9, 9)
pixObj(box(ox, oy, 0, 0, 75, 75).difference(box(ox, oy, 15, 15, 60, 60)), 'Islet', 2, 0)
pixObj(box(ox, oy, 30, 30, 45, 45), 'EndocrineSingle', 1, 0)
scen('H_swallow', ox, oy, 75, 75, [n_objects: 1, classes: ['Islet'], n_cells: [17], n_pixel_objects: [2], route: [3], hole_cells: 0, dropped: 0])
// I coverage threshold: 60% of X, 40% of Y
ox = 900; oy = 500
cell(ox, oy, 0, 0); cell(ox, oy, 20, 0)
pixObj(box(ox, oy, 6, 0, 26, 15), 'EndocrineCluster', 2, 0)
scen('I_cover', ox, oy, 35, 15, [n_objects: 1, classes: ['EndocrineSingle'], n_cells: [1], min_cell_cover: 0.6, hole_cells: 0, dropped: 0])
// J islet of 10x10 cells
ox = 1100; oy = 500
for (int i = 0; i < 10; i++) for (int j = 0; j < 10; j++) cell(ox, oy, 15 * i, 15 * j)
pixObj(box(ox, oy, 0, 0, 150, 150), 'Islet', 3, 0)
scen('J_islet', ox, oy, 150, 150, [n_objects: 1, classes: ['Islet'], n_cells: [100], n_pixel_objects: [1], route: [3], border: [0], hole_cells: 0, dropped: 0])
// K near the Tissue edge (Tissue starts at x = 16 px; the cell at 22 px)
ox = 22; oy = 1000
cell(ox, oy, 0, 0)
pixObj(box(ox, oy, 0, 0, 15, 15), 'EndocrineSingle', 2, 0)
scen('K_border', ox, oy, 15, 15, [n_objects: 1, classes: ['EndocrineSingle'], n_cells: [1], border: [1], hole_cells: 0, dropped: 0])

// L one cell with a 3-part outline (InstanSeg label with detached pixels)
ox = 1300; oy = 800
cellGeom(gf.createMultiPolygon([box(ox, oy, 0, 0, 15, 15), box(ox, oy, 16, 5, 18, 8), box(ox, oy, 22, 5, 24, 8)] as org.locationtech.jts.geom.Polygon[]),
        box(ox, oy, 4, 4, 11, 11))
pixObj(box(ox, oy, 0, 0, 25, 15), 'EndocrineCluster', 2, 0)
scen('L_fragments', ox, oy, 25, 15, [n_objects: 1, classes: ['EndocrineSingle'], n_cells: [1], parts: [2], hole_cells: 0, dropped: 0])
// M a non-member arm in the 1 lp gap between two members
ox = 1500; oy = 800
cell(ox, oy, 0, 0); cell(ox, oy, 16, 0)
cellGeom(box(ox, oy, 8, 16, 23, 24).union(box(ox, oy, 15, 4, 16, 16)), box(ox, oy, 12, 18, 19, 22))
pixObj(box(ox, oy, 0, 0, 31, 15), 'EndocrineCluster', 2, 0)
scen('M_protrusion', ox, oy, 31, 24, [n_objects: 1, classes: ['EndocrineCluster'], n_cells: [2], parts: [1], hole_cells: 0, dropped: 0])

addObjects(cells)
addObjects(pix)
fireHierarchyUpdate()
new File(A.EXPECT).text = qupath.lib.io.GsonTools.getInstance(true).toJson([scenarios: expect, n_cells: cells.size(), n_pixel_objects: pix.size()])
println "[synthetic1c] ${cells.size()} cells, ${pix.size()} pixel objects, ${expect.size()} scenarios -> ${A.EXPECT}"
