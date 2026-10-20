/*
 * 1c_islet_objects_to_cells_v3.groovy   (v3 step 1c: pixel islet objects -> cell-based objects)
 *
 * Runs per image AFTER 1_panc_objects_v3.groovy (the pixel objects) and the cells pass (InstanSeg cell
 * detections), and BEFORE 4_islet_expander_v3.groovy. Script 1's objects are thresholded, blurred pixel
 * regions: a bright cell's spill-over halo covers parts of its neighbours (one GCG-bright cell becomes a
 * 5-7 cell "EndocrineCluster" with a mean caliper > 10 um), inflates islet edges and joins nearby islets.
 * This step resolves every object to whole InstanSeg cell outlines (scripts/v3/ISLETS.md, section 2b):
 *
 *   1. coverage: for every Islet / EndocrineCluster / EndocrineSingle annotation (the "pixel objects"), the
 *      cells whose outline intersects it (STRtree of the pixel objects, queried with each cell's bounding
 *      box: never all cells x all objects). coverage = area(cell n object) / area(cell). A cell is KEPT when
 *      its coverage by one pixel object is >= COVER_MIN.
 *   2. components, across ALL pixel objects: kept cells whose outlines are at most TOUCH_GAP_PX full-resolution
 *      pixels apart touch (InstanSeg's label grid is 0.985 px here: one background label pixel between two
 *      cells is a 0.985 px gap, 1.39 px diagonally). Each connected component is one new object: the union of
 *      its member cells' WHOLE outlines plus, for each touching pair that does not share an edge, the gap
 *      between them (the region within their distance + 0.05 px of both cells, minus every cell that is not a
 *      member, so a fill never takes part of another cell), with polygon holes filled. Parts of one component
 *      still within TOUCH_GAP_PX of each other (e.g. a cell's own detached fragment) are joined the same way;
 *      parts farther apart keep the object multi-part (rare; logged). A component inside another's filled hole
 *      joins it.
 *   3. class: SINGLE_RULE 'cells' (default): n_cells == 1 -> EndocrineSingle; otherwise mean caliper diameter
 *      >= ISLET_MIN_DIAM_UM -> Islet, else EndocrineCluster. SINGLE_RULE 'caliper': script 1's size rule on the
 *      new outline (<= SINGLE_MAX_DIAM_UM single, < ISLET_MIN_DIAM_UM cluster). The mean caliper diameter is
 *      computed exactly as script 1: (Max diameter + Min diameter) / 2 from QuPath's shape features ('Max
 *      diameter' = diameter of the minimum bounding circle, 'Min diameter' = minimum caliper width).
 *   4. measurements: QuPath's shape features (as script 1), Mean caliper / Equivalent diameter; 'Endocrine
 *      route' = union (bitwise OR: 1 CD99, 2 hormone, 3 both) of the constituent pixel objects' routes;
 *      'Bridge flag' = 1 if any constituent pixel object carried it; 'Tissue border flag' and 'Distance to
 *      tissue border' recomputed on the new outline with script 1's rule (analysis area = Tissue minus
 *      pLN/Artifact, 5 um inside the image edge; flag <= BORDER_UM). New: n_cells, n_pixel_objects (pixel
 *      objects that fed it), cover_min (the COVER_MIN used) and min_cell_cover (the lowest coverage among its
 *      member cells, >= cover_min by construction). Script 1's pixel-composition measurements (CD99+/hormone
 *      areas, CD99-only fraction, cores, neck split, positive channels, clamped-only) described the halo-
 *      inclusive pixel region and are NOT carried over.
 *   5. replacement: everything is built and verified first (every object has >= 1 cell; singles have exactly
 *      one under the 'cells' rule; every member cell lies inside its object; no other cell is partly inside;
 *      every member's coverage >= COVER_MIN). Only then are the pixel objects (and any IsletExpanded bands)
 *      removed and the new objects added. Tissue, pLN and Artifact are never touched: their annotation hash
 *      (the anchors' ANCHOR_CHECK hash, same canonical text as script 1) is compared before and after, and a
 *      difference stops the run. Detections left under the root are re-parented as script 1 does.
 *   6. log: counts per class before/after, n_cells per class, pixel objects dropped (no kept cell: halo only),
 *      merges (new objects fed by >= 2 pixel objects), splits (pixel objects feeding >= 2 new objects), the
 *      old class -> new class table, the largest objects before/after; all of it in V3_RESULT (provenance).
 *
 * Objects that already carry n_cells (this step ran) are refused: re-run script 1 first (run_v3 --rerun-from 1)
 * or pass ALLOW_CELL_INPUT=true. Nothing is saved here (the run_v3 footer or `QuPath script -s` saves); no
 * top-level `return`, so the footer is always reached.
 *
 * Parameters: edit the block below, pass script args KEY=VALUE (headless: -a COVER_MIN=0.6), or let run_v3
 * rewrite a declaration.
 */

import qupath.lib.analysis.features.ObjectMeasurements
import qupath.lib.io.GsonTools
import qupath.lib.io.PathIO
import qupath.lib.objects.PathObjects
import qupath.lib.objects.PathObjectTools
import qupath.lib.objects.classes.PathClass
import qupath.lib.regions.ImagePlane
import qupath.lib.roi.GeometryTools
import qupath.lib.roi.ROIs
import com.google.gson.JsonParser
import org.locationtech.jts.algorithm.MinimumBoundingCircle
import org.locationtech.jts.algorithm.MinimumDiameter
import org.locationtech.jts.geom.Coordinate
import org.locationtech.jts.geom.Envelope
import org.locationtech.jts.geom.Geometry
import org.locationtech.jts.geom.Location
import org.locationtech.jts.geom.Polygon
import org.locationtech.jts.geom.prep.PreparedGeometry
import org.locationtech.jts.geom.prep.PreparedGeometryFactory
import org.locationtech.jts.geom.util.PolygonExtracter
import org.locationtech.jts.index.strtree.STRtree
import org.locationtech.jts.operation.distance.IndexedFacetDistance
import org.locationtech.jts.operation.overlayng.OverlayNG
import org.locationtech.jts.operation.overlayng.OverlayNGRobust
import java.security.MessageDigest
import java.util.concurrent.ConcurrentLinkedQueue
import java.util.concurrent.atomic.AtomicLong
import java.util.function.IntConsumer
import java.util.stream.IntStream

evaluate(new File(buildFilePath(PROJECT_BASE_DIR, 'scripts', 'lib', 'common.groovy')))

// ==================== PARAMETERS ====================
double COVER_MIN          = 0.5      // keep a cell when >= this fraction of its area lies inside one pixel object
String SINGLE_RULE        = 'cells'  // 'cells': a 1-cell object is an EndocrineSingle | 'caliper': script 1's size rule
double TOUCH_GAP_PX       = 1.0      // kept cells whose outlines are <= this far apart (full-resolution px) touch
double ISLET_MIN_DIAM_UM  = 35.0     // Islet: mean caliper diameter >= this (script 1)
double SINGLE_MAX_DIAM_UM = 10.0     // SINGLE_RULE 'caliper' only: EndocrineSingle <= this (script 1)
double DEBRIS_MAX_UM2     = 0.0      // optional, OFF by default: cell outlines smaller than this never count
String EXCLUDE_CLASSES    = 'pLN,Artifact'         // script 1's analysis area, for the tissue-border distance
double BORDER_UM          = 10.0     // Tissue border flag: <= this from the analysis-area edge (script 1)
double EDGE_MARGIN_UM     = 5.0      // script 1: the analysis area ends this far inside the image edge
String HASH_CLASSES       = 'Tissue,pLN,Artifact'  // the anchors' annotation hash; must not change
boolean APPLY             = true     // false: build, verify and report only; the hierarchy is untouched
boolean ALLOW_CELL_INPUT  = false    // true: accept objects that already carry n_cells (1c output) as input
String WRITE_GEOJSON      = ''       // tests/review: <base>.geojson (new), <base>.pixel_objects.geojson,
                                     // <base>.members.tsv, <base>.pixel_objects.tsv, <base>.summary.json
String GAP_SENSITIVITY_PX = '0,1.5'  // report only: how many objects other touch gaps would give
int THREADS               = 0        // parallelism; 0 = the ForkJoin common pool (run_v3 sets it)
// ====================================================

def P = new LinkedHashMap([
        COVER_MIN: COVER_MIN, SINGLE_RULE: SINGLE_RULE, TOUCH_GAP_PX: TOUCH_GAP_PX, ISLET_MIN_DIAM_UM: ISLET_MIN_DIAM_UM,
        SINGLE_MAX_DIAM_UM: SINGLE_MAX_DIAM_UM, DEBRIS_MAX_UM2: DEBRIS_MAX_UM2, EXCLUDE_CLASSES: EXCLUDE_CLASSES,
        BORDER_UM: BORDER_UM, EDGE_MARGIN_UM: EDGE_MARGIN_UM, HASH_CLASSES: HASH_CLASSES, APPLY: APPLY,
        ALLOW_CELL_INPUT: ALLOW_CELL_INPUT, WRITE_GEOJSON: WRITE_GEOJSON, GAP_SENSITIVITY_PX: GAP_SENSITIVITY_PX,
        THREADS: THREADS])
// Script args KEY=VALUE override the block above; an unknown key is an error (catches typos).
def scriptArgs = (binding.hasVariable('args') && args != null) ? (args as List).collect { it.toString() } : []
scriptArgs.each { String a ->
    int eq = a.indexOf('=')
    if (eq <= 0) throw new IllegalArgumentException("script argument '${a}' is not KEY=VALUE")
    String k = a.substring(0, eq).trim(), v = a.substring(eq + 1).trim()
    if (!P.containsKey(k)) throw new IllegalArgumentException("unknown parameter '${k}'; known: ${P.keySet()}")
    def old = P[k]
    if (old instanceof Boolean && !(v.toLowerCase() in ['true', 'false']))
        throw new IllegalArgumentException("parameter ${k} takes true or false, not '${v}'")
    P[k] = (old instanceof Double) ? Double.parseDouble(v) : (old instanceof Integer) ? Integer.parseInt(v) :
            (old instanceof Boolean) ? Boolean.parseBoolean(v) : v
}
println "[1c] Parameters: ${P}"
if (!(P.SINGLE_RULE in ['cells', 'caliper']))
    throw new IllegalArgumentException("SINGLE_RULE must be 'cells' or 'caliper', not '${P.SINGLE_RULE}'")
double coverMin = P.COVER_MIN as double
if (!(coverMin > 0 && coverMin <= 1)) throw new IllegalArgumentException("COVER_MIN must be in (0, 1], not ${coverMin}")
double gap = P.TOUCH_GAP_PX as double
if (!(gap >= 0)) throw new IllegalArgumentException("TOUCH_GAP_PX must be >= 0, not ${gap}")
def sensGaps = (P.GAP_SENSITIVITY_PX as String).split(',').collect { it.trim() }.findAll { it }.collect { it as double }
if ((P.THREADS as int) > 0) qupath.lib.common.ThreadTools.setParallelism(P.THREADS as int)

long tStart = System.currentTimeMillis()
def timings = new LinkedHashMap()
def lap = { String k, long since -> timings[k] = (System.currentTimeMillis() - since) / 1000.0d }
def ENDOCRINE = ['Islet', 'EndocrineCluster', 'EndocrineSingle']
def imageData = getCurrentImageData()
def hierarchy = imageData.getHierarchy()
def server = imageData.getServer()
def cal = server.getPixelCalibration()
double sx = cal.getPixelWidthMicrons(), sy = cal.getPixelHeightMicrons()
if (Double.isNaN(sx) || Double.isNaN(sy)) throw new IllegalStateException('the image has no pixel size (script 1 sets it)')
double pxUm = Math.sqrt(sx * sy)   // um per pixel; ROI geometries are in pixel units
double pxArea = sx * sy            // um^2 per pixel^2
String UM = cal.getPixelWidthUnit()   // 'µm': measurement names match QuPath's (and script 1's)
def entry = getProjectEntry()
String imageName = entry?.getImageName() ?: server.getMetadata().getName()
String tag = imageName.replaceAll(/[^A-Za-z0-9._-]+/, '_')
def gf = GeometryTools.getDefaultFactory()
def cls = { o -> o?.getPathClass()?.toString() }
def mval = { o, String k -> o.getMeasurementList().get(k) }
def r2 = { double v -> Math.round(v * 100.0d) / 100.0d }

// JSON-safe copy (NaN/Inf -> null; GStrings -> String)
def plain
plain = { v ->
    if (v == null) return null
    if (v instanceof CharSequence) return v.toString()
    if (v instanceof Double || v instanceof Float) return (Double.isNaN(v as double) || Double.isInfinite(v as double)) ? null : v
    if (v instanceof Map) { def mm = new LinkedHashMap(); v.each { k, x -> mm[String.valueOf(k)] = plain(x) }; return mm }
    if (v instanceof Collection) return v.collect { plain(it) }
    return v
}

// ---- annotation hash of Tissue/pLN/Artifact: script 1's canonical text (islet_anchors.annotation_signature) ----
def annotationHash = { List<String> classes ->
    def anns = getAnnotationObjects().findAll { cls(it) in classes && it.getROI() != null }
    def sigLines = []
    if (!anns.isEmpty()) {
        def bos = new ByteArrayOutputStream()
        PathIO.exportObjectsAsGeoJSON(bos, anns, PathIO.GeoJsonExportOptions.FEATURE_COLLECTION)
        def fc = JsonParser.parseString(bos.toString('UTF-8')).getAsJsonObject()
        fc.getAsJsonArray('features').each { f ->
            def fo = f.getAsJsonObject()
            def cl = fo.getAsJsonObject('properties')?.getAsJsonObject('classification')
            String c = cl?.get('name')?.getAsString()
            if (!(c in classes)) return
            def g = fo.getAsJsonObject('geometry')
            String type = g?.get('type')?.getAsString()
            def polys = type == 'Polygon' ? [g.getAsJsonArray('coordinates')] :
                    type == 'MultiPolygon' ? g.getAsJsonArray('coordinates').collect { it.getAsJsonArray() } : []
            polys.each { poly ->
                def rings = poly.collect { it.getAsJsonArray() }
                if (rings.isEmpty() || rings[0].size() < 3) return
                String body = rings.collect { ring -> ring.collect { pt -> def xy = pt.getAsJsonArray()
                    String.format(Locale.ROOT, '%.2f,%.2f', xy.get(0).getAsDouble(), xy.get(1).getAsDouble()) }.join(';') }.join('/')
                sigLines << (c + '|' + body)
            }
        }
    }
    sigLines.sort()
    String text = 'classes=' + new ArrayList(classes).sort().join(',') + '\n' + sigLines.join('\n')
    MessageDigest.getInstance('SHA-256').digest(text.getBytes('UTF-8')).collect { String.format('%02x', it) }.join('').substring(0, 16)
}
def hashClasses = (P.HASH_CLASSES as String).split(',').collect { it.trim() }.findAll { it }
String hashBefore = annotationHash(hashClasses)
def hashCountsBefore = new TreeMap(getAnnotationObjects().findAll { cls(it) in hashClasses }.countBy { cls(it) })
println "[1c] ${imageName}: ${hashClasses.join('/')} annotation hash ${hashBefore} ${hashCountsBefore}"

// ---- inputs ----
long t1 = System.currentTimeMillis()
def pixObjs = getAnnotationObjects().findAll { cls(it) in ENDOCRINE }
def notArea = pixObjs.findAll { !(it.getROI()?.isArea()) }
if (!notArea.isEmpty())
    throw new IllegalStateException("${notArea.size()} islet-class annotation(s) without an area ROI: ${notArea.take(5)}")
def already = pixObjs.findAll { !Double.isNaN(mval(it, 'n_cells')) }
if (!already.isEmpty() && !(P.ALLOW_CELL_INPUT as boolean))
    throw new IllegalStateException("[1c] ${already.size()} of ${pixObjs.size()} islet-class objects already carry n_cells: " +
            "they are 1c output, not script-1 pixel objects. Re-run script 1 first (run_v3 --rerun-from 1), or pass ALLOW_CELL_INPUT=true.")
def cells = new ArrayList(hierarchy.getCellObjects())
if (cells.isEmpty()) throw new IllegalStateException('[1c] no cell detections: this step needs the cells pass (step 0) first')
// deterministic pixel-object order (script 1's own order is not guaranteed after a save/load round trip)
pixObjs.sort { a, b ->
    a.getROI().getBoundsX() <=> b.getROI().getBoundsX() ?: a.getROI().getBoundsY() <=> b.getROI().getBoundsY() ?:
            a.getROI().getArea() <=> b.getROI().getArea() ?: a.getID().toString() <=> b.getID().toString()
}
int nP = pixObjs.size(), nC = cells.size()
def plane = nP > 0 ? pixObjs[0].getROI().getImagePlane() : ImagePlane.getDefaultPlane()
println "[1c] inputs: ${nP} pixel objects ${new TreeMap(pixObjs.countBy { cls(it) })}, ${nC} cells"

Geometry[] pGeom = new Geometry[nP]
PreparedGeometry[] pPrep = new PreparedGeometry[nP]
double[] pAreaUm2 = new double[nP]
def pTree = new STRtree()
for (int i = 0; i < nP; i++) {
    Geometry g = pixObjs[i].getROI().getGeometry()
    if (!g.isValid()) g = g.buffer(0)
    pGeom[i] = g
    pPrep[i] = PreparedGeometryFactory.prepare(g)
    pAreaUm2[i] = g.getArea() * pxArea
    pTree.insert(g.getEnvelopeInternal(), i)
}
pTree.build()
lap('inputs_s', t1)

// ---- 1. coverage: (cell, pixel object) pairs with area(cell n object) > 0 ----
long t2 = System.currentTimeMillis()
double debrisPx2 = (P.DEBRIS_MAX_UM2 as double) / pxArea
def pairQ = new ConcurrentLinkedQueue()     // [cell index, pixel object index, coverage]
def nZeroArea = new AtomicLong(), nDebrisSkipped = new AtomicLong(), nOverlayRetry = new AtomicLong()
if (nP > 0) {
    IntStream.range(0, nC).parallel().forEach({ int ci ->
        def roi = cells[ci].getROI()
        if (roi == null) return
        def env = new Envelope(roi.getBoundsX(), roi.getBoundsX() + roi.getBoundsWidth(),
                roi.getBoundsY(), roi.getBoundsY() + roi.getBoundsHeight())
        def cand = pTree.query(env)
        if (cand.isEmpty()) return
        Geometry cg = roi.getGeometry()
        double ca = cg.getArea()
        if (!(ca > 0)) { nZeroArea.incrementAndGet(); return }
        if (debrisPx2 > 0 && ca < debrisPx2) { nDebrisSkipped.incrementAndGet(); return }
        for (def pi0 : cand) {
            int pi = pi0 as int
            if (!pPrep[pi].intersects(cg)) continue
            double cov
            if (pPrep[pi].covers(cg)) cov = 1.0d
            else {
                Geometry inter
                try { inter = OverlayNGRobust.overlay(cg, pGeom[pi], OverlayNG.INTERSECTION) }
                catch (Exception e) { nOverlayRetry.incrementAndGet(); inter = OverlayNGRobust.overlay(cg.buffer(0), pGeom[pi], OverlayNG.INTERSECTION) }
                cov = Math.min(1.0d, inter.getArea() / ca)
            }
            if (cov > 0) pairQ.add([ci, pi, cov] as Object[])
        }
    } as IntConsumer)
}
def pairs = new ArrayList(pairQ)
pairs.sort { a, b -> (a[0] as int) <=> (b[0] as int) ?: (a[1] as int) <=> (b[1] as int) }
// per cell: pixel objects that keep it (coverage >= COVER_MIN), best coverage, total coverage
def keptBy = new TreeMap()      // cell index -> [pixel object indices]
def bestCov = [:]               // cell index -> max coverage over pixel objects
def totalCov = [:].withDefault { 0.0d }
pairs.each { r ->
    int ci = r[0] as int, pi = r[1] as int
    double cov = r[2] as double
    totalCov[ci] = (totalCov[ci] as double) + cov
    if (cov > ((bestCov[ci] ?: 0.0d) as double)) bestCov[ci] = cov
    if (cov >= coverMin) { if (!keptBy.containsKey(ci)) keptBy[ci] = []; keptBy[ci] << pi }
}
double[] pBest = new double[nP]   // per pixel object: the best coverage of any cell (report: why objects drop)
pairs.each { r -> int pi = r[1] as int; if ((r[2] as double) > pBest[pi]) pBest[pi] = r[2] as double }
int nCellsTouching = totalCov.size()
int nStraddle = totalCov.count { ci, tc -> !keptBy.containsKey(ci) && (tc as double) >= coverMin }   // kept only if coverage were summed over objects
int nKeptByTwo = keptBy.count { ci, lst -> lst.size() >= 2 }
int[] keptCells = keptBy.keySet() as int[]
int nK = keptCells.length
def pKept = new int[nP]
keptBy.each { ci, lst -> lst.each { pKept[it as int]++ } }
lap('coverage_s', t2)
println String.format(Locale.ROOT, "[1c] coverage: %d cells intersect a pixel object, %d kept (coverage >= %.2f), %d kept by two objects, " +
        "%d straddle objects (kept only if summed), %d zero-area outlines skipped%s (%.1f s)", nCellsTouching, nK, coverMin, nKeptByTwo,
        nStraddle, nZeroArea.get(), (debrisPx2 > 0 ? ", ${nDebrisSkipped.get()} debris < ${P.DEBRIS_MAX_UM2} um^2 skipped" : ''), timings.coverage_s)

// ---- 2. components: kept cells within TOUCH_GAP_PX of each other ----
long t3 = System.currentTimeMillis()
Geometry[] kGeom = new Geometry[nK]
for (int k = 0; k < nK; k++) {
    Geometry g = cells[keptCells[k]].getROI().getGeometry()
    kGeom[k] = g.isValid() ? g : g.buffer(0)
}
double maxGap = ([gap] + sensGaps).max() as double
def kTree = new STRtree()
for (int k = 0; k < nK; k++) kTree.insert(kGeom[k].getEnvelopeInternal(), k)
kTree.build()
// links [k, j, distance, needsFill]; needsFill = touching within TOUCH_GAP_PX but not sharing an edge or overlapping
def linkQ = new ConcurrentLinkedQueue()
IntStream.range(0, nK).parallel().forEach({ int k ->
    def env = new Envelope(kGeom[k].getEnvelopeInternal())
    env.expandBy(maxGap)
    for (def j0 : kTree.query(env)) {
        int j = j0 as int
        if (j <= k) continue
        if (kGeom[k].getEnvelopeInternal().distance(kGeom[j].getEnvelopeInternal()) > maxGap) continue
        double d = kGeom[k].distance(kGeom[j])
        if (d > maxGap) continue
        boolean fill = false
        if (d <= gap) {
            if (d > 0) fill = true
            else {
                def im = kGeom[k].relate(kGeom[j])
                fill = !(im.get(Location.BOUNDARY, Location.BOUNDARY) >= 1 || im.get(Location.INTERIOR, Location.INTERIOR) >= 0)
            }
        }
        linkQ.add([k, j, d, fill] as Object[])
    }
} as IntConsumer)
def links = new ArrayList(linkQ)
links.sort { a, b -> (a[0] as int) <=> (b[0] as int) ?: (a[1] as int) <=> (b[1] as int) }
def componentsAt = { double g0 ->
    int[] parent = new int[nK]
    for (int i = 0; i < nK; i++) parent[i] = i
    def find = { int i -> while (parent[i] != i) { parent[i] = parent[parent[i]]; i = parent[i] }; i }
    links.each { l -> if ((l[2] as double) <= g0) { int a = find(l[0] as int), b = find(l[1] as int); if (a != b) parent[Math.max(a, b)] = Math.min(a, b) } }
    // cells kept by two pixel objects are one node already; groups ordered by their first member
    def groups = new TreeMap()
    for (int i = 0; i < nK; i++) { int r = find(i); if (!groups.containsKey(r)) groups[r] = []; groups[r] << i }
    new ArrayList(groups.values())
}
def comps = componentsAt(gap)
def sensitivity = new LinkedHashMap()
sensGaps.each { double g0 -> sensitivity[String.valueOf(g0)] = componentsAt(g0).size() }
int nFillLinks = links.count { it[3] as boolean }
lap('components_s', t3)
println "[1c] components: ${comps.size()} at TOUCH_GAP_PX=${gap} (${links.count { (it[2] as double) <= gap }} touching pairs, " +
        "${nFillLinks} joined across a gap or a corner); at other gaps ${sensitivity} (${timings.components_s} s)"

// ---- geometry of one component: member outlines + the gaps between touching members, holes filled ----
def linksByComp = [:]
def compOf = new int[nK]
comps.eachWithIndex { List mem, int c -> mem.each { compOf[it as int] = c } }
links.each { l -> if (l[3] as boolean) { int c = compOf[l[0] as int]; if (!linksByComp.containsKey(c)) linksByComp[c] = []; linksByComp[c] << l } }
final double FILL_EPS_PX = 0.05d
// every cell within reach of a gap fill (kept or not): a fill never covers a cell that is not a member of its component
long t4 = System.currentTimeMillis()
def nearQ = new ConcurrentLinkedQueue()
if (nK > 0) {
    IntStream.range(0, nC).parallel().forEach({ int ci ->
        def roi = cells[ci].getROI()
        if (roi == null) return
        def env = new Envelope(roi.getBoundsX(), roi.getBoundsX() + roi.getBoundsWidth(), roi.getBoundsY(), roi.getBoundsY() + roi.getBoundsHeight())
        env.expandBy(gap + 2 * FILL_EPS_PX)
        if (!kTree.query(env).isEmpty()) nearQ.add(ci)
    } as IntConsumer)
}
int[] nearCells = (new ArrayList(nearQ)).sort() as int[]
Geometry[] nearGeom = new Geometry[nearCells.length]
def nTree = new STRtree()
for (int q = 0; q < nearCells.length; q++) {
    Geometry g = cells[nearCells[q]].getROI().getGeometry()
    nearGeom[q] = g.isValid() ? g : g.buffer(0)
    nTree.insert(nearGeom[q].getEnvelopeInternal(), q)
}
nTree.build()
def fillHoles = { Geometry g ->
    def polys = new ArrayList(PolygonExtracter.getPolygons(g)).findAll { !it.isEmpty() && it.getArea() > 0 }
    def shells = polys.collect { Polygon p -> gf.createPolygon(p.getExteriorRing().getCoordinates()) }
    // a part inside another part's filled hole is dropped (it is covered)
    shells = shells.findAll { s -> !shells.any { o -> !o.is(s) && o.getArea() > s.getArea() && o.covers(s) } }
    shells
}
def nMultipart = new AtomicLong(), nPartFills = new AtomicLong(), nTrimmedFills = new AtomicLong(), nStrayParts = new AtomicLong()
def buildGeom = { int c ->
    List mem = comps[c]
    def memCells = new HashSet(mem.collect { keptCells[it as int] })
    // the region within (their distance + FILL_EPS_PX) of both a and b, minus every non-member cell; only pieces that
    // still touch both a and b are kept (the bridge), so a gap fill never takes part of another cell
    def bridge = { Geometry a, Geometry b, double d ->
        double r = d + FILL_EPS_PX
        Geometry f = OverlayNGRobust.overlay(a.buffer(r), b.buffer(r), OverlayNG.INTERSECTION)
        if (f.isEmpty()) return []
        def blockers = []
        for (def q0 : nTree.query(f.getEnvelopeInternal())) {
            int q = q0 as int
            if (!memCells.contains(nearCells[q]) && nearGeom[q].intersects(f)) blockers << nearGeom[q]
        }
        if (!blockers.isEmpty()) {
            nTrimmedFills.incrementAndGet()
            f = OverlayNGRobust.overlay(f, blockers.size() == 1 ? blockers[0] : OverlayNGRobust.union(blockers), OverlayNG.DIFFERENCE)
        }
        new ArrayList(PolygonExtracter.getPolygons(f)).findAll { Polygon fp -> fp.getArea() > 0 && fp.intersects(a) && fp.intersects(b) }
    }
    def parts = mem.collect { kGeom[it as int] }
    (linksByComp[c] ?: []).each { l -> parts.addAll(bridge(kGeom[l[0] as int], kGeom[l[1] as int], l[2] as double)) }
    def shells = fillHoles(parts.size() == 1 ? parts[0] : OverlayNGRobust.union(parts))
    // parts still apart but within TOUCH_GAP_PX of each other touch too (a cell's own detached fragment; a bridge cut
    // by a non-member cell): the same trimmed fill between the two parts
    int guard = 0
    while (shells.size() > 1 && guard++ < 5) {
        def extra = []
        for (int i = 0; i < shells.size(); i++) for (int j = i + 1; j < shells.size(); j++) {
            double d = shells[i].distance(shells[j])
            if (d <= gap) extra.addAll(bridge(shells[i], shells[j], d))
        }
        if (extra.isEmpty()) break
        nPartFills.addAndGet(extra.size())
        shells = fillHoles(OverlayNGRobust.union(shells + extra))
    }
    if (shells.size() > 1) {   // safety net: a part without member material (not expected) is dropped
        def memU = OverlayNGRobust.union(mem.collect { kGeom[it as int] })
        def keep = shells.findAll { s -> OverlayNGRobust.overlay(s, memU, OverlayNG.INTERSECTION).getArea() > 0 }
        nStrayParts.addAndGet(shells.size() - keep.size())
        shells = keep
    }
    // parts farther apart than TOUCH_GAP_PX (a cell's far fragment) stay one multi-part object: whole cells, never cut
    if (shells.size() > 1) { nMultipart.incrementAndGet(); return gf.createMultiPolygon(shells as Polygon[]) }
    shells[0]
}
def compGeom = new Geometry[comps.size()]
IntStream.range(0, comps.size()).parallel().forEach({ int c -> compGeom[c] = buildGeom(c) } as IntConsumer)
// a component inside another's filled hole (or, not expected, overlapping it) joins that component
def cTree = new STRtree()
compGeom.eachWithIndex { g, int c -> cTree.insert(g.getEnvelopeInternal(), c) }
cTree.build()
int[] into = new int[comps.size()]
for (int c = 0; c < comps.size(); c++) into[c] = c
int nSwallowed = 0, nOverlapMerged = 0
def byArea = (0..<comps.size()).toList().sort { a, b -> compGeom[b].getArea() <=> compGeom[a].getArea() ?: a <=> b }
def prepC = [:]
byArea.each { int c ->   // largest first: a covering component is always processed before the ones it covers
    if (into[c] != c) return
    for (def o0 : cTree.query(compGeom[c].getEnvelopeInternal())) {
        int o = o0 as int
        if (o == c || into[o] != o || compGeom[o].getArea() > compGeom[c].getArea()) continue
        if (compGeom[o].getArea() == compGeom[c].getArea() && o < c) continue
        def pc = prepC[c] ?: (prepC[c] = PreparedGeometryFactory.prepare(compGeom[c]))
        if (!pc.intersects(compGeom[o])) continue
        if (pc.covers(compGeom[o])) { into[o] = c; nSwallowed++ }
        else if (OverlayNGRobust.overlay(compGeom[c], compGeom[o], OverlayNG.INTERSECTION).getArea() > 1e-9) { into[o] = c; nOverlapMerged++ }
    }
}
def absorbedBy = [:]
for (int c = 0; c < comps.size(); c++) if (into[c] != c) { if (!absorbedBy.containsKey(into[c])) absorbedBy[into[c]] = []; absorbedBy[into[c]] << c }
def finals = []   // [geom, member k list]
for (int c = 0; c < comps.size(); c++) {
    if (into[c] != c) continue
    def mem = new ArrayList(comps[c])
    def absorbed = absorbedBy[c] ?: []
    absorbed.each { mem.addAll(comps[it as int]) }
    Geometry g = compGeom[c]
    if (absorbed.any { !PreparedGeometryFactory.prepare(compGeom[c]).covers(compGeom[it]) }) {
        def sh = fillHoles(OverlayNGRobust.union([compGeom[c]] + absorbed.collect { compGeom[it] }))
        g = sh.size() == 1 ? sh[0] : gf.createMultiPolygon(sh as Polygon[])
    }
    finals << [g, mem.sort()]
}
finals.sort { a, b ->   // deterministic object order
    def ea = a[0].getEnvelopeInternal(), eb = b[0].getEnvelopeInternal()
    ea.getMinX() <=> eb.getMinX() ?: ea.getMinY() <=> eb.getMinY() ?: a[0].getArea() <=> b[0].getArea()
}
lap('geometry_s', t4)
println "[1c] geometry: ${finals.size()} objects (${nSwallowed} components inside another's filled hole joined it, " +
        "${nOverlapMerged} overlapping merged; ${nTrimmedFills.get()} gap fills trimmed by a non-member cell, ${nPartFills.get()} fills " +
        "between parts, ${nStrayParts.get()} stray parts dropped, ${nMultipart.get()} multi-part) (${timings.geometry_s} s)"

// ---- 3./4. objects, measurements, classes ----
long t5 = System.currentTimeMillis()
def pcIslet = PathClass.fromString(CLASS_ISLET), pcCluster = PathClass.fromString('EndocrineCluster'), pcSingle = PathClass.fromString('EndocrineSingle')
def newObjs = finals.collect { PathObjects.createAnnotationObject(GeometryTools.geometryToROI(it[0] as Geometry, plane)) }
ObjectMeasurements.addShapeMeasurements(newObjs, cal,
        ObjectMeasurements.ShapeFeatures.AREA, ObjectMeasurements.ShapeFeatures.LENGTH,
        ObjectMeasurements.ShapeFeatures.CIRCULARITY, ObjectMeasurements.ShapeFeatures.SOLIDITY,
        ObjectMeasurements.ShapeFeatures.MAX_DIAMETER, ObjectMeasurements.ShapeFeatures.MIN_DIAMETER)
// analysis area of script 1 (Tissue minus EXCLUDE_CLASSES, EDGE_MARGIN_UM inside the image) for the border distance
def tissueAnns = getAnnotationObjects().findAll { cls(it) == CLASS_TISSUE && it.getROI()?.isArea() }
if (tissueAnns.isEmpty()) throw new IllegalStateException("[1c] no ${CLASS_TISSUE} annotation")
Geometry area = OverlayNGRobust.union(tissueAnns.collect { it.getROI().getGeometry() })
def excludeNames = (P.EXCLUDE_CLASSES as String).split(',').collect { it.trim() }.findAll { it }
def excluded = getAnnotationObjects().findAll { cls(it) in excludeNames && it.getROI()?.isArea() }
if (!excluded.isEmpty()) {
    def ex = OverlayNGRobust.union(excluded.collect { it.getROI().getGeometry() })
    try { area = area.difference(ex) } catch (Exception e) { area = area.buffer(0).difference(ex.buffer(0)) }
}
double mX = (P.EDGE_MARGIN_UM as double) / sx, mY = (P.EDGE_MARGIN_UM as double) / sy
def inner = ROIs.createRectangleROI(mX, mY, Math.max(0.0d, server.getWidth() - 2 * mX), Math.max(0.0d, server.getHeight() - 2 * mY), plane).getGeometry()
def clipGeom = toPolygonal(area.intersection(inner))
if (clipGeom == null || clipGeom.isEmpty()) throw new IllegalStateException('[1c] Tissue lies outside the analysis area')
def borderDist = new IndexedFacetDistance(clipGeom.getBoundary())
borderDist.distance(gf.createPoint(clipGeom.getInteriorPoint().getCoordinate()))   // build the index before parallel use

String kMax = "Max diameter ${UM}", kMin = "Min diameter ${UM}"
double[] cellCov = new double[nK]
for (int k = 0; k < nK; k++) cellCov[k] = bestCov[keptCells[k]] as double
def objPix = new ArrayList()   // per new object: sorted pixel object indices that fed it
finals.eachWithIndex { f, int i ->
    def s = new TreeSet()
    f[1].each { k -> keptBy[keptCells[k as int]].each { s << (it as int) } }
    objPix << new ArrayList(s)
}
IntStream.range(0, newObjs.size()).parallel().forEach({ int i ->
    def o = newObjs[i]
    Geometry g = finals[i][0] as Geometry
    def mem = finals[i][1] as List
    def ml = o.getMeasurementList()
    double dMax = ml.get(kMax), dMin = ml.get(kMin)
    if (Double.isNaN(dMax) || Double.isNaN(dMin)) {   // QuPath's own definitions (script 1's fallback)
        dMax = 2.0d * new MinimumBoundingCircle(g).getRadius() * pxUm
        dMin = new MinimumDiameter(g).getLength() * pxUm
    }
    double meanCal = (dMax + dMin) / 2.0d
    double areaUm2 = g.getArea() * pxArea
    int route = 0, nRoute = 0
    int bridge = 0
    objPix[i].each { int pi ->
        double r = mval(pixObjs[pi], 'Endocrine route')
        if (!Double.isNaN(r)) { route |= (int) r; nRoute++ }
        if (mval(pixObjs[pi], 'Bridge flag') == 1.0d) bridge = 1
    }
    double minCov = mem.collect { cellCov[it as int] }.min() as double
    double dBorder = borderDist.distance(g) * pxUm
    ml.put("Mean caliper diameter ${UM}", meanCal)
    ml.put("Equivalent diameter ${UM}", 2.0d * Math.sqrt(areaUm2 / Math.PI))
    ml.put('Endocrine route', nRoute > 0 ? (double) route : Double.NaN)
    ml.put('Bridge flag', bridge)
    ml.put('Tissue border flag', dBorder <= (P.BORDER_UM as double) ? 1 : 0)
    ml.put("Distance to tissue border ${UM}", dBorder)
    ml.put('n_cells', mem.size())
    ml.put('n_pixel_objects', objPix[i].size())
    ml.put('cover_min', coverMin)
    ml.put('min_cell_cover', minCov)
    ml.close()
    if (P.SINGLE_RULE == 'cells')
        o.setPathClass(mem.size() == 1 ? pcSingle : (meanCal >= (P.ISLET_MIN_DIAM_UM as double) ? pcIslet : pcCluster))
    else
        o.setPathClass(meanCal >= (P.ISLET_MIN_DIAM_UM as double) ? pcIslet : (meanCal > (P.SINGLE_MAX_DIAM_UM as double) ? pcCluster : pcSingle))
} as IntConsumer)
lap('measurements_s', t5)

// ---- verification, before anything in the hierarchy changes ----
long t6 = System.currentTimeMillis()
def problems = []
int nObjNoCell = newObjs.count { mval(it, 'n_cells') < 1 }
if (nObjNoCell > 0) problems << "${nObjNoCell} object(s) without cells"
if (P.SINGLE_RULE == 'cells') {
    int bad = newObjs.count { (cls(it) == 'EndocrineSingle') != (mval(it, 'n_cells') == 1.0d) }
    if (bad > 0) problems << "${bad} object(s) break the single rule (EndocrineSingle <=> n_cells == 1)"
}
int nLowCov = newObjs.count { mval(it, 'min_cell_cover') < coverMin }
if (nLowCov > 0) problems << "${nLowCov} object(s) with a member cell below COVER_MIN"
int nMembers = finals.sum(0) { (it[1] as List).size() } as int
if (nMembers != nK) problems << "members ${nMembers} != kept cells ${nK}"
// every cell against the new objects: members must lie inside; other cells either inside a filled hole or outside
def oTree = new STRtree()
Geometry[] oGeom = finals.collect { it[0] as Geometry } as Geometry[]
PreparedGeometry[] oPrep = oGeom.collect { PreparedGeometryFactory.prepare(it) } as PreparedGeometry[]
oGeom.eachWithIndex { g, int i -> oTree.insert(g.getEnvelopeInternal(), i) }
oTree.build()
def memberOf = new HashMap()   // cell index -> new object index
finals.eachWithIndex { f, int i -> (f[1] as List).each { k -> memberOf[keptCells[k as int]] = i } }
def nMemberOutside = new AtomicLong(), nHoleCells = new AtomicLong(), nPartial = new AtomicLong(), nSegOverlap = new AtomicLong()
def partialQ = new ConcurrentLinkedQueue()   // overlap fractions of non-member cells partly inside an object
// k indices of each object's members, to tell a segmentation overlap (a non-member overlapping a member) from a bug
def objOfK = new int[nK]
finals.eachWithIndex { f, int i -> (f[1] as List).each { objOfK[it as int] = i } }
if (!finals.isEmpty()) {
    IntStream.range(0, nC).parallel().forEach({ int ci ->
        def roi = cells[ci].getROI()
        if (roi == null) return
        def env = new Envelope(roi.getBoundsX(), roi.getBoundsX() + roi.getBoundsWidth(), roi.getBoundsY(), roi.getBoundsY() + roi.getBoundsHeight())
        def cand = oTree.query(env)
        if (cand.isEmpty()) return
        Geometry cg = roi.getGeometry()
        double ca = cg.getArea()
        Integer own = memberOf.get(ci)
        for (def o0 : cand) {
            int oi = o0 as int
            if (!oPrep[oi].intersects(cg)) continue
            double f
            if (oPrep[oi].covers(cg)) f = 1.0d
            else if (!(ca > 0)) f = 0.0d
            else f = OverlayNGRobust.overlay(cg.isValid() ? cg : cg.buffer(0), oGeom[oi], OverlayNG.INTERSECTION).getArea() / ca
            if (own != null && own == oi) { if (f < 1.0d - 1e-6) nMemberOutside.incrementAndGet() }
            else if (f >= 0.98d) nHoleCells.incrementAndGet()
            else if (f > 0.02d) {
                // allowed only where the two outlines themselves overlap (InstanSeg seam overlap): the member brings it in
                boolean seg = false
                for (def k0 : kTree.query(cg.getEnvelopeInternal())) {
                    int k = k0 as int
                    if (objOfK[k] == oi && keptCells[k] != ci && kGeom[k].intersects(cg) &&
                            OverlayNGRobust.overlay(kGeom[k], cg.isValid() ? cg : cg.buffer(0), OverlayNG.INTERSECTION).getArea() > 1e-6) { seg = true; break }
                }
                if (seg) nSegOverlap.incrementAndGet() else { nPartial.incrementAndGet(); partialQ.add(f) }
            }
        }
    } as IntConsumer)
}
if (nMemberOutside.get() > 0) problems << "${nMemberOutside.get()} member cell(s) not inside their object"
if (nPartial.get() > 0) problems << "${nPartial.get()} non-member cell(s) partly inside an object"
double worstPartial = partialQ.isEmpty() ? 0.0d : (partialQ.max() as double)
lap('verify_s', t6)
println "[1c] verification: ${problems.isEmpty() ? 'OK' : problems}; ${nHoleCells.get()} non-member cell(s) inside filled holes; " +
        "${nPartial.get()} non-member cell(s) partly inside an object (2-98%, worst ${String.format('%.3f', worstPartial)}); " +
        "${nSegOverlap.get()} partly inside because their outline overlaps a member's (segmentation overlap) (${timings.verify_s} s)"
if (!problems.isEmpty()) throw new IllegalStateException("[1c] verification FAILED for ${imageName}: ${problems}; nothing was changed")

// ---- 6. summary: before/after, drops, merges, splits, transitions ----
def classCountsBefore = new LinkedHashMap(ENDOCRINE.collectEntries { c -> [c, pixObjs.count { cls(it) == c }] })
def classCountsAfter = new LinkedHashMap(ENDOCRINE.collectEntries { c -> [c, newObjs.count { cls(it) == c }] })
def pNew = (0..<nP).collect { [] }      // pixel object -> new object indices it fed
objPix.eachWithIndex { List lst, int i -> lst.each { pNew[it as int] << i } }
def outcomes = ['dropped', 'Islet', 'EndocrineCluster', 'EndocrineSingle', 'split']
def transitions = new LinkedHashMap(ENDOCRINE.collectEntries { c -> [c, new LinkedHashMap(outcomes.collectEntries { [it, 0] })] })
def intoMerge = new LinkedHashMap(ENDOCRINE.collectEntries { c -> [c, 0] })   // 1:1 outcome, but the new object is a merge
def pOutcome = new String[nP]
for (int pi = 0; pi < nP; pi++) {
    String oc = pNew[pi].isEmpty() ? 'dropped' : (pNew[pi].size() >= 2 ? 'split' : cls(newObjs[pNew[pi][0] as int]))
    pOutcome[pi] = oc
    transitions[cls(pixObjs[pi])][oc] += 1
    if (pNew[pi].size() == 1 && objPix[pNew[pi][0] as int].size() >= 2) intoMerge[cls(pixObjs[pi])] += 1
}
def droppedByClass = new LinkedHashMap(ENDOCRINE.collectEntries { c -> [c, transitions[c].dropped] })
// dropped objects by the best coverage any cell reached (0.1 bins below COVER_MIN): how close they came
def droppedBest = new LinkedHashMap([('no cell'): 0])
int nBins = Math.max(1, (int) Math.ceil(coverMin / 0.1d - 1e-9))
for (int b = 0; b < nBins; b++) droppedBest[String.format(Locale.ROOT, '%.1f-%.1f', b * 0.1d, Math.min((b + 1) * 0.1d, coverMin))] = 0
def binKeys = new ArrayList(droppedBest.keySet())
for (int pi = 0; pi < nP; pi++) {
    if (!pNew[pi].isEmpty()) continue
    if (pBest[pi] == 0) { droppedBest['no cell'] += 1; continue }
    droppedBest[binKeys[1 + Math.min(nBins - 1, (int) Math.floor(pBest[pi] / 0.1d))]] += 1
}
def splitByClass = new LinkedHashMap(ENDOCRINE.collectEntries { c -> [c, transitions[c].split] })
int nMergedObjects = objPix.count { it.size() >= 2 }
int nPixInMerges = objPix.findAll { it.size() >= 2 }.sum(0) { it.size() } as int
int nSplitPix = pNew.count { it.size() >= 2 }
int nSplitPieces = pNew.findAll { it.size() >= 2 }.sum(0) { it.size() } as int
def bins = [[1, 1], [2, 2], [3, 5], [6, 10], [11, 25], [26, 50], [51, 100], [101, 250], [251, 500], [501, Integer.MAX_VALUE]]
def nCellsByClass = new LinkedHashMap()
ENDOCRINE.each { c ->
    def v = newObjs.findAll { cls(it) == c }.collect { mval(it, 'n_cells') as int }
    nCellsByClass[c] = [n: v.size(), mean: v ? (v.sum() / (double) v.size()) : null, min: v ? v.min() : null, max: v ? v.max() : null,
                        hist: new LinkedHashMap(bins.collectEntries { b -> [(b[0] == b[1] ? "${b[0]}" : (b[1] == Integer.MAX_VALUE ? ">${b[0] - 1}" : "${b[0]}-${b[1]}")),
                                                                     v.count { it >= b[0] && it <= b[1] }] })]
}
def pRows = (0..<nP).collect { pi -> [area_um2: r2(pAreaUm2[pi]), cls: cls(pixObjs[pi]), kept_cells: pKept[pi], new_objects: pNew[pi].size(), outcome: pOutcome[pi]] }
def oRows = []
for (int i = 0; i < newObjs.size(); i++) {
    def o = newObjs[i]
    oRows << [area_um2: r2(mval(o, "Area ${UM}^2".toString())), cls: cls(o), n_cells: mval(o, 'n_cells') as int,
              n_pixel_objects: objPix[i].size(), mean_caliper_um: r2(mval(o, "Mean caliper diameter ${UM}".toString()))]
}
def largestBefore = pRows.sort(false) { -(it.area_um2 as double) }.take(5)
def largestAfter = oRows.sort(false) { -(it.area_um2 as double) }.take(5)
double bigUm2 = 50000.0d
def meanArea = { List rows, String c -> def a = rows.findAll { it.cls == c }.collect { it.area_um2 as double }; a ? a.sum() / a.size() : null }
def summary = new LinkedHashMap([
        image: imageName, tag: tag, params: P, apply: P.APPLY,
        annotation_hash_before: hashBefore, annotation_hash_classes: hashClasses, annotation_counts: hashCountsBefore,
        n_cells_image: nC, n_cells_intersecting: nCellsTouching, n_cells_kept: nK, n_cells_kept_by_two_objects: nKeptByTwo,
        n_cells_straddling_unkept: nStraddle, n_zero_area_outlines_skipped: nZeroArea.get(),
        n_debris_skipped: nDebrisSkipped.get(), n_hole_cells: nHoleCells.get(), n_partial_nonmember_cells: nPartial.get(),
        n_nonmember_overlapping_member: nSegOverlap.get(),
        worst_partial_nonmember_fraction: worstPartial,
        n_pixel_objects: nP, n_objects: newObjs.size(),
        counts_before: classCountsBefore, counts_after: classCountsAfter,
        mean_area_um2_before: new LinkedHashMap(ENDOCRINE.collectEntries { c -> [c, meanArea(pRows, c)] }),
        mean_area_um2_after: new LinkedHashMap(ENDOCRINE.collectEntries { c -> [c, meanArea(oRows, c)] }),
        n_cells_by_class: nCellsByClass,
        dropped_halo_only: droppedByClass, n_dropped_halo_only: droppedByClass.values().sum(0),
        dropped_by_best_cell_cover: droppedBest,
        merges: [n_new_objects: nMergedObjects, n_pixel_objects_in_them: nPixInMerges],
        splits: [n_pixel_objects: nSplitPix, n_new_objects_from_them: nSplitPieces, by_class: splitByClass],
        transitions: transitions, one_to_one_into_a_merge: intoMerge,
        clusters_to_singles: transitions.EndocrineCluster.EndocrineSingle,
        n_objects_over_50000_um2: [before: pRows.count { (it.area_um2 as double) > bigUm2 }, after: oRows.count { (it.area_um2 as double) > bigUm2 }],
        largest_before: largestBefore, largest_after: largestAfter,
        components_at_other_gaps: sensitivity, n_fill_links: nFillLinks, n_hole_swallowed: nSwallowed,
        n_near_cells: nearCells.length, n_fills_trimmed: nTrimmedFills.get(), n_part_fills: nPartFills.get(),
        n_stray_parts: nStrayParts.get(),
        n_overlap_merged: nOverlapMerged, n_multipart: nMultipart.get(),
        n_overlay_retries: nOverlayRetry.get(),
        tissue_border_flag: newObjs.count { mval(it, 'Tissue border flag') == 1.0d },
        bridge_flag: newObjs.count { mval(it, 'Bridge flag') == 1.0d },
        route_counts: new TreeMap(newObjs.countBy { def r = mval(it, 'Endocrine route'); Double.isNaN(r) ? 'NaN' : String.valueOf((int) r) })])

// ---- optional files (tests, review) ----
def outputs = []
def baseOut = (P.WRITE_GEOJSON as String) ? (P.WRITE_GEOJSON as String).replaceAll(/(?i)\.geojson$/, '') : null
if (baseOut) {
    new File(baseOut).getAbsoluteFile().getParentFile().mkdirs()
    if (!newObjs.isEmpty()) { exportObjectsToGeoJson(newObjs, "${baseOut}.geojson", 'FEATURE_COLLECTION'); outputs << "${baseOut}.geojson".toString() }
    if (!pixObjs.isEmpty()) { exportObjectsToGeoJson(pixObjs, "${baseOut}.pixel_objects.geojson", 'FEATURE_COLLECTION'); outputs << "${baseOut}.pixel_objects.geojson".toString() }
    new File("${baseOut}.members.tsv").withWriter('UTF-8') { w ->
        w.write('cell_id\tobject_id\tobject_class\tcoverage\tkept_by\n')
        finals.eachWithIndex { f, int i ->
            (f[1] as List).each { k ->
                int ci = keptCells[k as int]
                w.write("${cells[ci].getID()}\t${newObjs[i].getID()}\t${cls(newObjs[i])}\t${cellCov[k as int]}\t" +
                        keptBy[ci].collect { pixObjs[it as int].getID().toString() }.join(';') + '\n')
            }
        }
    }
    new File("${baseOut}.pixel_objects.tsv").withWriter('UTF-8') { w ->
        w.write('pixel_object_id\tclass\tarea_um2\tkept_cells\toutcome\tnew_object_ids\n')
        for (int pi = 0; pi < nP; pi++)
            w.write("${pixObjs[pi].getID()}\t${cls(pixObjs[pi])}\t${pAreaUm2[pi]}\t${pKept[pi]}\t${pOutcome[pi]}\t" +
                    pNew[pi].collect { newObjs[it as int].getID().toString() }.join(';') + '\n')
    }
    outputs << "${baseOut}.members.tsv".toString() << "${baseOut}.pixel_objects.tsv".toString()
}

// ---- 5. replacement: only now does the hierarchy change ----
def reparented = new LinkedHashMap([(CLASS_TISSUE): 0, pLN: 0, left_at_root: 0])
if ((P.APPLY as boolean) && nP > 0) {
    long t7 = System.currentTimeMillis()
    def oldBands = getAnnotationObjects().findAll { cls(it) == CLASS_ISLET_EXPANDED }
    summary.n_removed = [pixel_objects: pixObjs.size(), bands: oldBands.size()]
    removeObjects(pixObjs + oldBands, true)
    addObjects(newObjs)
    // detections of a removed object move to ITS parent (the root, for script 1's objects): give them to the region
    // containing their nucleus centroid (Tissue first, then pLN), as script 1 does
    def root = hierarchy.getRootObject()
    def orphans = root.getChildObjects().findAll { it.isDetection() }
    if (!orphans.isEmpty()) {
        def claimed = Collections.newSetFromMap(new IdentityHashMap())
        [CLASS_TISSUE, 'pLN'].each { rc ->
            getAnnotationObjects().findAll { cls(it) == rc && it.getROI()?.isArea() }.each { region ->
                def prep = PreparedGeometryFactory.prepare(region.getROI().getGeometry())
                def kids = orphans.findAll { d -> if (claimed.contains(d)) return false
                    def r = PathObjectTools.getROI(d, true)
                    prep.contains(gf.createPoint(new Coordinate(r.getCentroidX(), r.getCentroidY()))) }
                if (!kids.isEmpty()) { region.addChildObjects(kids); claimed.addAll(kids); reparented[rc] = (reparented[rc] as int) + kids.size() }
            }
        }
        reparented.left_at_root = orphans.size() - claimed.size()
        println "[1c] Re-parented detections that were left under the root: ${reparented}"
    }
    fireHierarchyUpdate()
    lap('replace_s', t7)
}
summary.reparented_root_detections = reparented
String hashAfter = annotationHash(hashClasses)
summary.annotation_hash_after = hashAfter
if (hashAfter != hashBefore)
    throw new IllegalStateException("[1c] the ${hashClasses.join('/')} annotation hash changed (${hashBefore} -> ${hashAfter}): nothing may be saved")
def nowCounts = new LinkedHashMap(ENDOCRINE.collectEntries { c -> [c, getAnnotationObjects().count { cls(it) == c }] })
if ((P.APPLY as boolean) && nP > 0 && nowCounts != classCountsAfter)
    throw new IllegalStateException("[1c] hierarchy holds ${nowCounts} after the replacement, expected ${classCountsAfter}")
if (hierarchy.getCellObjects().size() != nC) throw new IllegalStateException('[1c] the number of cells changed')

timings.total_s = (System.currentTimeMillis() - tStart) / 1000.0d
summary.seconds = timings
try { summary.vm_hwm_kb = new File('/proc/self/status').readLines().find { it.startsWith('VmHWM:') }?.replaceAll(/[^0-9]/, '') as Long } catch (Throwable e) { }
if (baseOut) {
    def f = new File("${baseOut}.summary.json")
    f.text = GsonTools.getInstance(true).toJson(plain(summary))
    outputs << f.getPath()
}
summary.outputs = outputs

println "[1c] ${imageName}: pixel objects ${classCountsBefore} -> cell-based ${classCountsAfter}"
println "[1c]   dropped (halo only, no kept cell): ${droppedByClass}, by their best cell's coverage ${droppedBest}; merges: ${nMergedObjects} " +
        "new objects from ${nPixInMerges} pixel objects; splits: ${nSplitPix} pixel objects -> ${nSplitPieces} new objects"
println "[1c]   old class -> outcome: ${transitions}"
println "[1c]   n_cells per class: " + nCellsByClass.collect { c, v -> "${c} n=${v.n} mean=${v.mean == null ? '-' : String.format('%.1f', v.mean)} max=${v.max} ${v.hist}" }.join('; ')
println "[1c]   largest before: ${largestBefore}"
println "[1c]   largest after:  ${largestAfter}"
println "[1c]   objects > 50,000 um^2: ${summary.n_objects_over_50000_um2}; annotation hash ${hashBefore} unchanged; ${timings}" +
        (summary.vm_hwm_kb ? String.format(' VmHWM %.1f GB', (summary.vm_hwm_kb as long) / 1048576.0d) : '')
boolean changed = (P.APPLY as boolean) && nP > 0
binding.setVariable('V3_RESULT', [ok: true, changed: changed,
                                  message: "1c ${changed ? '' : '(not applied) '}${classCountsBefore} -> ${classCountsAfter}; dropped ${summary.n_dropped_halo_only}, " +
                                          "merges ${nMergedObjects}, splits ${nSplitPix}",
                                  details: plain(summary)])
