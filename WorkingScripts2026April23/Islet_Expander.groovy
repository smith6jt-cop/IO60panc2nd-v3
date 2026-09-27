/*
 * 02_islet_bands.groovy
 *
 * Assigns stable UIDs to Islet annotations and generates non-overlapping
 * IsletExpanded bands of width EXPAND_MARGIN_UM um. Overlap between adjacent
 * islets is split at the midline via a densified-boundary Voronoi partition:
 *   - each core boundary is sampled at DENSIFY_SPACING_UM um,
 *   - one global Voronoi diagram is built over all sample points,
 *   - cells are merged by their source UID to form a per-core Voronoi region,
 *   - each core's raw band is intersected with its Voronoi region.
 *
 * Result: bands that meet at the equidistant locus between adjacent cores.
 *
 * Requires: Islet annotations present (stage 01 output).
 * Produces: IsletExpanded annotations, one per Islet core, each carrying
 *           an "Islet: UID" measurement matching its parent core.
 */

evaluate(new File(buildFilePath(PROJECT_BASE_DIR, 'scripts', 'lib', 'common.groovy')))

import org.locationtech.jts.geom.Coordinate
import org.locationtech.jts.geom.Envelope
import org.locationtech.jts.geom.GeometryFactory
import org.locationtech.jts.index.strtree.STRtree
import org.locationtech.jts.densify.Densifier
import org.locationtech.jts.triangulate.VoronoiDiagramBuilder
import org.locationtech.jts.operation.buffer.BufferOp
import org.locationtech.jts.operation.buffer.BufferParameters
import org.locationtech.jts.operation.union.CascadedPolygonUnion
import org.locationtech.jts.simplify.TopologyPreservingSimplifier

import qupath.lib.objects.PathObjects
import qupath.lib.objects.classes.PathClass
import qupath.lib.roi.GeometryTools

import java.util.concurrent.ConcurrentHashMap
import java.util.concurrent.atomic.AtomicInteger

// ---------- Calibration ----------
def server = getCurrentServer()
if (server == null) { println "ERROR: No image open."; return }
def cal = server.getPixelCalibration()
double pxSize = cal?.getAveragedPixelSizeMicrons()
if (pxSize == null || Double.isNaN(pxSize)) {
    println "WARN: pixel size unknown; assuming 1.0 um/px"; pxSize = 1.0
}

double expandPx      = EXPAND_MARGIN_UM   / pxSize
double densifyPx     = DENSIFY_SPACING_UM / pxSize
double simplifyTolPx = Math.max(0.5, expandPx / SIMPLIFY_DIVISOR)
println String.format("[02] expand=%.1f um (%.2f px), densify=%.1f um (%.2f px), simpTol=%.2f px",
        EXPAND_MARGIN_UM, expandPx, DENSIFY_SPACING_UM, densifyPx, simplifyTolPx)

// ---------- Collect cores ----------
def coresAll = getAnnotationObjects().findAll { a ->
    a.getROI() && a.getPathClass()?.toString() == CLASS_ISLET
}
if (coresAll.isEmpty()) { println "No '${CLASS_ISLET}' annotations found."; return }
println "[02] Found ${coresAll.size()} '${CLASS_ISLET}' cores."

// Re-run safety: remove any existing bands so the script is idempotent.
def existingBands = getAnnotationObjects().findAll { it.getPathClass()?.toString() == CLASS_ISLET_EXPANDED }
if (!existingBands.isEmpty()) {
    println "[02] Removing ${existingBands.size()} existing '${CLASS_ISLET_EXPANDED}' annotations."
    removeObjects(existingBands, true)
}

// ---------- Assign stable UIDs ----------
// Parse UID from an existing name if it matches /Islet[_\s-]*([0-9]+)/, otherwise
// assign a fresh sequential UID ordered by (x, y) of the core's bounding box.
def coresOrdered = coresAll.sort { a -> [a.getROI().getBoundsX(), a.getROI().getBoundsY()] }
int nextUID = 1
def uidMap = new IdentityHashMap()
coresOrdered.each { core ->
    int uid = -1
    def nm = core.getName() ?: ""
    def m = (nm =~ /(?:^|\s)Islet[_\s-]*([0-9]+)/)
    if (m.find()) uid = (m.group(1) as int)
    if (uid <= 0) { uid = nextUID++; if (!core.getName()) core.setName("Islet_${uid}") }
    putMeas(core, "Islet: UID", uid as double)
    uidMap[core] = uid
}
fireHierarchyUpdate()

// ---------- Voronoi partition from densified core boundaries ----------
println "[02] Building Voronoi partition..."
long t0 = System.currentTimeMillis()

def gf = new GeometryFactory()
def coresByPlane = coresAll.groupBy { it.getROI().getImagePlane() }
def voronoiByPlaneUID = [:]   // plane -> (uid -> Voronoi region Geometry)

// Quantize floating-point coordinates to a stable long key so we can look up
// the source UID from the Coordinate JTS attaches to each Voronoi cell
// (Polygon.getUserData()). 0.001 px resolution is far below any densification
// spacing, so no two distinct sites collide.
def quantize = { double x, double y ->
    ((long) Math.round(x * 1000.0)) * 10_000_000_000L + ((long) Math.round(y * 1000.0))
}

coresByPlane.each { plane, coresOnPlane ->
    def sites    = new ArrayList<Coordinate>()
    def keyToUID = new HashMap<Long, Integer>()

    coresOnPlane.each { core ->
        int uid = (int) uidMap[core]
        def geom = core.getROI().getGeometry()
        if (geom == null || geom.isEmpty()) return
        // Densifier.densify preserves geometry type and adds vertices so no
        // segment exceeds densifyPx. All rings of all polygons contribute.
        def densified = Densifier.densify(geom, densifyPx)
        densified.getCoordinates().each { c ->
            long k = quantize(c.x, c.y)
            if (!keyToUID.containsKey(k)) {
                sites.add(c)
                keyToUID.put(k, uid)
            }
        }
    }

    if (sites.isEmpty()) return
    println "[02]   plane=${plane}: ${coresOnPlane.size()} cores, ${sites.size()} boundary sites"

    def vBuilder = new VoronoiDiagramBuilder()
    vBuilder.setSites(sites)
    // Clip slightly outside the image so boundary cores still get full cells.
    double margin = expandPx * 3.0
    vBuilder.setClipEnvelope(new Envelope(
            -margin, server.getWidth()  + margin,
            -margin, server.getHeight() + margin))
    def diagram = vBuilder.getDiagram(gf)

    // Group Voronoi cells by source UID. JTS sets each cell's userData to the
    // Coordinate of its generating site; we look that Coordinate up in keyToUID.
    def cellsByUID = [:]
    int n = diagram.getNumGeometries()
    int matched = 0, unmatched = 0
    for (int i = 0; i < n; i++) {
        def cell = diagram.getGeometryN(i)
        def sd = cell.getUserData()
        if (!(sd instanceof Coordinate)) { unmatched++; continue }
        Integer uid = keyToUID.get(quantize(sd.x, sd.y))
        if (uid == null) { unmatched++; continue }
        def lst = cellsByUID.get(uid)
        if (lst == null) { lst = []; cellsByUID.put(uid, lst) }
        lst.add(cell)
        matched++
    }
    if (unmatched > 0)
        println "[02]     matched=${matched}, unmatched=${unmatched}"

    // Union the cells belonging to each UID into a single region per core.
    def perUID = [:]
    cellsByUID.each { uid, cells ->
        def unioned = (cells.size() == 1) ? cells[0] : CascadedPolygonUnion.union(cells)
        perUID.put(uid, unioned)
    }
    voronoiByPlaneUID.put(plane, perUID)
}
println String.format("[02] Voronoi built in %.1f s", (System.currentTimeMillis() - t0) / 1000.0)

// ---------- STRtree of cores per plane (safety net) ----------
// Even with good densification, Voronoi approximation may leave thin slivers
// that overlap a neighbouring core. A localized STRtree query + difference
// removes them cheaply.
def coreIndexByPlane = [:]
coresByPlane.each { plane, coresOnPlane ->
    def tree = new STRtree()
    coresOnPlane.each { c -> tree.insert(c.getROI().getGeometry().getEnvelopeInternal(), c) }
    tree.build()
    coreIndexByPlane[plane] = tree
}

// ---------- Generate bands in parallel ----------
println "[02] Generating bands (parallelStream)..."
long t1 = System.currentTimeMillis()

def bufParams = new BufferParameters(8, BufferParameters.CAP_FLAT, BufferParameters.JOIN_MITRE, 5.0)
def results = new ConcurrentHashMap()
def failed  = new AtomicInteger(0)

coresAll.parallelStream().forEach { core ->
    try {
        def coreROI  = core.getROI()
        def plane    = coreROI.getImagePlane()
        def coreGeom = coreROI.getGeometry()
        if (coreGeom == null || coreGeom.isEmpty()) return
        int uid = (int) uidMap[core]

        // Expand the core by the band thickness.
        def expanded = new BufferOp(coreGeom, bufParams).getResultGeometry(expandPx)
        if (expanded == null || expanded.isEmpty()) return
        expanded = TopologyPreservingSimplifier.simplify(expanded.buffer(0), simplifyTolPx)

        // Band = expanded \ own core.
        def band = expanded.difference(coreGeom)
        if (band == null || band.isEmpty()) return

        // Split overlap at midline: intersect with this UID's Voronoi region.
        def vRegion = voronoiByPlaneUID.get(plane)?.get(uid)
        if (vRegion != null && !vRegion.isEmpty()) {
            band = band.intersection(vRegion)
            if (band == null || band.isEmpty()) return
        }

        // Safety: subtract any other nearby core (STRtree-local).
        def tree = coreIndexByPlane[plane]
        def nearby = tree.query(band.getEnvelopeInternal())
        for (def other : nearby) {
            if (other.is(core)) continue
            def og = other.getROI().getGeometry()
            if (og.getEnvelopeInternal().intersects(band.getEnvelopeInternal())) {
                band = band.difference(og)
                if (band.isEmpty()) return
            }
        }

        // Final cleanup.
        band = toPolygonal(band)
        if (band == null || band.isEmpty()) return
        band = TopologyPreservingSimplifier.simplify(band.buffer(0), simplifyTolPx)
        band = toPolygonal(band)
        if (band == null || band.isEmpty()) return

        results.put(core, [geom: band, plane: plane, uid: uid])
    } catch (Throwable t) {
        failed.incrementAndGet()
        println "WARN: band gen failed for UID ${uidMap[core]}: ${t.message}"
    }
}
println String.format("[02] Band geometry complete in %.1f s (%d failures)",
        (System.currentTimeMillis() - t1) / 1000.0, failed.get())

// ---------- Create annotations (sequential; QuPath hierarchy is not thread-safe) ----------
def newAnns = []
coresAll.each { core ->
    def r = results.get(core)
    if (r == null) return
    def roi = GeometryTools.geometryToROI(r.geom, r.plane)
    def ann = PathObjects.createAnnotationObject(roi, PathClass.fromString(CLASS_ISLET_EXPANDED))
    int uid = r.uid
    ann.setName("Islet_${uid}_exp${(int) EXPAND_MARGIN_UM}um")
    putMeas(ann, "Islet: UID", uid as double)
    newAnns << ann
}
addObjects(newAnns)
fireHierarchyUpdate()
println "[02] Created ${newAnns.size()} '${CLASS_ISLET_EXPANDED}' annotations. Done."