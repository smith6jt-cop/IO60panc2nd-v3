/*
 * 4_islet_expander.groovy
 *
 * Generate non-overlapping IsletExpanded bands of width EXPAND_MARGIN_UM um around each Islet, split at the midline
 * between adjacent islets via a densified-boundary Voronoi partition.
 *
 * Must run AFTER Signed Distance to Annotations so that all cells have "Signed distance to Islet"
 * measurements computed against the islet cores only, before IsletExpanded
 * annotations exist in the hierarchy.
 *
 * Produces: IsletExpanded annotations with "Islet: UID" measurement matching
 *           parent core UID.
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

def server = getCurrentServer()
if (server == null) { println "[03] ERROR: no image open."; return }
def cal = server.getPixelCalibration()
double pxSize = cal?.getAveragedPixelSizeMicrons()
if (pxSize == null || Double.isNaN(pxSize)) {
    println "[03] WARN: pixel size unknown; assuming 1.0 um/px"; pxSize = 1.0
}

double expandPx      = EXPAND_MARGIN_UM   / pxSize
double densifyPx     = DENSIFY_SPACING_UM / pxSize
double simplifyTolPx = Math.max(0.5, expandPx / SIMPLIFY_DIVISOR)
println String.format("[03] expand=%.1f um (%.2f px), densify=%.1f um (%.2f px), simpTol=%.2f px",
        EXPAND_MARGIN_UM, expandPx, DENSIFY_SPACING_UM, densifyPx, simplifyTolPx)

def coresAll = getAnnotationObjects().findAll { a ->
    a.getROI() && a.getPathClass()?.toString() == CLASS_ISLET
}
if (coresAll.isEmpty()) { println "[03] No ${CLASS_ISLET} annotations found."; return }
println "[03] Found ${coresAll.size()} ${CLASS_ISLET} cores."

// Idempotent rerun
def existingBands = getAnnotationObjects().findAll { it.getPathClass()?.toString() == CLASS_ISLET_EXPANDED }
if (!existingBands.isEmpty()) {
    println "[03] Removing ${existingBands.size()} existing ${CLASS_ISLET_EXPANDED} annotations."
    removeObjects(existingBands, true)
}

// Assign stable, UNIQUE UIDs.
// - Sort with an explicit comparator: sort { [x, y] } compares the key Lists by hashCode in Groovy
//   (List is not Comparable), so it was not x-then-y order.
// - Two passes: parsed names claim their UIDs first; new UIDs start above the maximum. A UID that is
//   already taken (e.g. two split parts both named "Islet_7", or an unnamed core numbered 1..n while a
//   named core already holds that number) gets a fresh UID and a new name, because duplicate UIDs
//   merge the islets' Voronoi regions and produce overlapping bands.
def coresOrdered = new ArrayList(coresAll).sort { a, b ->
    a.getROI().getBoundsX() <=> b.getROI().getBoundsX() ?: a.getROI().getBoundsY() <=> b.getROI().getBoundsY()
}
def parseUID = { core ->
    def m = ((core.getName() ?: "") =~ /(?:^|\s)Islet[_\s-]*([0-9]+)/)
    m.find() ? (m.group(1) as int) : -1
}
def uidMap = new IdentityHashMap()
def usedUIDs = new HashSet<Integer>()
int nRenamed = 0
coresOrdered.each { core ->
    int uid = parseUID(core)
    if (uid > 0 && usedUIDs.add(uid)) uidMap[core] = uid
}
int nextUID = (usedUIDs.isEmpty() ? 0 : usedUIDs.max()) + 1
coresOrdered.each { core ->
    if (uidMap.containsKey(core)) return
    int uid = nextUID++
    usedUIDs.add(uid)
    if (parseUID(core) > 0) { core.setName("Islet_${uid}"); nRenamed++ }   // duplicate Islet_N name
    else if (!core.getName()) core.setName("Islet_${uid}")               // custom names are kept
    uidMap[core] = uid
}
if (nRenamed > 0) println "[03] Renamed ${nRenamed} core(s) whose Islet_N name duplicated another core."
coresOrdered.each { core -> putMeas(core, "Islet: UID", uidMap[core] as double) }
fireHierarchyUpdate()

println "[03] Building Voronoi partition..."
long t0 = System.currentTimeMillis()

def gf = new GeometryFactory()
def coresByPlane = coresAll.groupBy { it.getROI().getImagePlane() }
def voronoiByPlaneUID = [:]

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
    println "[03]   plane=${plane}: ${coresOnPlane.size()} cores, ${sites.size()} boundary sites"

    def vBuilder = new VoronoiDiagramBuilder()
    vBuilder.setSites(sites)
    double margin = expandPx * 3.0
    vBuilder.setClipEnvelope(new Envelope(
            -margin, server.getWidth()  + margin,
            -margin, server.getHeight() + margin))
    def diagram = vBuilder.getDiagram(gf)

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
    if (unmatched > 0) println "[03]     matched=${matched}, unmatched=${unmatched}"

    def perUID = [:]
    cellsByUID.each { uid, cells ->
        def unioned = (cells.size() == 1) ? cells[0] : CascadedPolygonUnion.union(cells)
        perUID.put(uid, unioned)
    }
    voronoiByPlaneUID.put(plane, perUID)
}
println String.format("[03] Voronoi built in %.1f s", (System.currentTimeMillis() - t0) / 1000.0)

def coreIndexByPlane = [:]
coresByPlane.each { plane, coresOnPlane ->
    def tree = new STRtree()
    coresOnPlane.each { c -> tree.insert(c.getROI().getGeometry().getEnvelopeInternal(), c) }
    tree.build()
    coreIndexByPlane[plane] = tree
}

println "[03] Generating bands (parallelStream)..."
long t1 = System.currentTimeMillis()

// Round joins: every point of the expanded outline is exactly expandPx from the core. JOIN_MITRE
// reached sqrt(2)*d at right-angle corners and up to mitreLimit*d (5 x 20 um) at acute tips.
def bufParams = new BufferParameters(8, BufferParameters.CAP_ROUND, BufferParameters.JOIN_ROUND, BufferParameters.DEFAULT_MITRE_LIMIT)
def results = new ConcurrentHashMap()
def failed  = new AtomicInteger(0)

coresAll.parallelStream().forEach { core ->
    try {
        def coreROI  = core.getROI()
        def plane    = coreROI.getImagePlane()
        def coreGeom = coreROI.getGeometry()
        if (coreGeom == null || coreGeom.isEmpty()) return
        int uid = (int) uidMap[core]

        def expanded = new BufferOp(coreGeom, bufParams).getResultGeometry(expandPx)
        if (expanded == null || expanded.isEmpty()) return
        // Simplify only the outer envelope, BEFORE the exact clipping steps below, so simplification
        // cannot push the band back into its own core, another core, or a neighbour's Voronoi region.
        expanded = TopologyPreservingSimplifier.simplify(expanded.buffer(0), simplifyTolPx)

        def band = expanded.difference(coreGeom)
        if (band == null || band.isEmpty()) return

        def vRegion = voronoiByPlaneUID.get(plane)?.get(uid)
        if (vRegion != null && !vRegion.isEmpty()) {
            band = band.intersection(vRegion)
            if (band == null || band.isEmpty()) return
        }

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

        band = toPolygonal(band)
        if (band == null || band.isEmpty()) return

        results.put(core, [geom: band, plane: plane, uid: uid])
    } catch (Throwable t) {
        failed.incrementAndGet()
        println "[03] WARN: band gen failed for UID ${uidMap[core]}: ${t.message}"
    }
}
println String.format("[03] Band geometry complete in %.1f s (%d failures)",
        (System.currentTimeMillis() - t1) / 1000.0, failed.get())

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
println "[03] Created ${newAnns.size()} ${CLASS_ISLET_EXPANDED} annotations. Done."
