/*
 * 4_islet_expander_v3.groovy   (v3 copy of scripts/4_islet_expander.groovy; the original is unchanged)
 *
 * Generate non-overlapping IsletExpanded bands of width EXPAND_MARGIN_UM um around each Islet, split at the
 * midline between adjacent islets via a densified-boundary Voronoi partition.
 *
 * v3 (1_panc_objects_v3.groovy, drotar_cd99) creates three endocrine classes by mean caliper diameter:
 *   Islet (>= 35 um)       -> "Islet_N",   'Islet: UID',   20 um IsletExpanded band (Heide: watershed-split
 *                             20 um; Bodenmiller: <= 20 um). Bands are made for Islets ONLY.
 *   EndocrineCluster       -> "Cluster_N", 'Cluster: UID', no band
 *   EndocrineSingle        -> "Single_N",  'Single: UID',  no band
 * Each class has its own UID space. A cluster lying inside a band keeps its own cells: script 8 parents
 * clusters and singles before bands.
 *
 * Islet spacing (Drotar 2025): for each Islet
 *   'Nearest islet centroid distance um'  centroid to centroid, among Islets whose centroid lies within
 *                                         +-NEAREST_WINDOW_UM in x and y (a 6 x 6 mm window); NaN if none
 *   'Nearest islet edge distance um'      boundary to boundary, among Islets within NEAREST_WINDOW_UM
 * Only Islet-class objects count as neighbours (clusters and singles do not).
 *
 * v3 also:
 *   - clips every band to Tissue minus (BAND_EXCLUDE_CLASSES = Artifact, pLN): a band never claims a fold,
 *     debris or lymph-node area (script 8 parents Artifact and pLN before IsletExpanded as well);
 *   - writes 'Artifact overlap um^2' (area inside Artifact annotations) on every Islet, EndocrineCluster and
 *     EndocrineSingle, so objects drawn/kept over an Artifact at review can be excluded downstream;
 *   - removes keys that belong to another class after a reclassification: 'Islet: UID', 'Cluster: UID' and
 *     'Single: UID' stay only on their own class ('Islet: UID' also on bands), the nearest-islet distances
 *     only on Islets;
 *   - gives detections that the band removal left directly under the root (bands live at the root) to the
 *     Tissue (else pLN) annotation containing their nucleus centroid.
 *
 * Order: signed distances (script 3) exclude IsletExpanded explicitly, so this may run before or after it.
 * Produces: IsletExpanded annotations with "Islet: UID" measurement matching parent core UID.
 */

import org.locationtech.jts.geom.Coordinate
import org.locationtech.jts.geom.Envelope
import org.locationtech.jts.geom.prep.PreparedGeometryFactory
import org.locationtech.jts.operation.union.UnaryUnionOp
import org.locationtech.jts.geom.GeometryFactory
import org.locationtech.jts.index.strtree.STRtree
import org.locationtech.jts.densify.Densifier
import org.locationtech.jts.triangulate.VoronoiDiagramBuilder
import org.locationtech.jts.operation.buffer.BufferOp
import org.locationtech.jts.operation.buffer.BufferParameters
import org.locationtech.jts.operation.distance.IndexedFacetDistance
import org.locationtech.jts.operation.union.CascadedPolygonUnion
import org.locationtech.jts.simplify.TopologyPreservingSimplifier

import qupath.lib.objects.PathObjects
import qupath.lib.objects.PathObjectTools
import qupath.lib.objects.classes.PathClass
import qupath.lib.roi.GeometryTools

import java.util.concurrent.ConcurrentHashMap
import java.util.concurrent.atomic.AtomicInteger

evaluate(new File(buildFilePath(PROJECT_BASE_DIR, 'scripts', 'lib', 'common.groovy')))

// ==================== PARAMETERS ====================
String CLASS_CLUSTER       = 'EndocrineCluster'
String CLASS_SINGLE        = 'EndocrineSingle'
double NEAREST_WINDOW_UM   = 3000.0    // Drotar: nearest islet among islets within +-3 mm
String BAND_EXCLUDE_CLASSES = 'Artifact,pLN'   // bands are clipped to Tissue minus these classes
// ====================================================

def server = getCurrentServer()
if (server == null) { throw new IllegalStateException("[04v3] no image open") }
def cal = server.getPixelCalibration()
double pxSize = cal?.getAveragedPixelSizeMicrons()
if (pxSize == null || Double.isNaN(pxSize)) {
    println "[04v3] WARN: pixel size unknown; assuming 1.0 um/px"; pxSize = 1.0
}
String UNIT = cal.getPixelWidthUnit()   // 'µm', matching QuPath's own measurement names

double expandPx      = EXPAND_MARGIN_UM   / pxSize
double densifyPx     = DENSIFY_SPACING_UM / pxSize
double simplifyTolPx = Math.max(0.5, expandPx / SIMPLIFY_DIVISOR)
println String.format("[04v3] expand=%.1f um (%.2f px), densify=%.1f um (%.2f px), simpTol=%.2f px",
        EXPAND_MARGIN_UM, expandPx, DENSIFY_SPACING_UM, densifyPx, simplifyTolPx)

def coresAll = getAnnotationObjects().findAll { a ->
    a.getROI() && a.getPathClass()?.toString() == CLASS_ISLET
}
def clusters = getAnnotationObjects().findAll { it.getROI() && it.getPathClass()?.toString() == CLASS_CLUSTER }
def singles  = getAnnotationObjects().findAll { it.getROI() && it.getPathClass()?.toString() == CLASS_SINGLE }
// No early `return`: under run_v3 the provenance footer is appended to this script and must be reached;
// with no Islets the islet_bands postcondition fails and records why.
if (coresAll.isEmpty()) println "[04v3] WARN: no ${CLASS_ISLET} annotations found."
println "[04v3] Found ${coresAll.size()} ${CLASS_ISLET}, ${clusters.size()} ${CLASS_CLUSTER}, ${singles.size()} ${CLASS_SINGLE}."

// Idempotent rerun
def existingBands = getAnnotationObjects().findAll { it.getPathClass()?.toString() == CLASS_ISLET_EXPANDED }
if (!existingBands.isEmpty()) {
    println "[04v3] Removing ${existingBands.size()} existing ${CLASS_ISLET_EXPANDED} annotations."
    removeObjects(existingBands, true)
}

// Keys that belong to another class (a reviewed object may have been reclassified, e.g. Islet -> EndocrineCluster):
// each UID key stays only on its own class ('Islet: UID' also on bands), the nearest-islet distances only on Islets.
def UID_KEYS = ['Islet: UID', 'Cluster: UID', 'Single: UID']
def OWN_UID = [(CLASS_ISLET): 'Islet: UID', (CLASS_CLUSTER): 'Cluster: UID', (CLASS_SINGLE): 'Single: UID',
               (CLASS_ISLET_EXPANDED): 'Islet: UID']
def NEAREST_KEYS = ["Nearest islet centroid distance ${UNIT}".toString(), "Nearest islet edge distance ${UNIT}".toString()]
def staleRemoved = new TreeMap()
getAnnotationObjects().each { a ->
    String c = a.getPathClass()?.toString()
    def ml = a.getMeasurementList()
    def drop = UID_KEYS.findAll { it != OWN_UID[c] && ml.containsKey(it) }
    if (c != CLASS_ISLET) drop += NEAREST_KEYS.findAll { ml.containsKey(it) }
    if (!drop.isEmpty()) {
        ml.removeAll(drop as String[])
        drop.each { k -> staleRemoved["${c}: ${k}".toString()] = (staleRemoved["${c}: ${k}".toString()] ?: 0) + 1 }
    }
    ml.close()
}
if (!staleRemoved.isEmpty()) println "[04v3] Removed keys of another class (reclassified objects): ${staleRemoved}"

// Band clipping area: Tissue minus BAND_EXCLUDE_CLASSES; Artifact overlap on endocrine objects
def unionOf = { List<String> classes ->
    def gs = getAnnotationObjects().findAll { it.getPathClass()?.toString() in classes && it.getROI()?.isArea() }.collect { it.getROI().getGeometry() }
    gs.isEmpty() ? null : UnaryUnionOp.union(gs)
}
def bandExclude = BAND_EXCLUDE_CLASSES.split(',').collect { it.trim() }.findAll { it }
def tissueU = unionOf([CLASS_TISSUE])
def bandExclU = unionOf(bandExclude)
def bandArea = tissueU == null ? null : (bandExclU == null ? tissueU : tissueU.difference(bandExclU))
def bandAreaPrep = bandArea == null ? null : PreparedGeometryFactory.prepare(bandArea)
if (bandArea == null) println "[04v3] WARN: no ${CLASS_TISSUE} annotation; bands are not clipped"
else println String.format("[04v3] bands are clipped to %s minus %s (%.2f mm^2)", CLASS_TISSUE, bandExclude,
        bandArea.getArea() * pxSize * pxSize / 1e6)
def artifactU = unionOf(['Artifact'])
def artifactPrep = artifactU == null ? null : PreparedGeometryFactory.prepare(artifactU)
String kArtifact = "Artifact overlap ${UNIT}^2".toString()
int nOnArtifact = 0
(coresAll + clusters + singles).each { o ->
    double a = 0
    if (artifactPrep != null) {
        def g = o.getROI().getGeometry()
        if (artifactPrep.intersects(g)) a = (artifactPrep.covers(g) ? g.getArea() : g.intersection(artifactU).getArea()) * pxSize * pxSize
    }
    if (a > 0) nOnArtifact++
    putMeas(o, kArtifact, a)
}
println "[04v3] '${kArtifact}' written on ${coresAll.size() + clusters.size() + singles.size()} endocrine objects (${nOnArtifact} overlap an Artifact)"

// Assign stable, UNIQUE UIDs per class.
// - Sort with an explicit comparator: sort { [x, y] } compares the key Lists by hashCode in Groovy
//   (List is not Comparable), so it was not x-then-y order.
// - Two passes: parsed names claim their UIDs first; new UIDs start above the maximum. A UID that is
//   already taken (e.g. two split parts both named "Islet_7", or an unnamed core numbered 1..n while a
//   named core already holds that number) gets a fresh UID and a new name, because duplicate UIDs
//   merge the islets' Voronoi regions and produce overlapping bands.
// - v3: an object carrying ANOTHER class's automatic name (e.g. "Cluster_12" re-classified to Islet during
//   review) is renamed to its own class's pattern; custom names are kept.
def AUTO = [(CLASS_ISLET): 'Islet', (CLASS_CLUSTER): 'Cluster', (CLASS_SINGLE): 'Single']
def parseAuto = { obj, String prefix ->
    def m = ((obj.getName() ?: "") =~ /(?:^|\s)${prefix}[_\s-]*([0-9]+)/)
    m.find() ? (m.group(1) as int) : -1
}
def isForeignAuto = { obj, String prefix ->
    AUTO.values().any { p -> p != prefix && parseAuto(obj, p) > 0 }
}
def assignUIDs = { List objs, String prefix ->
    def ordered = new ArrayList(objs).sort { a, b ->
        a.getROI().getBoundsX() <=> b.getROI().getBoundsX() ?: a.getROI().getBoundsY() <=> b.getROI().getBoundsY()
    }
    def uids = new IdentityHashMap()
    def used = new HashSet<Integer>()
    int renamed = 0
    ordered.each { o ->
        int uid = parseAuto(o, prefix)
        if (uid > 0 && used.add(uid)) uids[o] = uid
    }
    int next = (used.isEmpty() ? 0 : used.max()) + 1
    ordered.each { o ->
        if (uids.containsKey(o)) return
        int uid = next++
        used.add(uid)
        if (parseAuto(o, prefix) > 0 || isForeignAuto(o, prefix)) { o.setName("${prefix}_${uid}"); renamed++ }
        else if (!o.getName()) o.setName("${prefix}_${uid}")    // custom names are kept
        uids[o] = uid
    }
    ordered.each { o -> putMeas(o, "${prefix}: UID", uids[o] as double) }
    if (renamed > 0) println "[04v3] Renamed ${renamed} ${prefix} object(s) whose name duplicated a UID or another class's name."
    uids
}
def uidMap = assignUIDs(coresAll, 'Islet')
assignUIDs(clusters, 'Cluster')
assignUIDs(singles, 'Single')
fireHierarchyUpdate()

// ---- Drotar nearest-islet distances (Islets only) ----
long tN = System.currentTimeMillis()
double windowPx = NEAREST_WINDOW_UM / pxSize
String kCentroid = "Nearest islet centroid distance ${UNIT}"
String kEdge = "Nearest islet edge distance ${UNIT}"
coresAll.groupBy { it.getROI().getImagePlane() }.each { plane, list ->
    def cTree = new STRtree()
    def gTree = new STRtree()
    list.each { o ->
        def r = o.getROI()
        cTree.insert(new Envelope(r.getCentroidX(), r.getCentroidX(), r.getCentroidY(), r.getCentroidY()), o)
        gTree.insert(r.getGeometry().getEnvelopeInternal(), o)
    }
    cTree.build(); gTree.build()
    list.parallelStream().forEach { o ->
        def r = o.getROI()
        double cx = r.getCentroidX(), cy = r.getCentroidY()
        double best = Double.NaN
        for (def other : cTree.query(new Envelope(cx - windowPx, cx + windowPx, cy - windowPx, cy + windowPx))) {
            if (other.is(o)) continue
            double dx = other.getROI().getCentroidX() - cx, dy = other.getROI().getCentroidY() - cy
            if (Math.abs(dx) > windowPx || Math.abs(dy) > windowPx) continue
            double d = Math.hypot(dx, dy)
            if (Double.isNaN(best) || d < best) best = d
        }
        def g = r.getGeometry()
        def genv = new Envelope(g.getEnvelopeInternal())
        genv.expandBy(windowPx)
        def fd = null
        double bestEdge = Double.NaN
        for (def other : gTree.query(genv)) {
            if (other.is(o)) continue
            if (fd == null) fd = new IndexedFacetDistance(g)
            double d = fd.distance(other.getROI().getGeometry())
            if (d <= windowPx && (Double.isNaN(bestEdge) || d < bestEdge)) bestEdge = d
        }
        def ml = o.getMeasurementList()
        ml.put(kCentroid, Double.isNaN(best) ? Double.NaN : best * pxSize)
        ml.put(kEdge, Double.isNaN(bestEdge) ? Double.NaN : bestEdge * pxSize)
        ml.close()
    }
}
int nNoNeighbour = coresAll.count { Double.isNaN(it.getMeasurementList().get(kCentroid)) }
println String.format("[04v3] Nearest-islet distances for %d islets in %.1f s (%d without a neighbour within +-%.0f um)",
        coresAll.size(), (System.currentTimeMillis() - tN) / 1000.0, nNoNeighbour, NEAREST_WINDOW_UM)

println "[04v3] Building Voronoi partition..."
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
    println "[04v3]   plane=${plane}: ${coresOnPlane.size()} cores, ${sites.size()} boundary sites"

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
    if (unmatched > 0) println "[04v3]     matched=${matched}, unmatched=${unmatched}"

    def perUID = [:]
    cellsByUID.each { uid, cells ->
        def unioned = (cells.size() == 1) ? cells[0] : CascadedPolygonUnion.union(cells)
        perUID.put(uid, unioned)
    }
    voronoiByPlaneUID.put(plane, perUID)
}
println String.format("[04v3] Voronoi built in %.1f s", (System.currentTimeMillis() - t0) / 1000.0)

def coreIndexByPlane = [:]
coresByPlane.each { plane, coresOnPlane ->
    def tree = new STRtree()
    coresOnPlane.each { c -> tree.insert(c.getROI().getGeometry().getEnvelopeInternal(), c) }
    tree.build()
    coreIndexByPlane[plane] = tree
}

println "[04v3] Generating bands (parallelStream)..."
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

        // v3: never outside Tissue, never over an Artifact or a pLN
        if (bandAreaPrep != null && !bandAreaPrep.covers(band)) {
            band = band.intersection(bandArea)
            if (band == null || band.isEmpty()) return
        }

        band = toPolygonal(band)
        if (band == null || band.isEmpty()) return

        results.put(core, [geom: band, plane: plane, uid: uid])
    } catch (Throwable t) {
        failed.incrementAndGet()
        println "[04v3] WARN: band gen failed for UID ${uidMap[core]}: ${t.message}"
    }
}
println String.format("[04v3] Band geometry complete in %.1f s (%d failures)",
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

// Detections of the removed bands moved to the bands' parent, the root. Give them to the region containing their
// nucleus centroid (Tissue first, then pLN: the test point script 8 uses); script 8 parents them into the new bands.
def gfR = GeometryTools.getDefaultFactory()
def orphans = getCurrentHierarchy().getRootObject().getChildObjects().findAll { it.isDetection() }
if (!orphans.isEmpty()) {
    def claimed = Collections.newSetFromMap(new IdentityHashMap())
    def moved = new LinkedHashMap()
    [CLASS_TISSUE, 'pLN'].each { rc ->
        getAnnotationObjects().findAll { it.getPathClass()?.toString() == rc && it.getROI()?.isArea() }.each { region ->
            def prep = PreparedGeometryFactory.prepare(region.getROI().getGeometry())
            def kids = orphans.findAll { d -> if (claimed.contains(d)) return false
                def r = PathObjectTools.getROI(d, true)
                prep.contains(gfR.createPoint(new Coordinate(r.getCentroidX(), r.getCentroidY()))) }
            if (!kids.isEmpty()) { region.addChildObjects(kids); claimed.addAll(kids); moved[rc] = (moved[rc] ?: 0) + kids.size() }
        }
    }
    println "[04v3] Re-parented detections left under the root: ${moved}; ${orphans.size() - claimed.size()} lie in no Tissue/pLN and stay"
}
fireHierarchyUpdate()
println "[04v3] Created ${newAnns.size()} ${CLASS_ISLET_EXPANDED} annotations for ${coresAll.size()} ${CLASS_ISLET}; " +
        "named ${clusters.size()} clusters and ${singles.size()} singles. Done."
