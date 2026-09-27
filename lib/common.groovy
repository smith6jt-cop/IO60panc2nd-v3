/*
 * lib/common.groovy
 * Shared helpers and constants for the QuPath pancreas pipeline.
 *
 * Load at the top of each stage with:
 *     evaluate(new File(buildFilePath(PROJECT_BASE_DIR, 'scripts', 'lib', 'common.groovy')))
 *
 * All names assigned below (without `def`) become binding variables in the
 * caller script. The helpers remain in-scope because they are assigned to
 * script-level closures, not def'd locals.
 */

import org.locationtech.jts.geom.Geometry
import org.locationtech.jts.geom.util.PolygonExtracter
import org.locationtech.jts.operation.union.UnaryUnionOp
import org.locationtech.jts.operation.union.CascadedPolygonUnion
import org.locationtech.jts.geom.util.GeometryCombiner
import org.locationtech.jts.simplify.TopologyPreservingSimplifier

// ---------- PathClass names ----------
CLASS_TISSUE         = 'Tissue'
CLASS_ISLET          = 'Islet'
CLASS_ISLET_EXPANDED = 'IsletExpanded'
CLASS_NERVE          = 'Nerve'
CLASS_LYMPHATIC      = 'Lymphatic'
CLASS_BLOOD_VESSEL   = 'Blood_vessel'
CLASS_TISSUE_SAMPLE  = 'TissueSample'

// Legacy, retained for any back-compatible scripts. New code should not use.
CLASS_CAPILLARY      = 'Capillary'

// ---------- Pipeline parameters ----------
EDGE_MARGIN_UM       = 10.0    // stage 00: edge exclusion margin for inner-area clip
EXPAND_MARGIN_UM     = 20.0    // stage 01: islet band thickness
DENSIFY_SPACING_UM   = 5.0     // stage 01: boundary sample spacing for Voronoi
SIMPLIFY_DIVISOR     = 6.0     // stage 01: simplify tolerance = expandPx / divisor
TISSUE_SAMPLE_UM     = 200.0   // sampling square size (if used downstream)
TISSUE_SAMPLE_COUNT  = 5
RNG_SEED             = 42L

// ---------- Channels used by the object classifier ----------
// Stage 05 writes (norm) columns for each of these. Add CD68 alongside the
// three originally requested channels because the composite classifier
// "Iba1CD68_2D" reads Cell: Iba1: Mean (norm) AND Cell: CD68: Mean (norm).
OBJECT_CLF_CHANNELS  = ['M2Gal3', 'CD34', 'Iba1', 'CD68', 'B3TUBB', 'Podoplanin']

// ---------- Geometry helpers ----------

// Keep only polygonal parts, repair topology, union polygons.
// Returns null if input is empty or has no polygonal parts.
toPolygonal = { Geometry g ->
    if (g == null || g.isEmpty()) return null
    def fixed = g.buffer(0)
    def polys = []
    PolygonExtracter.getPolygons(fixed, polys)
    if (polys.isEmpty()) return null
    return new UnaryUnionOp(polys).union()
}

// Light cleanup: fix topology, simplify, fix again.
cleanPoly = { Geometry g, double simpTolPx = 0.2 ->
    if (g == null || g.isEmpty()) return null
    try {
        def a = g.buffer(0)
        def b = TopologyPreservingSimplifier.simplify(a, simpTolPx)
        return b.buffer(0)
    } catch (Throwable t) {
        try { return g.buffer(0) } catch (Throwable t2) { return null }
    }
}

// Robust union of many polygons. Prefers CascadedPolygonUnion (fastest for
// many adjacent polygons) and falls back to GeometryCombiner on topology errors.
safeUnionPolys = { List<Geometry> polys ->
    def ls = polys.collect { toPolygonal(it) }.findAll { it != null && !it.isEmpty() }
    if (ls.isEmpty()) return null
    try {
        return CascadedPolygonUnion.union(ls)
    } catch (Throwable t) {
        try { return toPolygonal(GeometryCombiner.combine(ls).buffer(0)) }
        catch (Throwable t2) {
            println "WARN: safeUnionPolys failed (${ls.size()} polys)"
            return null
        }
    }
}

// Robust difference with topology repair.
safeDiff = { Geometry a, Geometry b ->
    def aa = cleanPoly(a); def bb = cleanPoly(b)
    if (aa == null || aa.isEmpty()) return null
    if (bb == null || bb.isEmpty()) return aa
    try { return aa.difference(bb).buffer(0) }
    catch (Throwable t) {
        try { return aa.buffer(0).difference(bb.buffer(0)).buffer(0) }
        catch (Throwable t2) {
            println "WARN: safeDiff failed"
            return null
        }
    }
}

// Measurement writer. Uses the public `pathObject.measurements` map view,
// which is the documented QuPath scripting API since v0.4 and is the only
// path verified to write reliably across all 0.4-0.7 builds. Earlier
// implementations using getMeasurementList().put(...) silently failed on
// some installs because the returned list view did not always write back
// to the underlying persistent store.
putMeas = { obj, String key, double val ->
    obj.measurements[key] = val
}

// Measurement reader matching the writer. Returns Double.NaN when missing.
getMeas = { obj, String key ->
    def v = obj.measurements[key]
    return (v == null) ? Double.NaN : (v as double)
}

println "[common] Loaded: expand=${EXPAND_MARGIN_UM}um, densify=${DENSIFY_SPACING_UM}um, edge=${EDGE_MARGIN_UM}um"
println "[common] Object-classifier channels: ${OBJECT_CLF_CHANNELS}"