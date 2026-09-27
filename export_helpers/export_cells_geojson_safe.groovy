/*
 * export_cells_geojson_safe.groovy
 *
 * Robust variant of export_cells_geojson.groovy for images that fail the batch
 * exporter with a GSON "Nesting problem" (a malformed cell geometry breaks QuPath's
 * atomic exportObjectsToGeoJson — seen on donor 6623). Same output/format, but it
 * FILTERS OUT invalid/empty/non-polygon cell geometries first (and reports the count),
 * so a handful of degenerate segmentation objects can't kill the whole export. Dropped
 * cells are simply absent from the mask -> they get no REDSEA spillover correction
 * (acceptable; REDSEA already tolerates a few missing/NaN edge cells, joined by object_id).
 *
 * Output (one file per image): <OUT_DIR>/cells__<imageTag>.geojson  (id = object_id UUID)
 * Run headless with the QuPath GUI closed, same invocation as export_cells_geojson.groovy.
 */

String DEFAULT_OUT_DIR =
    '/home/smith6jt/IO60panc2nd/Islet-Explorer-Senior/data/redsea_scratch/geojson'
String outDir = (args != null && args.size() > 0 && args[0]?.trim()) ? args[0].trim()
                                                                      : DEFAULT_OUT_DIR

def server = getCurrentServer()
if (server == null) { println "ERROR: no image open."; return }

def all = getDetectionObjects()
def cells = all.findAll {
    def r = it.getROI()
    if (r == null || r.isEmpty() || !r.isArea()) return false
    def g
    try { g = r.getGeometry() } catch (Throwable t) { return false }
    if (g == null || g.isEmpty() || !g.isValid()) return false
    def gt = g.getGeometryType()
    return gt == 'Polygon' || gt == 'MultiPolygon'
}
int dropped = all.size() - cells.size()
println "[safe-export] kept ${cells.size()} / ${all.size()} cells (dropped ${dropped} invalid/empty/non-polygon)"
if (cells.isEmpty()) { println "ERROR: no valid cell geometries."; return }

mkdirs(outDir)
def imageTag = (getProjectEntry()?.getImageName() ?: server.getMetadata()?.getName() ?: 'image')
        .replaceAll('[^A-Za-z0-9._-]', '_')
def outPath = buildFilePath(outDir, "cells__${imageTag}.geojson")

println "[safe-export] ${cells.size()} cells -> ${outPath}"
long t0 = System.currentTimeMillis()
exportObjectsToGeoJson(cells, outPath, "FEATURE_COLLECTION", "EXCLUDE_MEASUREMENTS")
def sizeMB = new File(outPath).length() / (1024.0 * 1024.0)
println String.format("[safe-export] Wrote %s (%.1f MB) in %.1fs",
        outPath, sizeMB, (System.currentTimeMillis() - t0) / 1000.0)
println "[safe-export] Done."
