/*
 * export_cells_geojson_v3.groovy -- per-cell polygons (GeoJSON) for pixel-level REDSEA and the
 * Parquet join. Fixed copy of scripts/export_helpers/export_cells_geojson.groovy (left untouched).
 *
 * Each feature carries "id" = the QuPath detection UUID (== "Object ID" in the cells Parquet) and
 * the full-resolution pixel polygon; measurements are excluded (they are in the Parquet).
 *
 * Fixes relative to the original:
 *   - duplicate check: the original built its buckets with `[:].withDefault { [] }` and then looked up
 *     neighbours with buckets[[bx+dx, by+dy]] while iterating the same map -- withDefault's get()
 *     INSERTS missing keys, so the lookup mutated the map mid-iteration (ConcurrentModificationException
 *     risk). Now: computeIfAbsent to fill, getOrDefault to look up, iteration over a snapshot.
 *   - pair ordering: the original skipped a pair when identityHashCode(b) <= identityHashCode(a);
 *     identity hashes are not unique, so two objects with equal hashes were never compared and a
 *     duplicate pair could hide. Now each cell has a unique list index and pairs use ia < ib.
 *   - duplicates are a failure (v3 script 0 guarantees none), not a warning followed by an export.
 *   - the output is written to <path>.part, stream-verified (every detection UUID exactly once, no
 *     extra IDs) and only then renamed; the verdict goes to V3_RESULT for the run_v3 footer.
 *
 * Output path: args[0]. If it ends in .geojson it is the file; otherwise a directory, and the file is
 * <dir>/cells__<imageTag>.geojson as before. Default: the original's REDSEA scratch directory.
 * Headless (read-only; never saves):
 *   QuPath script -p <project> -i "<image>" -a <out.geojson> export_cells_geojson_v3.groovy
 */

// ==================== PARAMETERS ====================
String DEFAULT_OUT_DIR = '/home/smith6jt/IO60panc2nd/Islet-Explorer-Senior/data/redsea_scratch/geojson'
double DUP_RADIUS_UM = 1.0     // same rule as 0_cell_detection_panc.groovy and phenocycler.cell_qc
double DUP_AREA_TOL  = 0.15
boolean FAIL_ON_DUPLICATES = true
// ====================================================

String target = (binding.hasVariable('args') && args && args.length > 0 && args[0]?.trim()) ? args[0].trim() : DEFAULT_OUT_DIR
def server = getCurrentServer()
if (server == null) throw new IllegalStateException('no image open')
def imageTag = (getProjectEntry()?.getImageName() ?: server.getMetadata()?.getName() ?: 'image').replaceAll('[^A-Za-z0-9._-]', '_')
def outFile = target.endsWith('.geojson') ? new File(target) : new File(target, "cells__${imageTag}.geojson")
def partFile = new File(outFile.getPath() + '.part')
outFile.getParentFile().mkdirs()
partFile.delete()

def cells = new ArrayList(getDetectionObjects().findAll { it.getROI() != null })
def details = [cells: cells.size(), output: outFile.getPath()]
if (cells.isEmpty()) {
    V3_RESULT = [ok: false, changed: false, message: 'no detections with a ROI', details: details]
    println '[export_cells_geojson] ERROR: no detections'
} else {
    // ---- duplicate-detection guard (grid buckets of the radius: O(n)) ----
    def cal = server.getPixelCalibration()
    double pxSize = cal.hasPixelSizeMicrons() ? cal.getAveragedPixelSizeMicrons() : 1.0
    double rPx = DUP_RADIUS_UM / pxSize
    def buckets = new HashMap()
    for (int i = 0; i < cells.size(); i++) {
        def roi = cells[i].getROI()
        buckets.computeIfAbsent([(int) Math.floor(roi.getCentroidX() / rPx), (int) Math.floor(roi.getCentroidY() / rPx)]) { [] } << i
    }
    long dupPairs = 0
    for (entry in new ArrayList(buckets.entrySet())) {
        int bx = entry.key[0], by = entry.key[1]
        def near = []
        for (dx in -1..1) for (dy in -1..1) near.addAll(buckets.getOrDefault([bx + dx, by + dy], []))
        for (int ia in entry.value) {
            def ra = cells[ia].getROI()
            for (int ib in near) {
                if (ib <= ia) continue            // unique indices: every unordered pair exactly once
                def rb = cells[ib].getROI()
                if (Math.hypot(ra.getCentroidX() - rb.getCentroidX(), ra.getCentroidY() - rb.getCentroidY()) > rPx) continue
                double aa = ra.getArea(), ab = rb.getArea()
                if (Math.abs(aa - ab) / Math.max(Math.max(aa, ab), 1e-9) < DUP_AREA_TOL) dupPairs++
            }
        }
    }
    details.duplicate_pairs = dupPairs
    println "[export_cells_geojson] duplicate check: ${dupPairs} pair(s) within ${DUP_RADIUS_UM} um"

    if (dupPairs > 0 && FAIL_ON_DUPLICATES) {
        V3_RESULT = [ok: false, changed: false, details: details,
                     message: "${dupPairs} duplicate detection pairs; re-run script 0 (nothing exported)".toString()]
    } else {
        long t0 = System.currentTimeMillis()
        exportObjectsToGeoJson(cells, partFile.getPath(), 'FEATURE_COLLECTION', 'EXCLUDE_MEASUREMENTS')
        details.export_seconds = (System.currentTimeMillis() - t0) / 1000.0

        // ---- verify by streaming the file: every UUID exactly once, nothing else ----
        def want = new HashSet(cells.size() * 2)
        cells.each { want.add(it.getID().toString()) }
        def seen = new HashSet(cells.size() * 2)
        long nFeatures = 0, nRepeated = 0, nUnknown = 0, nNoId = 0
        def reader = new com.google.gson.stream.JsonReader(new InputStreamReader(new BufferedInputStream(new FileInputStream(partFile), 1 << 20), 'UTF-8'))
        try {
            reader.beginObject()
            while (reader.hasNext()) {
                if (reader.nextName() != 'features') { reader.skipValue(); continue }
                reader.beginArray()
                while (reader.hasNext()) {
                    reader.beginObject()
                    String id = null
                    while (reader.hasNext()) {
                        if (reader.nextName() == 'id') id = reader.nextString() else reader.skipValue()
                    }
                    reader.endObject()
                    nFeatures++
                    if (id == null) nNoId++
                    else if (!want.contains(id)) nUnknown++
                    else if (!seen.add(id)) nRepeated++
                }
                reader.endArray()
            }
            reader.endObject()
        } finally {
            reader.close()
        }
        long nMissing = want.size() - seen.size()
        details += [features: nFeatures, missing_ids: nMissing, repeated_ids: nRepeated, unknown_ids: nUnknown,
                    features_without_id: nNoId, bytes: partFile.length()]
        boolean ok = nFeatures == cells.size() && nMissing == 0 && nRepeated == 0 && nUnknown == 0 && nNoId == 0
        if (ok) {
            java.nio.file.Files.move(partFile.toPath(), outFile.toPath(), java.nio.file.StandardCopyOption.REPLACE_EXISTING)
            details.outputs = [outFile.getPath()]
        }
        V3_RESULT = [ok: ok, changed: false, details: details,
                     message: "${nFeatures} features for ${cells.size()} cells (missing ${nMissing}, repeated ${nRepeated}, unknown ${nUnknown})".toString()]
    }
    println "[export_cells_geojson] ${V3_RESULT.ok ? 'OK' : 'FAILED'}: ${V3_RESULT.message} -> ${outFile}"
}
