/*
 * exports.groovy -- read-only exports for run_v3.sh (passes C and J). Never saves the image.
 *
 *   args[0] = mode, args[1] = output path
 *   cells        every detection measurement (saveDetectionMeasurements) EXCEPT 'Signed distance…'
 *                (step 3's columns, exported by 'context') -> TSV on /scratch; v3ctl.py post converts it
 *                to Parquet, checks rows, Object IDs and 'Cell: DAPI: Mean', deletes it. The column set is
 *                therefore the same whether this export runs before or after step 3.
 *   context      Object ID, Parent, centroids and every 'Signed distance…' column -> TSV -> Parquet
 *   annotations  annotation measurement table (saveAnnotationMeasurements) -> TSV, kept
 *   islets       every annotation except Tissue and pLN as a GeoJSON FeatureCollection
 *                (object IDs, names, classes, measurements), re-read to verify
 *
 * Each file is written as <path>.part and renamed only after it passed its check, so a file with the
 * final name is always complete. The verdict and the files go to V3_RESULT for the run_v3 footer.
 */

import qupath.lib.io.PathIO

if (args == null || args.length < 2) throw new IllegalArgumentException('usage: exports.groovy <mode> <output path>')
String mode = args[0]
def outFile = new File(args[1])
def partFile = new File(outFile.getPath() + '.part')
outFile.getParentFile().mkdirs()
partFile.delete()

def imageData = getCurrentImageData()
def hierarchy = imageData.getHierarchy()
String unit = imageData.getServer().getPixelCalibration().getPixelWidthUnit()
def detections = hierarchy.getDetectionObjects()
def annotations = hierarchy.getAnnotationObjects()
def clsOf = { o -> o.getPathClass()?.toString() ?: 'Unclassified' }
def details = [mode: mode, output: outFile.getPath(), detections: detections.size(), annotations: annotations.size()]
boolean ok = false
String message = ''
long t0 = System.currentTimeMillis()

// header of a tab-separated file, as written by QuPath's MeasurementExporter
def headerOf = { File f -> f.withReader('UTF-8') { it.readLine() }?.split('\t', -1) as List }
def countLines = { File f -> long n = 0; f.withReader('UTF-8') { r -> while (r.readLine() != null) n++ }; n }

if (mode == 'cells') {
    // the export's own column list (metadata columns + every measurement), minus the columns that later
    // steps add: signed distances belong to the context export
    def table = new qupath.lib.gui.measure.ObservableMeasurementTableData()
    table.setImageData(imageData, detections)
    def allNames = new ArrayList(table.getAllNames())
    def later = allNames.findAll { it.startsWith('Signed distance') }
    def wanted = allNames.findAll { !it.startsWith('Signed distance') }
    details.excluded_columns = later
    saveDetectionMeasurements(partFile.getPath(), wanted as String[])
    def header = headerOf(partFile)
    details.columns = header.size()
    details.header_sample = header.take(12)
    def leaked = header.findAll { it.startsWith('Signed distance') }
    // rows are counted by v3ctl.py post (DuckDB), which reads the whole file anyway
    ok = partFile.length() > 0 && header.contains('Object ID') && leaked.isEmpty()
    message = ok ? "${header.size()} columns (${later.size()} signed-distance columns left to the context export)" :
            (leaked ? "signed-distance columns in the cells export: ${leaked}" : "no 'Object ID' column in ${header.take(8)}")
} else if (mode == 'context') {
    def distanceCols = qupath.lib.objects.PathObjectTools.getAvailableFeatures(detections)
            .findAll { it.startsWith('Signed distance') }.sort()
    def wanted = ['Object ID', 'Parent', "Centroid X ${unit}".toString(), "Centroid Y ${unit}".toString()] + distanceCols
    details.distance_columns = distanceCols
    if (distanceCols.isEmpty()) {
        message = "no 'Signed distance…' measurements (run step 3 first)"
    } else {
        saveDetectionMeasurements(partFile.getPath(), wanted as String[])
        def header = headerOf(partFile)
        def absent = wanted.findAll { !(it in header) }
        details.columns = header
        ok = absent.isEmpty()
        message = ok ? "${header.size()} columns incl. ${distanceCols.size()} distance columns" :
                "requested columns missing from the export: ${absent}"
    }
} else if (mode == 'annotations') {
    saveAnnotationMeasurements(partFile.getPath())
    long rows = countLines(partFile) - 1
    details.rows = rows
    details.columns = headerOf(partFile)?.size()
    ok = rows == annotations.size()
    message = "${rows} rows for ${annotations.size()} annotations"
} else if (mode == 'islets') {
    def objs = annotations.findAll { !(clsOf(it) in ['Tissue', 'pLN']) }
    details.exported = objs.size()
    details.by_class = new TreeMap(objs.countBy { clsOf(it) })
    exportObjectsToGeoJson(objs, partFile.getPath(), 'FEATURE_COLLECTION')
    // readObjects(File) picks its reader from the extension, and '.part' is not one of them
    def back = partFile.withInputStream { PathIO.readObjectsFromGeoJSON(it) }
    def wantIds = objs.collect { it.getID().toString() } as Set
    def gotIds = back.collect { it.getID().toString() } as Set
    details.reread = back.size()
    ok = back.size() == objs.size() && gotIds == wantIds
    message = "${objs.size()} objects exported, ${back.size()} re-read, IDs ${gotIds == wantIds ? 'match' : 'DIFFER'}"
} else {
    throw new IllegalArgumentException("unknown export mode '${mode}'")
}

if (ok) {
    java.nio.file.Files.move(partFile.toPath(), outFile.toPath(), java.nio.file.StandardCopyOption.REPLACE_EXISTING)
    details.bytes = outFile.length()
}
details.seconds = (System.currentTimeMillis() - t0) / 1000.0
// cells/context TSVs are intermediates (v3ctl.py post deletes them); the rest are deliverables
details[(mode in ['cells', 'context']) ? 'transient_outputs' : 'outputs'] = ok ? [outFile.getPath()] : []
V3_RESULT = [ok: ok, changed: false, message: message.toString(), details: details]
println "[exports] ${mode}: ${ok ? 'OK' : 'FAILED'} -- ${message} -> ${outFile} (${details.seconds} s)"
