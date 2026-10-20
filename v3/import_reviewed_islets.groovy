/*
 * import_reviewed_islets.groovy -- run_v3.sh step 'R_import' (pass R, used by hpg/j_cpu.sbatch).
 *
 * Pass J on HiPerGator works on the .qpdata that script 2 produced there, but the islets were
 * reviewed locally (gate 1). The reviewed set travels as GeoJSON (step 'R_export', pass X, writes
 * results/v3/reviewed_islets/<tag>.geojson: every annotation except Tissue and pLN). This step makes
 * that set authoritative: it removes every annotation except Tissue and pLN (detections under them
 * move up; script 8 re-parents them) and imports the file, keeping object IDs, names, classes and
 * measurements such as 'Islet: UID'. Verified by counting and matching IDs; saved by the footer.
 *
 *   args[0] = reviewed GeoJSON
 */

import qupath.lib.io.PathIO

if (args == null || args.length < 1) throw new IllegalArgumentException('usage: import_reviewed_islets.groovy <reviewed.geojson>')
def src = new File(args[0])
def hierarchy = getCurrentHierarchy()
def clsOf = { o -> o.getPathClass()?.toString() ?: 'Unclassified' }
def keep = ['Tissue', 'pLN'] as Set

if (!src.isFile()) {
    V3_RESULT = [ok: false, changed: false, message: "reviewed GeoJSON not found: ${src}".toString(), details: [:]]
    println "[import_reviewed_islets] ERROR: ${V3_RESULT.message}"
} else {
    def incoming = PathIO.readObjects(src)
    def badClass = incoming.findAll { clsOf(it) in keep }
    def old = hierarchy.getAnnotationObjects().findAll { !(clsOf(it) in keep) }
    // the file's hash lets j_cpu.sbatch re-import when a newer review export arrives
    def details = [source: src.getPath(), source_sha256: java.security.MessageDigest.getInstance('SHA-256').digest(src.bytes).encodeHex().toString(),
                   removed: old.size(), removed_by_class: new TreeMap(old.countBy { clsOf(it) }),
                   incoming: incoming.size(), incoming_by_class: new TreeMap(incoming.countBy { clsOf(it) })]
    if (!badClass.isEmpty()) {
        V3_RESULT = [ok: false, changed: false, details: details,
                     message: "the file contains ${badClass.size()} Tissue/pLN objects; it must hold only reviewed islet-type objects".toString()]
    } else {
        hierarchy.removeObjects(old, true)          // keepChildren: detections move up, script 8 re-parents them
        // added at the root, NOT insertPathObjects: placing each object by containment tests it against every
        // detection under its parent -- 7,000 objects x 1 M cells did not finish in 8 min (measured 2026-10-01).
        // Nothing here needs the nesting: script 8 parents the cells, and steps 3 and 4 work on geometry.
        hierarchy.addObjects(incoming)
        fireHierarchyUpdate()
        def now = hierarchy.getAnnotationObjects().findAll { !(clsOf(it) in keep) }
        boolean idsMatch = (now.collect { it.getID().toString() } as Set) == (incoming.collect { it.getID().toString() } as Set)
        details.after = now.size()
        details.ids_match = idsMatch
        boolean ok = now.size() == incoming.size() && idsMatch
        V3_RESULT = [ok: ok, changed: true, details: details,
                     message: "replaced ${old.size()} annotations with ${incoming.size()} reviewed ones (IDs ${idsMatch ? 'match' : 'DIFFER'})".toString()]
    }
    println "[import_reviewed_islets] ${V3_RESULT.ok ? 'OK' : 'FAILED'}: ${V3_RESULT.message}"
}
