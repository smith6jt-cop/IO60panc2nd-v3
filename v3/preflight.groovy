/*
 * preflight.groovy -- run_v3.sh step 'preflight': is this image safe to process?
 *
 * Runs first in passes I and C. Every check runs and is reported; nothing throws, so the
 * provenance JSON records the full picture even when the image fails. The verdict goes to
 * V3_RESULT (ok / changed / details), which the run_v3 footer turns into pass/fail; the footer
 * saves only when every check passed AND something changed.
 *
 * Checks
 *   1. the channels scripts 0 (InstanSeg's 13 inputs) and 1 (INS/GCG/SST) need, each exactly once
 *      (a duplicated name makes ColorTransforms pick the first match silently);
 *   2. pixel size 0.5078 um in x and y, strictly. Scripts 0 and 1 force 0.5078 with
 *      setPixelSizeMicrons(), but a headless run cannot write that back to project.qpproj, and on the
 *      next load QuPath takes the calibration from the project entry. An entry with a different value
 *      would therefore flip between steps (6450panc's entry says 0.5077663810243286 um). Fix the entry
 *      in the GUI, not here;
 *   3. exactly one Tissue annotation -- so none is nested in another, which script 0 would segment twice;
 *   4. every pLN annotation lies >= 99% (by area) inside the Tissue union: script 0 only detects
 *      inside Tissue, so pLN area outside it would silently have no cells.
 * Change
 *   5. removes the stale 'Other: Gal3Norm' annotations (derived objects of the retired Gal3 density
 *      map; archived by export_annotations.sh) and lists what it removed.
 */

import org.locationtech.jts.operation.union.UnaryUnionOp

// ==================== PARAMETERS ====================
List<String> REQUIRED_CHANNELS = ['DAPI', 'CD68', 'CD44', 'Vimentin', 'CD99', 'E-cadherin', 'M2Gal3', 'HLA-A',
                                  'Pan-Cytokeratin', 'CD79a', 'CD3e', 'Beta-actin', 'Iba1', 'INS', 'GCG', 'SST']
double EXPECTED_PIXEL_UM = 0.5078
double PIXEL_TOLERANCE   = 1e-6      // strict on purpose, see check 2
double PLN_MIN_INSIDE    = 0.99
String STALE_CLASS       = 'Other: Gal3Norm'
// ====================================================

def imageData = getCurrentImageData()
def server = imageData.getServer()
def hierarchy = imageData.getHierarchy()
def clsOf = { o -> o.getPathClass()?.toString() ?: 'Unclassified' }
def checks = [:]

// 1. channels
def names = server.getMetadata().getChannels().collect { it.getName() }
def missing = REQUIRED_CHANNELS.findAll { !(it in names) }
def duplicated = REQUIRED_CHANNELS.findAll { req -> names.count { it == req } > 1 }
checks.channels = [ok: missing.isEmpty() && duplicated.isEmpty(), n_channels: names.size(),
                   missing: missing, duplicated: duplicated]

// 2. pixel size
def cal = server.getPixelCalibration()
double pw = cal.getPixelWidth().doubleValue(), ph = cal.getPixelHeight().doubleValue()
boolean unitOk = cal.getPixelWidthUnit() == 'µm' && cal.getPixelHeightUnit() == 'µm'
checks.pixel_size = [ok   : unitOk && Math.abs(pw - EXPECTED_PIXEL_UM) <= PIXEL_TOLERANCE &&
                            Math.abs(ph - EXPECTED_PIXEL_UM) <= PIXEL_TOLERANCE,
                     width: pw, height: ph, unit: cal.getPixelWidthUnit(), expected: EXPECTED_PIXEL_UM]

// 3. exactly one Tissue
def anns = hierarchy.getAnnotationObjects()
def tissues = anns.findAll { clsOf(it) == 'Tissue' }
def tissueSet = Collections.newSetFromMap(new IdentityHashMap())
tissueSet.addAll(tissues)
def hasTissueAncestor = { o ->
    def p = o.getParent()
    while (p != null) { if (tissueSet.contains(p)) return true; p = p.getParent() }
    return false
}
double umPerPx2 = pw * ph
checks.tissue = [ok     : tissues.size() == 1,
                 count  : tissues.size(),
                 nested : tissues.count { hasTissueAncestor(it) },
                 details: tissues.collect { [id: it.getID().toString(), name: it.getName(),
                                             parent: it.getParent()?.isRootObject() ? 'ROOT' : clsOf(it.getParent()),
                                             area_um2: it.getROI().getArea() * umPerPx2] }]

// 4. pLN inside the Tissue union
def plns = anns.findAll { clsOf(it) == 'pLN' }
def plnRows = []
boolean plnOk = true
if (!plns.isEmpty()) {
    def tissueUnion = tissues.isEmpty() ? null :
            UnaryUnionOp.union(tissues.collect { it.getROI().getGeometry().buffer(0) })
    plns.each { p ->
        def row = [id: p.getID().toString(), name: p.getName(), area_um2: p.getROI().getArea() * umPerPx2]
        try {
            def g = p.getROI().getGeometry().buffer(0)
            double frac = (tissueUnion == null || g.getArea() <= 0) ? 0.0 : g.intersection(tissueUnion).getArea() / g.getArea()
            row.fraction_inside_tissue = frac
            row.ok = frac >= PLN_MIN_INSIDE
        } catch (Throwable t) {
            row.error = t.toString()
            row.ok = false
        }
        plnOk &= row.ok
        plnRows << row
    }
}
checks.pln_inside_tissue = [ok: plnOk, count: plns.size(), min_fraction: PLN_MIN_INSIDE, details: plnRows]

// 5. stale Gal3Norm annotations
def stale = anns.findAll { clsOf(it) == STALE_CLASS }
def removed = stale.collect { [id: it.getID().toString(), name: it.getName(),
                               area_um2: it.getROI()?.getArea() * umPerPx2, n_children: it.nChildObjects()] }
if (!stale.isEmpty()) {
    hierarchy.removeObjects(stale, true)   // keepChildren: any detections under them move up to their parent
    fireHierarchyUpdate()
}

boolean ok = checks.values().every { it.ok }
def failed = checks.findAll { k, v -> !v.ok }.keySet()
V3_RESULT = [ok     : ok,
             changed: !stale.isEmpty(),
             message: ok ? 'all checks passed' : "failed: ${failed.join(', ')}".toString(),
             details: [checks             : checks,
                       image_type         : imageData.getImageType()?.toString(),
                       removed_class      : STALE_CLASS,
                       removed_count      : removed.size(),
                       removed            : removed]]
println "[preflight] ${ok ? 'PASS' : 'FAIL'}: " + checks.collect { k, v -> "${k}=${v.ok}" }.join(', ') +
        "; removed ${removed.size()} '${STALE_CLASS}'" + (ok ? '' : " -- ${V3_RESULT.message}")
checks.findAll { k, v -> !v.ok }.each { k, v -> println "[preflight]   ${k}: ${v}" }
