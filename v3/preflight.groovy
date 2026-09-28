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
 *   4. every pLN annotation is classified by exactly v3 script 0's rule, so anything step 0 would refuse
 *      fails here, before any GPU time:
 *        disjoint  shares <= STEP0_MAX_OVERLAP_UM2 with the Tissue union: segmented as a region of its
 *                  own (6533, HDL112, HDL115 and 6551 draw their lymph node this way);
 *        contained >= 99% of its area inside Tissue: segmented together with the Tissue;
 *        partial   anything in between: fails (step 0 could neither segment it on its own -- the shared
 *                  area would be segmented twice -- nor skip it -- the rest would get no cells).
 *      Two disjoint pLNs sharing more than STEP0_MAX_OVERLAP_UM2 also fail (step 0 refuses them).
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
double PLN_CONTAINED_MIN = 0.99      // >= 99% inside Tissue: contained
double STEP0_MAX_OVERLAP_UM2 = 10.0  // = v3 script 0's MAX_REGION_OVERLAP_UM2: disjoint means sharing at most this
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

// 4. each pLN: disjoint from Tissue or contained in it; a partial overlap fails
def plns = anns.findAll { clsOf(it) == 'pLN' }
def plnRows = []
def warnings = []
boolean plnOk = true
def geomOf = { o -> def g = o.getROI().getGeometry(); g.isValid() ? g : g.buffer(0) }
if (!plns.isEmpty()) {
    def tissueUnion = tissues.isEmpty() ? null : UnaryUnionOp.union(tissues.collect { geomOf(it) })
    plns.each { p ->
        def row = [id: p.getID().toString(), name: p.getName(), area_um2: p.getROI().getArea() * umPerPx2]
        try {
            if (!p.getROI().isArea()) throw new IllegalArgumentException('not an area ROI')
            def g = geomOf(p)
            double inside = (tissueUnion == null || g.getArea() <= 0) ? 0.0 : g.intersection(tissueUnion).getArea()
            double frac = g.getArea() > 0 ? inside / g.getArea() : 0.0
            row.fraction_inside_tissue = frac
            row.inside_tissue_um2 = inside * umPerPx2
            // the same rule as v3 script 0 (section 2 there), so what step 0 would refuse fails here
            row.kind = row.inside_tissue_um2 <= STEP0_MAX_OVERLAP_UM2 ? 'disjoint' :
                       (frac >= PLN_CONTAINED_MIN ? 'contained' : 'partial')
            row.ok = row.kind != 'partial'
            if (row.kind == 'partial')
                warnings << "pLN '${p.getName() ?: row.id}' shares ${String.format('%.1f', row.inside_tissue_um2)} um^2 " +
                        "(${String.format('%.1f', frac * 100)}%) with Tissue: trim it off the Tissue (<= ${STEP0_MAX_OVERLAP_UM2} um^2) " +
                        "or draw it inside (>= ${PLN_CONTAINED_MIN * 100}%)".toString()
            if (row.kind == 'contained' && frac < 1.0)
                warnings << "pLN '${p.getName() ?: row.id}' is contained in Tissue; its ${String.format('%.1f', (1 - frac) * row.area_um2)} um^2 " +
                        "outside Tissue get no cells".toString()
        } catch (Throwable t) {
            row.error = t.toString()
            row.kind = 'error'
            row.ok = false
        }
        plnOk &= row.ok
        plnRows << row
    }
    // two disjoint pLNs are both segmented on their own by step 0, so an overlap between them fails
    for (int i = 0; i < plns.size(); i++) {
        for (int j = i + 1; j < plns.size(); j++) {
            if (plnRows[i].kind != 'disjoint' || plnRows[j].kind != 'disjoint') continue
            try {
                double shared = geomOf(plns[i]).intersection(geomOf(plns[j])).getArea() * umPerPx2
                if (shared > STEP0_MAX_OVERLAP_UM2) {
                    plnOk = false
                    warnings << "pLN '${plns[i].getName() ?: i}' and '${plns[j].getName() ?: j}' overlap by ${String.format('%.1f', shared)} um^2; " +
                            "step 0 refuses overlapping regions".toString()
                }
            } catch (Throwable t) {
                plnOk = false
                warnings << "pLN overlap check failed: ${t}".toString()
            }
        }
    }
}
checks.pln_vs_tissue = [ok: plnOk, count: plns.size(), disjoint_max_overlap_um2: STEP0_MAX_OVERLAP_UM2,
                        contained_min_fraction: PLN_CONTAINED_MIN, by_kind: plnRows.countBy { it.kind },
                        details: plnRows]

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
                       warnings           : warnings,
                       image_type         : imageData.getImageType()?.toString(),
                       removed_class      : STALE_CLASS,
                       removed_count      : removed.size(),
                       removed            : removed]]
println "[preflight] ${ok ? 'PASS' : 'FAIL'}: " + checks.collect { k, v -> "${k}=${v.ok}" }.join(', ') +
        "; removed ${removed.size()} '${STALE_CLASS}'" + (ok ? '' : " -- ${V3_RESULT.message}")
checks.findAll { k, v -> !v.ok }.each { k, v -> println "[preflight]   ${k}: ${v}" }
plnRows.each { println "[preflight]   pLN '${it.name ?: it.id}': ${it.kind}" + (it.fraction_inside_tissue != null ?
        String.format(' (%.4f of its area inside Tissue)', it.fraction_inside_tissue) : '') + (it.error ? " ${it.error}" : '') }
warnings.each { println "[preflight]   WARNING: ${it}" }
