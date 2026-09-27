/*
 * 06_composite_classifier.groovy
 *
 * Workflow step 7. Apply the composite object classifier
 *   "M2Gal3->Iba1CD68 2D"
 * to all detections, using the normalized (stage 05) measurement columns.
 *
 * QuPath composite classifiers chain individual object classifiers in the
 * order saved. The preferred setup is:
 *   - Classifier A: "Iba1CD68_2D_norm" (trained on Cell: Iba1 Mean (norm) and
 *     Cell: CD68 Mean (norm)) splits cells into Iba1+CD68+ macrophages vs
 *     the rest.
 *   - Classifier B: "M2Gal3_norm" (trained on Cell: M2Gal3 Mean (norm))
 *     further splits the macrophage class into M2-like (M2Gal3+) vs other.
 *   - Composite "M2Gal3->Iba1CD68_2D" = A then B, saved via
 *     Classify > Object classification > Create composite classifier.
 *
 * Training is interactive (Classify > Object classification > Train object
 * classifier). Train each component classifier on a representative image
 * whose norm_stats (from stage 05) sit near the middle of the range across
 * all 15 images; see NORMALIZATION.md for details.
 *
 * Precondition: stage 05 has added "(norm)" measurement columns for the
 * object-classifier channels. If CD68 is also used by the composite
 * classifier but is not in OBJECT_CLF_CHANNELS, extend that list in
 * common.groovy before training.
 *
 * This script is intentionally minimal. Run via: Run -> Run for project.
 */

evaluate(new File(buildFilePath(PROJECT_BASE_DIR, 'scripts', 'lib', 'common.groovy')))

// ==================== PARAMETERS ====================
// Name of the composite classifier saved in
//   classifiers/object_classifiers/<CLASSIFIER_NAME>.json
String CLASSIFIER_NAME = 'M2Gal3->Iba1CD68_2D'
// ====================================================

def classifierPath = buildFilePath(PROJECT_BASE_DIR, 'classifiers', 'object_classifiers', "${CLASSIFIER_NAME}.json")
if (!new File(classifierPath).exists()) {
    println "[06] ERROR: composite classifier not found at ${classifierPath}"
    println "[06]        Train component classifiers and create the composite interactively."
    return
}

int nDetBefore = getDetectionObjects().size()
if (nDetBefore == 0) { println "[06] No detections - nothing to classify."; return }

runObjectClassifier(CLASSIFIER_NAME)

// Report class distribution after classification
def counts = [:]
getDetectionObjects().each { det ->
    def c = det.getPathClass()?.toString() ?: '(unclassified)'
    counts[c] = (counts[c] ?: 0) + 1
}
println "[06] Applied '${CLASSIFIER_NAME}' to ${nDetBefore} detections. Class distribution:"
counts.sort { -it.value }.each { cls, n ->
    println String.format("[06]   %-30s %d  (%.1f%%)", cls, n, 100.0 * n / nDetBefore)
}

fireHierarchyUpdate()
