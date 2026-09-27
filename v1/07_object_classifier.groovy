/*
 * 07_object_classifier.groovy
 *
 * Stage 07: apply the trained object classifier to all detections.
 *
 * Precondition: the classifier must have been trained against the normalized
 * measurement columns "<comp>: <chan>: Mean (norm)" produced by stage 06.
 * Training is done interactively in QuPath: select detections, open Classify
 * > Object classification > Train object classifier, choose the "(norm)"
 * features for M2Gal3, CD34, and Iba1, save the classifier by name.
 *
 * This script is intentionally minimal. Run via: Run -> Run for project.
 */

evaluate(new File(buildFilePath(PROJECT_BASE_DIR, 'scripts', 'lib', 'common.groovy')))

// ==================== PARAMETERS ====================
String CLASSIFIER_NAME = 'M2Gal3_CD34_Iba1_norm'
// ====================================================

def classifierPath = buildFilePath(PROJECT_BASE_DIR, 'classifiers', 'object_classifiers', "${CLASSIFIER_NAME}.json")
if (!new File(classifierPath).exists()) {
    println "[07] ERROR: object classifier not found at ${classifierPath}"
    println "[07]        Train it interactively first, saving as '${CLASSIFIER_NAME}'."
    return
}

int nDetBefore = getDetectionObjects().size()
if (nDetBefore == 0) { println "[07] No detections - nothing to classify."; return }

runObjectClassifier(CLASSIFIER_NAME)

// Report class distribution after classification
def counts = [:]
getDetectionObjects().each { det ->
    def c = det.getPathClass()?.toString() ?: '(unclassified)'
    counts[c] = (counts[c] ?: 0) + 1
}
println "[07] Applied '${CLASSIFIER_NAME}' to ${nDetBefore} detections. Class distribution:"
counts.sort { -it.value }.each { cls, n ->
    println String.format("[07]   %-30s %d  (%.1f%%)", cls, n, 100.0 * n / nDetBefore)
}

fireHierarchyUpdate()
