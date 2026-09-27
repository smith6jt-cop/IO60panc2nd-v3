/*
 * Simple script to apply the composite object classifier
 *
 */

// ==================== PARAMETERS ====================
// Name of the composite classifier saved in
//   classifiers/object_classifiers/<CLASSIFIER_NAME>.json
String CLASSIFIER_NAME = 'NormPancMacsGal3'

// ====================================================

// Reset detections is necessary
//resetDetectionClassifications();


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
