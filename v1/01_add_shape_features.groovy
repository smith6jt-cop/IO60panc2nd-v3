/*
 * 01_add_shape_features.groovy
 *
 * Workflow step 3. Scripted equivalent of QuPath UI:
 *   Analyze -> Calculate features -> Add shape features...
 *
 * Adds area, length, circularity, solidity, diameter, and nucleus/cell-ratio
 * measurements to all detection objects. Running this as a script lets
 * step 3 be driven by Run -> Run for project instead of clicking through
 * 15 images in the UI.
 *
 * Requires: detection objects from step 00b (cell detection).
 */

evaluate(new File(buildFilePath(PROJECT_BASE_DIR, 'scripts', 'lib', 'common.groovy')))

def n = getDetectionObjects().size()
if (n == 0) {
    println "[01] No detection objects. Run cell detection (step 00b) first."
    return
}

selectDetections()
addShapeMeasurements(
    "AREA",
    "LENGTH",
    "CIRCULARITY",
    "SOLIDITY",
    "MAX_DIAMETER",
    "MIN_DIAMETER",
    "NUCLEUS_CELL_RATIO"
)
detectionToAnnotationDistancesSigned(false)
fireHierarchyUpdate()

println "[01] Added shape measurements to ${n} detections."
