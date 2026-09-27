/*
 * 03_cell_detection.groovy  —  SUPERSEDED
 *
 * This stage is no longer a placeholder. Cell detection is implemented in
 *
 *     /home/smith6jt/IO60panc2nd/scripts/0_cell_detection.groovy
 *
 * Do not run detection from the QuPath GUI. Detection being an unscripted manual step is exactly
 * what caused the duplicate-cell corruption: runs that used selectAnnotations() selected the Tissue
 * annotation AND every nested Islet_N, and InstanSeg segments each selected parent independently
 * without deduplicating across them. Every islet nucleus was therefore segmented twice, with two
 * distinct UUIDs, in 7 of 20 donors (155,669 duplicate pairs). Nothing downstream caught it.
 *
 * 0_cell_detection.groovy selects only the Tissue annotation, refuses to run if anything else is
 * selected, clears prior detections, and verifies afterwards that no two detections share a
 * centroid. Use it.
 */

println "[03] SUPERSEDED — run scripts/0_cell_detection.groovy instead."
println "[03] Current detection count: ${getDetectionObjects().size()}"
