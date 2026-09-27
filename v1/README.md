# QuPath 59-Channel Pancreas Pipeline

Pipeline for detecting, measuring, and classifying cells in 59-channel multiplex immunofluorescence whole-slide pancreas images in QuPath v0.7. Designed for a batch of 15 images, restart-safe at the tile level, and with per-image channel normalization so a single trained object classifier generalizes across all images.

## Workflow (as specified)

| # | Step | How it runs | Script |
|---|------|-------------|--------|
| 1 | Check/create Islets pixel classifier | QuPath UI | — |
| 2 | Panc_Objects: channel norm + Tissue + Islet + erode/dilate | script, Run for project | `00_panc_objects.groovy` |
| 2.5 | **Cell detection** (not listed but required before step 3) | StarDist/Cellpose, UI or script | `00b_cell_detection.groovy` (placeholder) |
| 3 | Add shape features | script (wraps Analyze → Calculate features → Add shape features...) | `01_add_shape_features.groovy` |
| 4 | Signed Distance to Islet annotations | script (wraps Analyze → Spatial analysis → Signed Distance to Annotations 2D) | `02_signed_distance.groovy` |
| 5 | Islet Expander: generate IsletExpanded bands | script, Run for project | `03_islet_expander.groovy` |
| 6 | FastCellMeasurements: shape + intensity measurements | script, Run for project | `04_fast_cell_measurements.groovy` |
| 6.5 | **Channel normalization** (inserted; see NORMALIZATION.md) | script, Run for project | `05_normalize_channels.groovy` |
| 7 | Composite classifier `M2Gal3->Iba1CD68 2D` | trained interactively, then applied by script | `06_composite_classifier.groovy` |
| 8 | Check thresholders for B3TUBB, CD34, Podoplanin | QuPath UI inspection | — |
| 9 | `01b run one classifier tiled` for Nerve, CD34, Lymphatic | script, Run for project | `07_tiled_pixel_classifiers.groovy` |
| 10 | FastHierarchyAnnotationsDetections: parent detections | script, Run for project | `08_hierarchy_parenting.groovy` |

Two deliberate insertions: **2.5** (cell detection) is required by steps 3, 4, and 6 but was not in the original list; and **6.5** (channel normalization) is the reason for this restructuring and is detailed in `NORMALIZATION.md`.

## Ordering constraints

1. **Cell detection must precede step 3.** Shape features and signed distance both operate on detections. If your cell-detection step is already part of the UI workflow between Panc_Objects and the Analyze steps, `00b_cell_detection.groovy` can be skipped; fill it in only if you want the whole pipeline driven by Run-for-project.

2. **Signed distance must precede Islet Expander (step 5).** `02_signed_distance.groovy` passes only `CLASS_ISLET` to `DistanceTools.detectionToAnnotationDistancesSigned`, so the measurement is specifically "distance to Islet core" regardless of when it runs. However, running it before IsletExpanded exists matches the specified workflow and guarantees there is no alternative class to pick accidentally in the UI version.

3. **Stage 05 (normalization) must run between FastCellMeasurements (step 6) and the composite classifier (step 7).** The classifier reads `(norm)` columns that do not exist until stage 05 writes them. If you add channels to `OBJECT_CLF_CHANNELS` (e.g. `CD68` to match the Iba1CD68 2D classifier feature set), stage 05 must be re-run before step 7.

4. **Step 8 (thresholder checks) is pure UI inspection.** No script is provided; verify the three thresholders load and look correct in `Classify → Pixel classification → Load classifier` for each of `B3TUBB`, `CD34`, and `Podoplanin`.

## Directory layout

```
<project root>/
  scripts/
    lib/
      common.groovy
    00_panc_objects.groovy
    00b_cell_detection.groovy
    01_add_shape_features.groovy
    02_signed_distance.groovy
    03_islet_expander.groovy
    04_fast_cell_measurements.groovy
    05_normalize_channels.groovy
    06_composite_classifier.groovy
    07_tiled_pixel_classifiers.groovy
    08_hierarchy_parenting.groovy
  classifiers/
    pixel_classifiers/
      Islets.json
      Tissue.json           (if used; see 00_panc_objects.groovy)
      CD34dil1430.json
      B3Tub_Ero0.5_Thr600.json
      PDPNmorphClos0-5_Thr300.json
    object_classifiers/
      Iba1CD68_2D_norm.json           (component, trained interactively)
      M2Gal3_norm.json                (component, trained interactively)
      M2Gal3->Iba1CD68_2D.json        (composite, created via UI)
  results/
    Blood_vessel__<image>.geojson        (GeoJSON only; CD34 output class label)
    Nerve__<image>.geojson               (GeoJSON only)
    Lymphatic__<image>.geojson           (GeoJSON only)
    *.scratch.jsonl                      (tile-level checkpoints)
    norm_stats.csv                       (per-image normalization anchors)
```

## Summary of changes from the original five scripts

**`common.groovy`** — added `CLASS_BLOOD_VESSEL = 'Blood_vessel'` and `OBJECT_CLF_CHANNELS = ['M2Gal3','CD34','Iba1']`. Added `getMeas` reader companion to the existing `putMeas` writer. Load path in the header comment now matches the loader used in the other scripts (`scripts/lib/`).

**00 (was `Panc_Objects.groovy`)** — expanded the channel alias map so `M2Gal3`, `CD34`, and `Iba1` are canonicalized before any classifier reads them. Post-alias verification fails fast if any of the three is absent. Stage is idempotent: existing Islet annotations are removed before recreation. Header now points at `00b_cell_detection.groovy` for the dependency.

**00b (new)** — cell-detection placeholder. Not part of the original five scripts; inserted because steps 3, 4, 6 need detection objects that Panc_Objects does not create.

**01 (new)** — thin script wrapper for `Analyze → Calculate features → Add shape features`. Calls `addShapeMeasurements(...)` with the same seven features the original `FastCellMeasurements` set.

**02 (new)** — thin script wrapper for `Analyze → Spatial analysis → Signed Distance to Annotations 2D`. Calls `DistanceTools.detectionToAnnotationDistancesSigned(imageData, [CLASS_ISLET])`.

**03 (was `Islet_Expander.groovy`)** — functionally unchanged; header and log prefix updated.

**04 (was `FastCellMeasurements.groovy`)** — adds an optional `MEASURE_CHANNELS` parameter that, when non-null, restricts intensity measurement to a channel subset via `TransformedServerBuilder.extractChannels`. Leave `null` for full 59-channel coverage. Shape-measurement call is retained (harmless if 01 already ran).

**05 (new)** — per-image channel normalization. Uses full Tissue as the reference mask (no Islet or Blood_vessel exclusion). See `NORMALIZATION.md`.

**06 (new)** — thin wrapper for `runObjectClassifier('M2Gal3->Iba1CD68_2D')`. Training of the two component classifiers and creation of the composite are interactive UI steps; the script only runs the already-saved composite.

**07 (was `01b_run_one_classifier_tiled.groovy`)** — runs `B3Tubb → Nerve`, `CD34 → Blood_vessel`, and `Podoplanin → Lymphatic` in a single invocation per image. Caveolin is removed per the current workflow. All outputs are written to per-marker GeoJSON under `results/` and are removed from the hierarchy after tile serialization because the CD34 coverage is too large to retain loaded. The throwaway-class safety model and preserved-class snapshot/verification are retained. Resume and `TEST_MODE` behavior unchanged.

**08 (was `FastHierarchyAnnotationsDetections.groovy`)** — references `CLASS_BLOOD_VESSEL` from `common.groovy` (fixes the previous `'Blood_Vessel'` vs `'Blood_vessel'` case mismatch). Under the current pipeline the Blood_vessel parenting bucket is empty because stage 07 does not write those annotations into the hierarchy; the priority slot is retained for the case of a manual, spatially-subsetted re-import. Optional area-fraction parenting rule included as commented code.

## Batch execution

All scripts respect QuPath's `Run → Run for project`. Recommended sequence:

1. Step 1: verify the `Islets` pixel classifier exists in `classifiers/pixel_classifiers/`.
2. Run 00 for project. Inspect two images, confirm Islet counts and channel names (all three classifier channels should appear in the post-alias log).
3. Run 00b for project — or run cell detection interactively for each image if you prefer.
4. Run 01 for project.
5. Run 02 for project. Verify each detection has a `Signed distance to annotation with Islet` measurement before proceeding.
6. Run 03 for project. Confirm `IsletExpanded` annotations appear in each image.
7. Run 04 for project. This is the long-running stage on a 59-channel image; plan compute time accordingly.
8. Run 05 for project. Review `results/norm_stats.csv` — `p_low` should be stable within ~2× across images per channel; `p_high` may vary more (see `NORMALIZATION.md`).
9. Train the two component object classifiers interactively on a representative image whose `norm_stats` sit near the middle of the range. Save `Iba1CD68_2D_norm` and `M2Gal3_norm`. Create the composite `M2Gal3->Iba1CD68_2D` via `Classify → Object classification → Create composite classifier`.
10. Run 06 for project.
11. Step 8: open each thresholder (`B3TUBB`, `CD34`, `Podoplanin`) in the UI on one test image and confirm it produces sensible output.
12. Run 07 for project in `TEST_MODE = true` on one image first; inspect one tile of each marker. Then set `TEST_MODE = false` and run for project.
13. Run 08 for project.

## Assumptions and open questions

- **Cell detection tool is not specified.** The placeholder `00b_cell_detection.groovy` lists StarDist and Cellpose as the likely options but does not commit to either. Confirm which tool, channel, and parameters you use, and I can fill in the script.
- **CD68 is not in `OBJECT_CLF_CHANNELS`.** The composite classifier name `M2Gal3->Iba1CD68 2D` implies CD68 is used. If CD68 measurements are used as a classifier feature, add `'CD68'` to `OBJECT_CLF_CHANNELS` in `common.groovy` and re-run stage 05 so `(norm)` columns exist for it. The current list (`M2Gal3, CD34, Iba1`) is what you stated in the original request; adding CD68 is a one-line change.
- QuPath v0.7 is assumed for `ObjectMeasurements.addIntensityMeasurements`, `TransformedServerBuilder.extractChannels`, and `DistanceTools.detectionToAnnotationDistancesSigned`. All three are present in v0.6 as well.
- The pixel classifiers (`CD34dil1430`, `B3Tub_Ero0.5_Thr600`, `PDPNmorphClos0-5_Thr300`) and the `Islets` classifier must already exist in `classifiers/pixel_classifiers/`. They are not regenerated.
