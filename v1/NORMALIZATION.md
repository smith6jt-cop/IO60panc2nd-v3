# Channel Normalization for the Object Classifier

## Purpose

The object classifier in QuPath reads per-detection measurements (`Cell: CD34: Mean`, etc.) and does not see raw pixel values. Display and visualization transforms do not affect those measurement values. If the same object classifier is to be applied across all 15 images, the measurement values it reads must be on a comparable scale across images — which raw fluorescence values generally are not, due to exposure drift, staining-batch variation, and slide-to-slide autofluorescence differences.

Stage 05 (`05_normalize_channels.groovy`) addresses this by computing image-specific anchors from tissue pixels and adding a normalized copy of each per-detection intensity measurement alongside the raw one. The classifier is trained and applied on the normalized copies.

## Method

For each image, for each channel in `OBJECT_CLF_CHANNELS = ['M2Gal3', 'CD34', 'Iba1']`:

1. **Build the reference mask.** Let `R = Tissue` — the union of all `Tissue` annotations at full extent, with no exclusions. Islet and Blood_vessel regions are **not** subtracted. This is deliberate in this pipeline: Blood_vessel is not present in the hierarchy (stage 07 writes CD34 output only to GeoJSON because the coverage is too large to retain), so a mask definition that depends on it would be inconsistent across runs. Including Islet also keeps the reference stable when the Islet segmentation changes.

2. **Sample pixels.** Rasterize `R` to a boolean mask at `SAMPLE_DOWNSAMPLE = 16` (default). Read the channel image at the same downsample via a single-channel `TransformedServerBuilder` view. Keep pixel values where the mask is true.

3. **Compute anchors.**
   - `p_low  = quantile(sampled pixels, 0.50)` — median of tissue intensity.
   - `p_high = quantile(sampled pixels, 0.99)` — 99th percentile of tissue intensity.

4. **Normalize each detection measurement.** For each `compartment ∈ {Cell, Cytoplasm, Membrane, Nucleus}` and each `stat ∈ STATS_TO_NORMALIZE` (default: `Mean`), write
   ```
   norm = (raw - p_low) / max(p_high - p_low, 1e-6)
   ```
   to a new measurement column named `"<compartment>: <channel>: <stat> (norm)"`.

5. **Persist anchors.** Append `(image, channel, p_low, p_high, span, n_pixels, p_low_quantile, p_high_quantile)` to `results/norm_stats.csv`.

## Interpretation of normalized values

With the full-Tissue reference, `p_high` sits near the brightest tissue pixels on the slide, which for a positively-staining marker means inside the stained structures (vessels for CD34, macrophages for M2Gal3 and Iba1). That shifts the meaning of the normalized scale slightly compared to a background-only reference:

- `norm ≈ 0` — the detection sits at the image's tissue median (roughly background-to-low-signal tissue).
- `norm ≈ 1` — at the 99th percentile of tissue intensity. Most positive cells land below 1; only the very brightest exceed it.
- `norm > 1` — within the brightest ~1% of tissue pixels.
- `norm < 0` — below the tissue median.

The classifier is trained on these values and learns whatever threshold separates the positive from the negative class in this scale. Because the same rescaling is applied image by image, a trained threshold transfers across images even though the raw intensity at that threshold may differ.

## Why this method

**Per-image, pixel-level anchors.** Computing anchors from pixels inside a per-image mask is the core property that makes the normalization work across images. Each image gets its own `p_low` and `p_high`, so image-to-image exposure, laser-power, and staining-batch differences are absorbed before the classifier ever sees a value.

**Robust statistics.** Median and 99th percentile are not shifted by a small number of extremely bright pixels (saturation, imaging artifacts), unlike mean and standard deviation. This matters for 16-bit fluorescence data where outliers are common.

**No cross-stage dependency.** The reference mask depends only on `Tissue`, which is produced in stage 00. Stage 05 does not require stage 07 output, stage 08 parenting, or any information about Islet/Blood_vessel to compute anchors. This is a simplification over an exclusion-based mask and avoids situations where partial pipeline runs would produce inconsistent normalization.

**No change to display or raw measurements.** Channel display settings, pixel classifier JSONs, and raw measurement columns are untouched. The `(norm)` columns sit alongside the raw ones. At any point you can compare raw vs normalized values to diagnose a regression.

## Caveat: prevalence sensitivity

Because the mask is full Tissue, `p_high` is shifted toward positive staining wherever the positive population is large. An image where, say, 10% of pixels fall inside vessels will have a higher `p_high` for CD34 than an image where 2% do, purely because of prevalence.

For a classifier trained on a single representative image, this means:

- If the training image has **near-average prevalence**, the classifier generalizes well.
- If the training image has **unusually high or low** prevalence of a marker, the classifier's decision threshold may shift the opposite direction on images with different prevalence.

Two practical mitigations:

1. **Select the training image carefully.** Pick one whose `p_high` values in `results/norm_stats.csv` sit in the middle of the range across all 15 images — not an outlier on either end. This is the primary mitigation.

2. **Use multiple training images.** Train the object classifier on detections pooled from two or three images that span the `p_high` range in `norm_stats.csv`. QuPath supports this: annotate representative positive and negative cells in each of those images, then train with all of them as inputs.

If prevalence variation turns out to be too large for the classifier to handle with either mitigation, the fallback is to tighten `P_HIGH` to `0.995` so the anchor sits further into the tail and is less affected by the middle of the positive population.

## Alternatives considered

| Method | Fit for this pipeline |
|--------|------------------------|
| **Full-Tissue reference, robust anchors** (current) | Matches user requirements: no dependency on Blood_vessel in hierarchy, no exclusions. Caveat above applies. |
| Tissue minus (Islet ∪ Blood_vessel) | Requires Blood_vessel in hierarchy; rejected for resource reasons. Also adds cross-stage dependency. |
| Per-image z-score `(raw - mean) / std` | `mean` and `std` are sensitive to positive-cell prevalence and saturated pixels. Worse than the current method on exactly the skewed distributions encountered here. |
| Cross-image quantile matching to a reference image | Depends on the chosen reference and distorts the upper tail. Adds a cross-image dependency that complicates re-runs when images are added or removed. |
| Divide by a reference channel (e.g. DAPI) | Only valid if the reference channel is uniform across tissue. Not true in this 59-plex panel; DAPI varies with nuclear density. |
| Fit a GMM to background vs signal per channel | More principled but requires per-channel model fitting and assumes a bimodal structure that is not always present. The `p_low`/`p_high` rule captures the same idea with no model assumptions. |

## Diagnostic checks

After running stage 05 on all 15 images, inspect `results/norm_stats.csv`:

- **`p_low` should be stable within a factor of ~2** across images for a given channel. Large variation suggests autofluorescence differences or a tissue-mask problem.
- **`p_high` may vary more widely** than `p_low` because of prevalence sensitivity. A factor of 3–5× across images is expected; a factor of 20× suggests either a staining protocol change or a dramatic difference in positive-cell prevalence on one slide.
- **`span = p_high - p_low` should be positive and non-trivial.** `span < 10` on a 16-bit channel is a sign that the channel is saturated or empty; normalization will still run (clamped at 1e-6) but the `(norm)` values will be dominated by noise.
- **`n_pixels` should be millions.** A low count suggests the reference mask is tiny and the percentile estimates are noisy; consider a lower `SAMPLE_DOWNSAMPLE`.

A simple sanity check after stage 05 on one image: sort detections by `Cell: CD34: Mean (norm)` descending. The top ~1% should visibly overlap with vessel regions in the image; below `norm = 0` should be detections in clearly non-vessel areas. Same pattern for `M2Gal3` (macrophage-rich areas) and `Iba1` (microglia/macrophages).

## Retraining

If the classifier later underperforms on a new batch of images:

1. Re-run stage 05 on the new images — anchors are computed per image, so the normalization does not need to be redefined.
2. Inspect `norm_stats.csv` for the new images against the training images. If `p_high` has shifted by more than a factor of 5 consistently, the staining protocol may have changed and retraining is warranted.
3. If retraining is needed, follow the two-mitigation guidance in **Caveat: prevalence sensitivity** and train with 2–3 images spanning the range.

## Parameters in `05_normalize_channels.groovy`

| Parameter | Default | Effect of changing |
|-----------|---------|--------------------|
| `SAMPLE_DOWNSAMPLE` | `16.0` | Lower → more pixels sampled → more precise percentiles, more memory. `16` gives millions of in-mask pixels on a typical WSI, which is ample for a 99th percentile estimate. |
| `P_LOW` | `0.50` | Lower anchor quantile. `0.50` is the median. Values in `[0.25, 0.75]` are defensible; higher values shift the zero of the normalized scale toward positive staining. |
| `P_HIGH` | `0.99` | Upper anchor quantile. `0.99` is standard. `0.995` is more conservative against prevalence-driven shifts; `0.999` is less reliable unless `SAMPLE_DOWNSAMPLE` is lowered. |
| `COMPARTMENTS` | all four | Restrict to `['Cell']` only if the classifier uses only whole-cell intensities; reduces the number of measurements written. |
| `STATS_TO_NORMALIZE` | `['Mean']` | Add `'Max'` if the classifier also uses max intensity. `Std_Dev` is not meaningful to rescale with these anchors. |
