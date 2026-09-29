# Batch handling (v3)

**Status:** simplified 2026-09-28 per the user. It replaces the earlier gated draft.

## Approach (user ruling, 2026-09-28)

- **Correct batch at the expression level.** UniFORM per-donor alignment is the normalization of record. Where embeddings or clustering are used, an off-the-shelf integration method such as Harmony handles batch.
- **Escalate only on clear evidence.** A pixel- or segmentation-level correction is added only if batch clearly changes cell-type positivity or negativity. Minor shifts in pixel-intensity distributions are not enough.
- **Segmentation just runs.** Every image is re-segmented with InstanSeg on the 13 v3 input channels, with no batch gate in front of it.
- **No extra slides.** There will be no additional PhenoCycler runs to bridge batches.

## Batch structure (facts)

The batch unit is the staining and imaging run; stain groups come from `donor_acquisition_runs.csv`. Runs nest within acquisition groups, each a panel version plus a Fusion software version. The scanner's periodic reference captures are maintenance, not batches.

| Group | Fusion | Stain runs | Panel | ND | Aab+ | T1D |
|---|---|---|---|---|---|---|
| E1 | 2.3.1 | S0–S2 | batch 1 (59 ch) | 6476, HDL115 | 6521 | 6579, 6380, 6551 |
| E2 | 2.4.0 | S3–S5 | batch 1 (59 ch) | 6482, 6374, 6493, 6539, 6442, 6461 | — | 6593, 6436, (6457, excluded) |
| E3 | 2.4.1 | S6–S8 | batch 2 (58/59 ch; no CD38 or b-Catenin1) | none | 6505, 6450 | 6623, 6534, 6566, 6591, 6523 |
| E4 | 2.4.1 | pLNpanc run 1 (HDL112, 6533); run 2 (6347, 6580) | metabolic (70 ch) | HDL112 | 6347 | 6533, 6580 |

- **Disease and batch overlap, but not completely.** T1D is present in every group, and E2 is mostly ND but has two T1D donors. So disease and batch effects can be compared within groups. E3's lack of ND is handled in the statistics.

## What each stage does

| Stage | Batch handling |
|---|---|
| Islet detection | Thresholds are set per image, relative to that image's own background. That is part of the islet method, so detection works across donors and disease. It is checked by the benchmark (representative positives and negatives across donors and disease) and by the user's review. |
| Segmentation (InstanSeg) | None. InstanSeg already normalizes each channel per region. |
| Measurements (script 2) | None: the full raw measurements are kept. |
| Expression and typing | REDSEA, then UniFORM, then positivity. Harmony is used if embeddings or clusters are involved. |
| Check | The residual-batch audit (Senior PR #18) runs after UniFORM. It compares each marker's positive fraction across groups within the same disease. A marker is examined further only if its calls differ by group beyond donor-to-donor spread. |
| Statistics | Donor-level, with acquisition group as a covariate where it can be estimated. |

## Known acquisition artifacts (documented; nothing to correct unless a call changes)

Measured 2026-09-28 on all 26 images. Tissue intensity was folded onto the 960 × 720 px camera field, which every image size is an exact multiple of.

- **E2 (Fusion 2.4.0), all 9 images: residual DAPI shading.**
  - About 10% top-to-bottom within every camera field; the centre is +4.1% vs the border.
  - It is absent in 2.3.1 and 2.4.1.
  - E-cadherin (−2.5%) and Vimentin (+3%) show weaker versions.
  - All 9 images share both Fusion 2.4.0 and one reference-capture date (2026-03-12), so the cause can't be separated.
- **E-cadherin, all versions:** about 8% top-to-bottom within fields. It is instrument-wide.
- **E1 (Fusion 2.3.1) metadata:** the files lack the `ImageProcessingApplied = OffsetAndShading` entry and use signal-units code 64 (69 elsewhere). The pixels are corrected anyway: DAPI is flat within fields (+0.5%) and the off-tissue background is at zero.
- **Scanner responsivity drift:** from the reference captures, over Sep 2025 to Sep 2026. This is instrument QC only.
