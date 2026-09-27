# v3 islets: endocrine objects (Drotar 2025 / ADIIDAS, CD99 primary)

The v3 islet step replaces the trained `Islets2` pixel classifier and its size-adaptive erode/dilate
with a rule-based, batch-comparable definition: endocrine pixels are **folds over each image's own
tissue background**; objects are those pixels merged where they touch, then sized by **mean caliper
diameter** into `Islet` (>= 35 um), `EndocrineCluster` (10-35 um) and `EndocrineSingle` (<= 10 um).
CD99 is the primary marker; INS, GCG and SST are the hormone union. IAPP (failed QC) and HLA-A
(disease-associated) are never used. PPY, CHGA and SYP are not in the panel.

| File | Role |
|---|---|
| `islet_anchors.py` | per-image, per-channel background anchors -> `results/v3/islet_anchors.csv` (tests: `tests/test_islet_anchors.py`) |
| `qupath/1_panc_objects_v3.groovy` | the detector (`ISLET_METHOD='drotar_cd99'`, default) or the legacy `Islets2` baseline |
| `qupath/4_islet_expander_v3.groovy` | names/UIDs (`Islet_N`, `Cluster_N`, `Single_N`), 20 um bands for Islets only, nearest-islet distances |
| `qupath/3_signed_distance_v3.groovy` | signed distances; always writes Islet, EndocrineCluster, EndocrineSingle, Artifact, pLN, Tissue columns |
| `qupath/8_hierarchy_parenting_v3.groovy` | parenting in the order Islet -> EndocrineCluster -> EndocrineSingle -> IsletExpanded -> Artifact -> pLN -> Tissue |
| `qupath/islet_benchmark_tiles.groovy` | blinded 1 mm tiles + key + clipped predictions for calibration |

The originals in `scripts/` are unchanged.

## 1. Anchors (`islet_anchors.py`)

Rule (Phenocycler_Analysis `INTENSITY_SCALE.md`): *standardize the unit, never the distribution.* An
anchor is one scalar per image and channel, estimated from that image in log space and applied as a
**linear divisor**. Nothing matches distributions between images, and no anchor is a quantile of the
endocrine population.

Per channel (CD99, INS, GCG, SST; index by name from `project.qpproj`), on the pyramid level nearest
2 um/px, inside the manual `Tissue` polygon(s) from `results/v3/manual_annotations/<tag>.geojson`:

- `n_pixels`, `zero_fraction`: the scanner subtracts offset and autofluorescence and stores UINT16, so
  negatives clamp to 0. **Zeros are excluded from every estimate** and reported, both at the 2 um level
  and sampled at full resolution (`zero_fraction_l0`: INS 0.20-0.81 across the cohort). The pyramid is
  block-averaged, which is why the 2 um level has fewer exact zeros.
- `bg_smo_mean`, `bg_smo_std`: a port of the `smo` 2.0.2 package that KINTSUGI `signal/smo.py` wraps,
  used as `SMO(sigma=1.0, size=7)` with threshold 0.05, then `bg_rv().mean()/.std()`. The operator
  sees the level image; zeros are dropped from the background **sample**, not from the operator input.
  SMO dilates every masked pixel to about a 15 x 15 block, so masking 30-80% of zeros would discard
  almost every pixel.
- `bg_logmode`: the mode of the log-intensity density of nonzero tissue pixels. Integer counts are
  spread over [k-0.5, k+0.5] first, so low counts are not a comb. When the density is bimodal (CD99 and
  GCG always; SST in some panels), the **upper substantial peak** is taken. The lower peak (about 5-13
  counts) is empty space inside the Tissue outline; the upper one is the exocrine tissue that endocrine
  signal must exceed.
- `bg_halfmax` = **`anchor`**: the right half-maximum of that peak in log space, which captures the
  tissue background level and its width together.
- `exposure_ms` (QPTIFF page XML `ExposureTime` is in microseconds), `responsivity`, `filter_date`,
  `signal_units`, and `biomarker_match` (the XML `<Biomarker>` against the project channel name: 96/96 match).

Why the half-maximum is the anchor. v2 islets were used only as a rough islet/exocrine reference mask.
On 8 images from 4 batches, with a 1.5 um blur at about 2 um, the median across images was:

| anchor candidate | H = max(INS,GCG,SST) at fold 3: exocrine positive % | v2-islet positive % | at fold 5: exo % | islet % |
|---|---|---|---|---|
| tissue log-mode | 9.1 (up to 26) | 99 | 2.9 | 97 |
| **right half-max** | **0.70** | **94** | **0.28** | **88** |
| tissue p90 | 0.66 | 94 | 0.26 | 88 |
| tissue p95 | 0.29 | 91 | 0.13 | 84 |

- **The mode ignores background width.** In the panc-batch panel, SST (channel 1) has a broad
  exocrine background: mode 5 counts, but the exocrine p90 is about 7x the mode.
- **Half-max matches p90 in performance** while remaining a property of the background peak rather
  than a quantile of the whole population. It is also scale-equivariant and moves by less than 3%
  when 3% bright pixels are added (unit tests).
- **SMO is recorded but not used.** At 2 um it accepts only 1-2% of tissue pixels, roughly uniformly
  across intensity, so its mean contains islet signal (6450panc INS: mean 21 vs median 3). At full
  resolution it does discriminate, but what it isolates is the empty-space noise floor.

Anchors vary 4-9x between images per channel. Exposure differs at most 1.9x, and anchors spread as
widely within a filter epoch as between epochs (2025-09-10, 2026-03-12, 2026-06-11, plus 2026-08-19 on
the two metabolic-panel slides). An exposure correction alone would not equalize the images.

Runtime: 7-40 s and about 2-10 GB per image. All 24 images took about 7 min. 6457 has no manual export
(its `.qpdata` is empty), so it was run with
`--tissue-geojson results/Geojson/Tissue/Tissue__6457_...geojson`.

```bash
/scratch/io60/venv-v3/bin/python scripts/v3/islet_anchors.py --all            # upserts (image, channel) rows
/scratch/io60/venv-v3/bin/python scripts/v3/islet_anchors.py --images 6450panc --out /scratch/x.csv
/scratch/io60/venv-v3/bin/python -m pytest -q scripts/v3/tests/test_islet_anchors.py   # 18 tests, synthetic only
```

## 2. Detector (`qupath/1_panc_objects_v3.groovy`, `ISLET_METHOD='drotar_cd99'`)

1. **Pixels, via ImageOps at `RESOLUTION_UM` (1.0156 um, pyramid level 1).** Each channel is divided by
   its anchor, then blurred with `Gaussian(SIGMA_UM)`, giving:
   - `C = CD99 fold >= TAU_CD99`
   - `H = max(INS, GCG, SST) fold >= TAU_H`
   - `INS`, `GCG` and `SST` fold `>= TAU_H` each, for composition.

   **One pixel pass.** A float data server carries one channel per (mask, threshold), and the contour
   tracer applies the thresholds. A multi-channel *probability* server would be wrong: `ContourTracing`
   collapses it to an argmax, but C and H must overlap. `ensureType(FLOAT32)` is required because the
   op server otherwise declares UINT16 and truncates folds. The same ops are saved per image as pixel
   classifiers in `OUT_CLASSIFIER_DIR/<tag>/` for review in the GUI:
   - `drotar_cd99_{C,H,INS,GCG,SST}.json`
   - `drotar_cd99_routes.json` (1 = CD99 only, 2 = hormone only, 3 = both)
2. **Geometry (JTS), inside Tissue minus `EXCLUDE_CLASSES` (pLN, Artifact)**, clipped 5 um inside the
   image edge:
   - `C_only = C - H`, then opened by `OPEN_R_UM`, **CD99-only pixels only**. This is the bridge
     guard: it removes CD99+ strands and scattered CD99-bright exocrine cells that could join islets.
   - `U = H + open(C_only)`. `C & H` is already inside H.
   - Parts sharing an **edge** form one object. A shared corner, i.e. diagonal pixels, does not join
     them: this is the tracer's own 4-connectivity.
   - Holes are filled. A part swallowed by a filled hole merges into its host.
   - Objects smaller than `MIN_AREA_UM2` are dropped.
   - **No erode/dilate:** it would delete single cells.
   - Optional: `OPEN_H_R_UM > 0` also opens hormone pieces. It is off by default; see limitations.
3. **Class** by `Mean caliper diameter um = (Max diameter + Min diameter) / 2`, using QuPath shape
   features in um.
4. **Measurements:**
   - `Mean caliper diameter um`, `Equivalent diameter um`, `Circularity`, `Solidity`, `Area um^2`
   - `Endocrine route` (1 = CD99 only, 2 = hormone only, 3 = both; at least one classifier pixel each)
   - `CD99+`, `INS+`, `GCG+`, `SST+` and `Hormone+ area um^2`
   - `CD99-only fraction`
   - `Bridge flag`: area >= `BRIDGE_AREA_UM2` and CD99-only fraction >= `BRIDGE_CD99ONLY_FRAC`
   - `Tissue border flag` (within `BORDER_UM` of the analysis-area edge; flagged, not removed) and
     `Distance to tissue border um`
5. **Image sanity:** endocrine area as a fraction of tissue. It is printed, written to
   `<WRITE_GEOJSON>.summary.json` and to run_v3's `V3_RESULT`, and flagged above 10% (ADIIDAS).
6. **Replacement and saving:**
   - Everything is computed first. Only then are the old Islet, EndocrineCluster and EndocrineSingle
     objects removed, along with their orphaned bands; their detections move up. A failure therefore
     leaves the hierarchy unchanged.
   - The script never saves. The run_v3 footer, or the user, decides.
   - No top-level `return`, so the provenance footer is always reached.

Legacy `ISLET_METHOD='Islets2'` runs the original code unchanged. It also removes drotar
clusters and singles, so the two methods never coexist.

## 3. Scripts 4, 3, 8

- **4:**
  - Two-pass stable UIDs per class. An object carrying another class's automatic name, e.g. a
    reviewed cluster reclassified as Islet, is renamed.
  - Voronoi 20 um bands are built for `Islet` only.
  - `Nearest islet centroid distance um`: centroid to centroid, among islets within +-3 mm in x and y;
    NaN if there are none.
  - `Nearest islet edge distance um`: boundary to boundary, within 3 mm.
  - The inherited outline simplification lets bands reach about 23 um (tolerance 3.3 um).
- **3:** unchanged except that the six columns above always exist (NaN when a class is absent), and
  a missing image or missing detections throw instead of returning silently.
- **8:** uses the explicit order. Other classes (Nerve, Blood_vessel, Lymphatic, Gal3Norm) are not
  used for parenting; their annotation and covered-cell counts are printed. A cluster inside a band
  keeps its cells.

Fixture chain 1 -> 4 -> 3 -> 8 (2048^2 crop of 6450, 13,516 cells): all checks pass. The checks were:
- UIDs unique and equal to the name numbers
- one band per Islet with no orphans
- nearest distances equal to brute force
- six signed-distance columns present
- the parent class agrees with the sign of the signed distance
- classes agree with the caliper diameter
- reruns of 1 and 4 are idempotent

## 4. Parameters (script 1)

| Parameter | Default | Meaning |
|---|---|---|
| `ISLET_METHOD` | `drotar_cd99` | or `Islets2` (legacy baseline) |
| `TAU_CD99` | 3 | CD99 fold threshold (provisional) |
| `TAU_H` | 5 | hormone fold threshold on max(INS, GCG, SST) (provisional) |
| `SIGMA_UM` | 1.5 | Gaussian sigma before thresholding |
| `OPEN_R_UM` | 4 | opening radius for CD99-only pixels (removes structures narrower than about 8 um) |
| `OPEN_H_R_UM` | 0 | optional opening of hormone pixels |
| `MIN_AREA_UM2` | 20 | smaller parts dropped |
| `RESOLUTION_UM` | 1.0156 | classifier resolution |
| `ANCHORS_CSV` | `results/v3/islet_anchors.csv` | |
| `OUT_CLASSIFIER_DIR` | `classifiers/generated` | per-image classifier JSONs (`SAVE_CLASSIFIERS=false` to skip) |
| `BRIDGE_AREA_UM2`, `BRIDGE_CD99ONLY_FRAC` | 5000, 0.5 | bridge flag |
| `BORDER_UM` | 10 | tissue-border flag |
| `EXCLUDE_CLASSES` | `pLN,Artifact` | cut from the analysis area |
| `THREADS` | 0 | tile workers; 0 = run_v3's ForkJoin parallelism. 16 workers: 48-67 GB RSS; the default 63: about 100 GB |
| `WRITE_GEOJSON` | '' | write objects plus `<name>.summary.json` |
| `TAU_CD99_GRID`, `TAU_H_GRID`, `OPEN_R_GRID`, `OPEN_H_GRID` | '' | sweeps (below) |

Pass parameters as `-a KEY=VALUE` script arguments; an unknown key is an error. run_v3 can also
override them by rewriting a single `<type> NAME = ...` declaration.

## 5. Calibration

```bash
# one pixel pass, every tau x opening combination; hierarchy untouched; objects + summary per combination
QuPath script -p <worker project> -i "<image>" \
  -a WRITE_GEOJSON=/scratch/io60/calib/<tag>/s1.5.geojson \
  -a TAU_CD99_GRID=2,3,5 -a TAU_H_GRID=3,5,7,10 -a OPEN_R_GRID=3,4,5 -a OPEN_H_GRID=0,2 \
  -a SIGMA_UM=1.5 scripts/v3/qupath/1_panc_objects_v3.groovy
# SIGMA_UM (and RESOLUTION_UM) change the pixels: one run per value (-a SIGMA_UM=1.0 / 2.0 / 3.0)
```

- **Cost:** about 25-70 s per image and sigma for 36 combinations. The sweep writes
  `<base>__tC.._tH.._r..[_rH..].geojson` plus `<base>.summary.json` (counts by class, route and
  class x route; bridge and border flags; endocrine fraction; timings).
- **Benchmark tiles:** `islet_benchmark_tiles.groovy`, one run per image. `OBJECTS_GEOJSON` accepts a
  sweep output, so tiles can be cut without saving any `.qpdata`. Tile kinds:
  - random tiles inside Tissue
  - tiles centred on detector objects, round-robin over class x route
  - `exo_cd99` tiles: the most CD99-bright pixels with no Islet and <= 0.5% hormone-positive area,
    the bridge test
- **Tile outputs:**
  - `tiles/<code>.ome.tif`: DAPI, CD99, INS, GCG and SST as raw counts. The file is checked to contain
    no source name.
  - `pred/<code>.geojson`: objects in tile coordinates.
  - `key.tsv`: the only link back to the image.
- **Scoring** (elsewhere): per-stratum macro-F1 per route and size class against annotators, matching
  each sweep's GeoJSON to the tiles through `key.tsv`.
- **Tuning:** choose `TAU_CD99`, `TAU_H`, `SIGMA_UM` and `OPEN_R_UM` by macro-F1 across strata
  (epoch x disease x size), with the Islets2 baseline alongside.

**Provisional values** (tau_CD99 = 3, tau_H = 5, sigma = 1.5 um, r = 4 um) come from two images:

| image (v2 Islets) | tau_H | Islet | Cluster | Single | route 1/2/3 | endocrine % tissue | v3 Islets on a v2 islet | v2 islets touched |
|---|---|---|---|---|---|---|---|---|
| 6450panc (453) | 2 | 535 | 8618 | 10210 | 1/18635/727 | 5.8 | 46% | 100% |
| | 3 | 288 | 3178 | 3876 | 2/6512/828 | 3.0 | 78% | 100% |
| | **5** | **192** | **1273** | **1200** | **10/1778/877** | **1.9** | **95%** | **100%** |
| | 10 | 141 | 739 | 536 | 36/565/815 | 1.3 | 98% | 100% |
| 6533 metabolic (1596) | 2 | 802 | 8317 | 9484 | 4/17900/699 | 5.0 | 47% | 98% |
| | 3 | 514 | 5577 | 6056 | 7/11360/780 | 3.2 | 62% | 97% |
| | **5** | **287** | **3144** | **3644** | **14/6115/946** | **1.9** | **91%** | **93%** |
| | 10 | 155 | 1409 | 2091 | 112/2168/1375 | 1.0 | 99% | 74% |

- **tau_H dominates.** At tau_H 2-3, hormone noise gives 3-6% endocrine area and about 20k specks.
  tau_H = 5 gives 1.9% on both images; v2 islets covered 1.8% of 6450panc. 91-95% of v3 Islets then
  sit on a v2 islet, and 93-100% of v2 islets are touched by some v3 object.
- **tau_CD99 changes counts by less than 3%** at r = 4 um; it mainly moves objects between routes 2
  and 3. The opening removes most CD99-only rims.
- **Bridge flags were 0** at every tau_H >= 3.
- **sigma mainly removes singles and clusters.** At tau_H = 5, going from sigma 1.0 to 3.0 takes
  6450panc singles from 1432 to 446, while Islets stay at about 195.
- **v2 counts are not comparable one to one.** 64% (6450) and 83% (6533) of v2 "Islets" are smaller
  than a 35 um disc, so they are clusters or singles under Drotar's rule.

## 6. Known limitations

- **A linear fold standardizes scale, not offset.** When the scanner's subtraction overshoots (INS: up
  to 81% zeros at full resolution), the anchor measures noise width rather than a level. This
  additive batch effect cannot be removed by a divisor; the benchmark strata must cover such images.
- **The anchor is estimated at 2 um but applied at 1 um after the blur.** This is consistent for block
  means but not for noise tails.
- **Hormone noise dominates the small objects.** Most singles and clusters are hormone-only (route 2).
  In the panc panel, SST (channel 1) has a broad, structured exocrine background. In the metabolic
  panel, INS (channel 1, 100 ms) shows exocrine speckle. Both produce small false objects at tau_H = 5.
  Singles are the least reliable class.
- **Hormone-positive bridges are not guarded.** On 6533 at tau_H = 5, the largest object (68,000 um^2)
  spans 3 v2 islets. The CD99-only opening cannot split it, and `OPEN_H_R_UM` = 1.5-2.5 um did not
  either (66,000 um^2) while it shrank the Islet count by 6-14%. It separates only at tau_H = 10.
  Whether this is one islet or two is for the benchmark to decide.
- **CD99 contributes little to object extent** at r = 4 um. CD99-only (route 1) objects are rare
  (10-50 per image), so hormone-negative islet remnants in T1D are found only when their CD99-only area
  is at least about 8 um thick.
- **`exo_cd99` tiles** are defined by no Islet and <= 0.5% hormone-positive area, because strictly
  "none" selects almost nothing.
- **Detections are not re-parented** by script 1 (script 8 does it), so between steps 1 and 8, cells
  sit under Tissue.
- **The 6450panc project entry** was fixed from 0.5077663810 to 0.5078 um on 2026-09-27 at 12:18, and a
  backup `project.qpproj.pre_pixelfix_20260927` exists. The first benchmark tiles were cut from a
  pre-fix copy, so their OME PhysicalSize is 0.5077664 (0.004% off; pixels are identical).
- **pLN is disjoint from Tissue** in all four images that have one (6533, HDL112, 6551, HDL115: 0%
  overlap), so run_v3's preflight check "pLN inside Tissue" will fail on them. The detector excludes
  pLN either way.

## 7. Decisions for you

1. Accept the half-max anchor, and SMO as a diagnostic only, or ask for exposure x responsivity
   scaling as an alternative unit.
2. Choose per-hormone thresholds (for example a separate tau for SST) or keep one tau_H for the union,
   if the benchmark shows SST- or INS-driven false singles.
3. Decide whether hormone-only singles count as endocrine singles, or whether singles should require
   CD99 co-positivity (route 3).
4. Decide whether 6533's merged 68k um^2 object is one islet or two; this sets how strict tau_H must be.
5. Set the benchmark size (N per kind per image; the default is 8) and the annotator protocol for the
   tile codes.
