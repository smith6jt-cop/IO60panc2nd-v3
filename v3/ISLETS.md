# v3 islets: endocrine objects (Drotar 2025 / ADIIDAS, CD99 primary)

The v3 islet step replaces the trained `Islets2` pixel classifier and its size-adaptive erode/dilate with a
rule-based definition that works across donors and disease: a pixel is endocrine when it is positive in CD99
(primary) or in the hormone union INS/GCG/SST, **relative to that image's own tissue background**; objects are
those pixels merged where they touch, sized by **mean caliper diameter** into `Islet` (>= 35 um),
`EndocrineCluster` (10-35 um) and `EndocrineSingle` (<= 10 um). IAPP (failed QC) and HLA-A (disease-associated)
are never used. PPY, CHGA and SYP are not in the panels.

Batch (user ruling 2026-09-28): batch is corrected at the expression level (UniFORM; Harmony for embeddings);
nothing gates segmentation. The per-image background reference below IS the islet method, not a batch gate.
Strata for checking it are acquisition group x disease (E1 Fusion 2.3.1 panel 1; E2 2.4.0 panel 1; E3 2.4.1
panel 2, no ND donor; E4 metabolic 70-channel panel). The qptiff Responsivity date is not a batch.

| File | Role |
|---|---|
| `islet_anchors.py` | per-image, per-channel background anchor + censored clamp model + annotation hash -> `results/v3/islet_anchors.csv` (tests: `tests/test_islet_anchors.py`) |
| `qupath/1_panc_objects_v3.groovy` | the detector (`ISLET_METHOD='drotar_cd99'`, default) or the legacy `Islets2` baseline |
| `qupath/4_islet_expander_v3.groovy` | names/UIDs, nearest-islet distances, `Artifact overlap um^2`; **no bands by default** (`MAKE_BANDS = false`: peri-islet regions dropped, 2026-09-29) |
| `qupath/3_signed_distance_v3.groovy` | signed distances; always writes Islet, EndocrineCluster, EndocrineSingle, Artifact, pLN, Tissue columns |
| `qupath/8_hierarchy_parenting_v3.groovy` | parenting: Islet -> EndocrineCluster -> EndocrineSingle -> Artifact -> pLN -> Tissue |
| `qupath/islet_benchmark_tiles.groovy` | blinded, quota-based benchmark tiles (design v2) + private key and predictions |
| `qupath/islet_benchmark_report.py` | benchmark coverage (category x stratum, gaps) and the tile leak check |

The originals in `scripts/` are unchanged.

## 1. Anchors and clamp models (`islet_anchors.py`, method 2.0)

Rule (Phenocycler_Analysis `INTENSITY_SCALE.md`): *standardize the unit, never the distribution.* Everything is a
property of the image's own tissue BACKGROUND, estimated in the detector's pixel space, and applied as a linear
scale (the anchor) or an additive floor (the clamp model). Nothing matches distributions between images, and
nothing is a quantile of the endocrine population.

**Pixel space = the detector's.** Pyramid level 1 (1.0156 um) blurred with the detector's Gaussian (sigma 1.5 um,
4-sigma kernel), so the anchor and the detector see the same values. Statistics use a seeded random sample of at
most 4 M mask pixels; the blur runs on the whole crop. (v1 estimated at ~2 um without blur: the effective INS tau
differed by up to 21% between images.)

**Analysis mask = the detector's analysis area:** union of the Tissue polygons minus the union of pLN and Artifact,
5 um inside the image edge. Polygons are rasterised one at a time and OR-ed (v1 filled every ring even-odd, so
overlapping Tissue polygons cancelled).

**Occupied tissue (added with the clamp model; flagged).** The Tissue outline encloses empty space (lumens, fat,
gaps: about half the outline on 6450panc and HDL112). A pixel counts as background tissue when DAPI blurred at
5 um exceeds the valley between the empty-space and tissue peaks of its log density (everything when DAPI is
unimodal: 9 of 24 images; otherwise 51-98% of the mask), and it is not "no data" (all four channels exactly 0). Without it the level-0 clamp
fit describes empty space (6450panc CD99 mu 5 instead of 93 counts, GCG 4 instead of 102) and the unimodal
INS/SST anchors drop 2-3x (6450panc INS 7.9 vs 16.6). CD99/GCG anchors do not change (their upper-peak rule
already skips empty space).

**Two passes (halo).** Pass 1 estimates every channel on all occupied pixels. The detector rule at tau_CD99 3,
tau_H 5, k 3 marks endocrine pixels, which are dilated by 25 um; pass 2 re-estimates on occupied pixels outside
that mask. Signal spill-over around islets fattens the background peak's right shoulder: over the 96 (image,
channel) pairs pass 2 is 3% lower at the median, >= 12% lower for a tenth, up to 26% (6534 SST 190 -> 140); 12-32%
of the mask is excluded. `anchor_pass1` is kept.

**The anchor** (`anchor` = `bg_halfmax`) is unchanged in definition: the right half-maximum of the tissue
background peak of the log density of the blurred values (the UPPER substantial peak when bimodal). It is
scale-equivariant and ignores a few percent of bright pixels.

**The clamp model (finding 1).** The scanner subtracts offset + autofluorescence and clamps negatives to 0: up to
84% of level-0 INS tissue pixels are exactly 0 (6593, 6533). Clamping shrinks the half-max anchor, so the INS anchor
tracks the clamp fraction (Spearman -0.96 over 24 images) and tau x anchor falls into the noise on the most clamped
images. The latent background behind the clamped counts is modelled as N(mu, sigma^2) in LINEAR counts, fitted by
maximum likelihood to level-0 integers from 64 random 512 px tiles of occupied tissue (pass 2), with **zeros
left-censored** (P(X = 0) = Phi((0.5 - mu)/sigma)) and the likelihood truncated above the median of the positives,
so only the zero mass and the lower half of the positives enter (no signal, no speckle). Recorded: `clamp_mu`,
`clamp_sigma`, `clamp_fraction`, `clamp_fit_ks`, `clamp_converged`, pass-1 values. A fit is NaN (floor off for that
channel, with a warning in the detector) when > 95% of values are zero or |ECDF - model| > 0.02 twice: CD99 on 3
images (0.2-0.5% zeros), nothing else.
- *Linear, not log.* The level-0 histograms are flat across 1-8 counts right after the zero spike (6533 INS:
  80% zeros, then 0.83, 0.81 ... 0.66% per count): a normal whose mean lies BELOW zero. A zero-censored log-normal
  fitted the same way runs to its parameter bounds with KS 0.035-0.22, against <= 0.018 for the linear model.
- Scale-equivariant (clamping at 0 is invariant to scale).

**Annotation hash (finding 6).** `annotation_hash` = SHA-256 (16 hex) over `classes=<sorted Tissue,pLN,Artifact>`
and one sorted line per polygon, `class|x,y;x,y;.../<hole ring>`, every vertex at 2 decimals, the precision of
QuPath's GeoJSON export. The detector serialises its current annotations with the same exporter, parses the text
and builds the same lines, so no floating-point arithmetic is compared. `--skip-existing` recomputes an image
whose hash or estimation parameters (`params_hash`) changed.

Diagnostics kept: `anchor_l2_unblurred` (the v1 estimator on the new mask), SMO at ~2 um (`bg_smo_*`, not
scale-robust; see v1 notes in git history), exposure / responsivity / filter date / biomarker check from the
QPTIFF page XML (96/96 biomarkers match).

Runtime (2026-09-28): all 24 images in 16.5 min (23 in 944 s + 6457 in 48 s; 21-71 s per image), max RSS 13.7 GB.

```bash
/scratch/io60/venv-v3/bin/python scripts/v3/islet_anchors.py --all --skip-existing   # recomputes changed exports only
/scratch/io60/venv-v3/bin/python -m pytest -q scripts/v3/tests/test_islet_anchors.py   # 35 tests, synthetic only
```

## 2. Detector (`qupath/1_panc_objects_v3.groovy`, `ISLET_METHOD='drotar_cd99'`)

1. **Anchor check.** The image's CSV rows must carry the hash of the current Tissue/pLN/Artifact annotations, the
   detector's classes (`EXCLUDE_CLASSES`), `SIGMA_UM` and `RESOLUTION_UM`. Any difference stops the run with the
   command to fix it (`ANCHOR_CHECK='error'`, default) or is printed as a loud warning (`'warn'`). A method-1.x
   CSV (no clamp columns) is a mismatch too.
2. **Pixels** (ImageOps at 1.0156 um; x = Gaussian(1.5 um) of the raw counts). Channel c is positive when
   `x_c >= max(tau * anchor_c, clamp_mu_c + FLOOR_K * clamp_sigma_c)`; tau = `TAU_CD99` for CD99, `TAU_H` for INS,
   GCG and SST. `C` = CD99 positive, `H` = any hormone positive (one traced channel: max over hormones of
   fold/threshold >= 1). Thresholds are traced by the contour tracer, one float channel per (mask, threshold), so a
   sweep reads the pixels once.
3. **Objects** (JTS, inside Tissue minus pLN and Artifact, 5 um inside the image edge): `U = H + (C - H)`, parts
   that share an edge are one object (a shared corner is not), holes filled, parts swallowed by a filled hole
   merged, < `MIN_AREA_UM2` dropped. **No opening of CD99-only pixels:** CD99-only cell-sized blobs are endocrine
   candidates (the orchestrator's 4 um "bridge guard" deleted 99.7% of the CD99-only parts; removed). No erode/dilate.
4. **Optional safeguards, all OFF by default** (named, sweepable; the benchmark decides):
   - `NECK_SPLIT` (`NECK_R_UM` 4, `CORE_MIN_AREA_UM2` 200): in an object with >= 2 hormone cores (connected
     hormone+ parts >= 200 um^2), the CD99-only pixels removed by an opening of radius 4 um are candidate necks;
     if removing them separates the cores, the object is split, thin parts touching one piece go back to it and
     those touching >= 2 pieces (the necks) are dropped. CD99-only regions thicker than 8 um are never cut.
   - `SINGLE_COPOS` (`CLAMP_HEAVY_FRAC` 0.6): an EndocrineSingle positive ONLY in one hormone channel whose
     level-0 clamp fraction is >= 0.6 (no CD99, no second hormone) is dropped. Heavily clamped at 0.6: INS on
     11/24 images (E1 6521, 6380, 6579; E2 6442, 6461, 6436, 6593; E3 6450, 6523, 6534; E4 6533), SST on 6523.
   - `OPEN_H_R_UM`: opening of the hormone pieces (kept from v1; did not split the 6533 merge).
5. **Class** by `Mean caliper diameter um = (Max diameter + Min diameter) / 2` with QuPath's shape features:
   **`Max diameter` = the diameter of the minimum bounding circle (not the longest chord); `Min diameter` = the
   minimum caliper width** (JTS `MinimumBoundingCircle` / `MinimumDiameter`, read from the 0.7 bytecode). The JTS
   fallback uses the same two definitions. Islet >= 35 um; EndocrineCluster > 10 and < 35 um; EndocrineSingle
   <= 10 um.
6. **Measurements:** Mean caliper diameter, Equivalent diameter, Area, Circularity, Solidity, Max/Min diameter,
   `Endocrine route` (1 CD99 only, 2 hormone only, 3 both; >= one classifier pixel each), CD99+/INS+/GCG+/SST+/
   Hormone+ area, CD99-only fraction, **`Hormone cores`** and **`Multi-core flag`** (>= 2 cores: a possible merge
   of separate islets through CD99-only tissue), `Neck split`, `Positive channels` (0-4), `Clamped-only flag`,
   Bridge flag (>= 5000 um^2 and CD99-only fraction >= 0.5), Tissue border flag (<= 10 um from the analysis-area
   edge; flagged, not removed) and Distance to tissue border.
7. **Image sanity:** endocrine area / tissue area, flagged above 10% (ADIIDAS); summary JSON and run_v3 V3_RESULT.
8. **Replacement:** everything is computed first; only then are old Islet/EndocrineCluster/EndocrineSingle objects
   and their IsletExpanded bands removed (their detections move up). **Detections left directly under the root**
   (bands live at the root; 4,307 of 13,516 fixture cells) are re-parented to the Tissue (else pLN) annotation
   containing their nucleus centroid. The script never saves (the run_v3 footer or the user decides); no top-level
   `return`, so the provenance footer is always reached. The per-image review classifiers
   (`drotar_cd99_{C,H,INS,GCG,SST,routes}.json`) use the same thresholds.

Legacy `ISLET_METHOD='Islets2'` runs the original code unchanged and removes drotar clusters and singles, so the two
methods never coexist.

## 3. Scripts 4, 3, 8

- **4:** two-pass stable UIDs per class; an object carrying another class's automatic name is renamed. **Keys of
  another class are removed** after a reclassification ('Islet: UID', 'Cluster: UID', 'Single: UID' only on their
  own class; the nearest-islet distances only on Islets). **No peri-islet bands** (user, 2026-09-29: peri-islet
  regions dropped entirely; the per-cell signed distances of script 3 carry that relation). `MAKE_BANDS = true`
  restores 20 um Voronoi bands for `Islet` only, clipped to Tissue minus (Artifact and pLN), with 'Islet: UID' on
  each band; old IsletExpanded bands are always removed. `Nearest islet centroid distance um` (Drotar,
  within +-3 mm) and `Nearest islet edge distance um`. **`Artifact overlap um^2`** on every Islet, cluster and single
  (area inside Artifact annotations; objects drawn or kept over an Artifact at review can be excluded downstream).
  Detections that the band removal left under the root are re-parented to Tissue/pLN.
- **3:** unchanged; the six columns above always exist (NaN when a class is absent).
- **8:** order **Islet -> EndocrineCluster -> EndocrineSingle -> Artifact -> pLN -> Tissue** (no IsletExpanded since
  2026-09-29), tested at the nucleus centroid as in script 3; an endocrine object
  keeps its cells even over an Artifact (and carries `Artifact overlap um^2`).

Fixture chain (2048^2 crop of 6450 with a pLN inside Tissue, 13,516 cells; `/scratch/io60/work/isletC2/fixture`):
1 -> 4 -> 3 -> 8, all checks pass, including: no detection under the root after script 1 (was 4,307 from the old
bands), anchor hash equal in Python and Groovy, bands inside Tissue and never on pLN/Artifact, pLN cells parented to
pLN, and, after an Artifact drawn over an islet and a reclassification (Cluster -> Islet, Islet -> Cluster):
`Artifact overlap` equal to the geometric overlap on 8 objects, 235 cells parented to Artifact, no foreign UID or
nearest-islet key; script 1 then refuses (hash changed), runs with `ANCHOR_CHECK=warn`, and passes after
re-export + `islet_anchors.py --skip-existing` (which recomputes only that image).

## 4. Parameters (script 1)

| Parameter | Default | Meaning |
|---|---|---|
| `ISLET_METHOD` | `drotar_cd99` | or `Islets2` (legacy baseline) |
| `TAU_CD99` | 3 | CD99 fold (x / anchor) threshold (provisional) |
| `TAU_H` | 5 | hormone fold threshold, INS, GCG and SST each (provisional) |
| `FLOOR_K` | **4** | clamp floor `clamp_mu + k * clamp_sigma` (0 = off; provisional, section 6) |
| `SIGMA_UM`, `RESOLUTION_UM` | 1.5, 1.0156 | must equal the anchors' (checked) |
| `MIN_AREA_UM2` | 20 | smaller parts dropped |
| `NECK_SPLIT`, `NECK_R_UM`, `CORE_MIN_AREA_UM2` | **false**, 4, 200 | optional targeted neck split |
| `SINGLE_COPOS`, `CLAMP_HEAVY_FRAC` | **false**, 0.6 | optional co-positivity for singles in heavily clamped channels |
| `OPEN_H_R_UM` | 0 | optional hormone opening |
| `ANCHORS_CSV`, `ANCHOR_CHECK` | `results/v3/islet_anchors.csv`, `error` | anchors; mismatch policy (`warn`) |
| `EXCLUDE_CLASSES` | `pLN,Artifact` | cut from the analysis area (and from the anchors' mask) |
| `BRIDGE_AREA_UM2`, `BRIDGE_CD99ONLY_FRAC`, `BORDER_UM` | 5000, 0.5, 10 | flags |
| `OUT_CLASSIFIER_DIR`, `SAVE_CLASSIFIERS` | `classifiers/generated`, true | per-image review classifiers |
| `THREADS` | 0 | tile workers (0 = run_v3's ForkJoin parallelism); 16 workers: 68-97 GB RSS |
| `WRITE_GEOJSON` | '' | objects + `<name>.summary.json` |
| `TAU_CD99_GRID`, `TAU_H_GRID`, `FLOOR_K_GRID`, `NECK_SPLIT_GRID`, `SINGLE_COPOS_GRID`, `OPEN_H_GRID` | '' | sweeps |

Pass `-a KEY=VALUE` script arguments (an unknown key or a non-boolean for a boolean is an error), or let run_v3
rewrite a declaration.

```bash
# one pixel pass, every combination; hierarchy untouched; <base>__tC.._tH.._k..[_rH..][_neck][_copos].geojson + summary
QuPath script -p <worker project> -i "<image>" -a WRITE_GEOJSON=/scratch/io60/calib/<tag>/s.geojson \
  -a FLOOR_K_GRID=0,3,4,5 -a NECK_SPLIT_GRID=false,true -a SINGLE_COPOS_GRID=false,true -a TAU_H_GRID=4,5,7 \
  scripts/v3/qupath/1_panc_objects_v3.groovy
# SIGMA_UM / RESOLUTION_UM: recompute anchors in that pixel space first (islet_anchors.py --sigma-um 2 --out x.csv),
# then -a SIGMA_UM=2 -a ANCHORS_CSV=x.csv (the anchor check enforces this)
```

## 5. Benchmark (design v2: `islet_benchmark_tiles.groovy`, `islet_benchmark_report.py`)

Aim: enough REPRESENTATIVE positives and negatives in every stratum (group x disease), not a tile count.

- **Categories** (a quota of `QUOTA_PER_STRATUM` = 3 tiles per category per stratum, split over the stratum's
  images, at least 1 per image when its pool is not empty):
  - positives: `pos_islet_small` (35-60 um), `_medium` (60-120), `_large` (>= 120), `pos_ins_rich` (INS >= 50% of
    hormone area), `pos_ins_deficient` (INS < 5%: GCG/SST-dominant), `pos_infiltrated` (>= 10% of the islet's 20 um
    periphery immune-bright: CD3e/CD8/CD20 above the image's own p99; no CD45 in these panels), `pos_cluster`,
    `pos_single`, `pos_cd99_only` (route 1), `pos_hormone_only` (route 2), `pos_legacy_only` (a legacy Islets2 islet
    no v3 setting found);
  - hard cases: `hard_near_threshold` (found by some but not all v3 settings);
  - negatives (window grid, top 10% by the rule): `neg_exocrine` (random, no islet/cluster candidate),
    `neg_exo_cd99` (most CD99+ single candidates), `neg_speckle` (most hormone-only single candidates; the dominant
    channel is recorded: SST in panel 2, INS in the metabolic panel), `neg_duct` (Pan-Cytokeratin-bright),
    `neg_vessel` (CD34-bright), `neg_immune` (immune-dense, no islet/cluster), `neg_edge` (tissue 30-70% of the
    tile), `neg_fold` (bright in every channel).
- **Candidates** = the UNION of several detector settings (`CANDIDATES`, e.g. default k4, k0, k3, neck, copos
  outputs of one sweep), the legacy Islets2 islets (`results/v3/annotations_all/<tag>.geojson`), the structure rules
  and random draws. Overlapping candidates of different sources are one group (one pool element). Each tile records
  its pool size, n drawn and the **inclusion probability** (n / pool size; seeded shuffle); rare categories are drawn
  first; tiles of one image overlap <= 25%.
- **Blinding:** `OUT_DIR/annotate/tiles/<code>.ome.tif` is the only annotator-visible output (400 um tiles; DAPI,
  CD99, INS, GCG, SST raw counts; fixed channel names and **fixed palette**; pixel size **0.5078 for every image**;
  OME name = code; one fixed file mtime). Tile centres are **jittered** (+-30% of the tile for objects, the target
  kept >= 10% inside; windows by half a grid step). Everything else is in `OUT_DIR/private/` (never share):
  `key.tsv` (code -> image, x, y, w, h, tile_kind = category, stratum, inclusion probability...), `runs/<tag>.json`,
  `pred/<label>/<code>.geojson` (the scorer's `--pred-dir label=...`; objects only) and `regions/<code>.geojson`
  (Tissue/pLN/Artifact outlines, kept out of `pred/` because the scorer reads every feature there).
- **Leak check** (`islet_benchmark_report.py`): annotate/ holds tiles only; every TIFF tag string of every page is
  scanned for image names, tags, stems, donor IDs, 'qptiff', 'Scan<N>', 'panc', 'pLN', 'HDL' and other pixel sizes;
  raw bytes for the long identifiers; OME name = code, PhysicalSize = 0.5078, channel names/colours identical in all
  tiles, no acquisition/instrument elements, XResolution consistent, one tile size, one mtime; key <-> files.
- **Coverage report** (`private/coverage.tsv`): tiles per category x stratum against the quota; gaps flagged
  `NOT_RUN` (stratum not tiled yet), `POOL_EMPTY` (no candidates in the tiled images) or `SHORT`.
- **Scoring** (PR #1 scorer): per stratum and category, macro-F1 per route and size class against the annotators,
  every candidate setting scored from `private/pred/<label>`; weight by 1/inclusion probability for prevalence.

```bash
QuPath script -p <worker project> -i "<image>" -a OUT_DIR=/scratch/io60/benchmark_v2 \
  -a CANDIDATES=default=<sweep>__tC3_tH5_k4.geojson,k0=<sweep>__tC3_tH5_k0.geojson,neck=<...>_neck.geojson \
  scripts/v3/qupath/islet_benchmark_tiles.groovy
/scratch/io60/venv-v3/bin/python scripts/v3/qupath/islet_benchmark_report.py --out-dir /scratch/io60/benchmark_v2
```

## 6. Results (2026-09-28, private worker copies; nothing saved to the project)

**Detector variants** (read-only sweeps, one pixel pass per image; tau_CD99 3, tau_H 5; `route 1/2/3` = CD99 only /
hormone only / both; "INS-only singles" = singles positive in INS and no other channel, per mm^2 of analysis area):

| image | run | Islet / Cluster / Single | route 1/2/3 | CD99-only obj | INS-only singles /mm2 | endocrine % | multi-core (necks cut) | copos dropped |
|---|---|---|---|---|---|---|---|---|
| 6476 (E1 ND) | default (k4) | 972 / 2126 / 1415 | 145/2647/1721 | 145 | 4.13 | 3.64 | 3 | 0 |
| 6476 (E1 ND) | k0 (no floor) | 972 / 2126 / 1415 | 145/2647/1721 | 145 | 4.13 | 3.64 | 3 | 0 |
| 6476 (E1 ND) | k3 | 972 / 2126 / 1415 | 145/2647/1721 | 145 | 4.13 | 3.64 | 3 | 0 |
| 6476 (E1 ND) | k4 + neck split | 972 / 2127 / 1415 | 145/2647/1722 | 145 | 4.13 | 3.64 | 2 (1) | 0 |
| 6476 (E1 ND) | k4 + singles co-pos | 972 / 2126 / 1415 | 145/2647/1721 | 145 | 4.13 | 3.64 | 3 | 0 |
| HDL112 (E4 ND) | default (k4) | 528 / 1806 / 904 | 19/2397/822 | 19 | 5.68 | 1.84 | 0 | 0 |
| HDL112 (E4 ND) | k0 (no floor) | 528 / 1806 / 904 | 19/2397/822 | 19 | 5.68 | 1.84 | 0 | 0 |
| HDL112 (E4 ND) | k3 | 528 / 1806 / 904 | 19/2397/822 | 19 | 5.68 | 1.84 | 0 | 0 |
| HDL112 (E4 ND) | k4 + neck split | 528 / 1806 / 904 | 19/2397/822 | 19 | 5.68 | 1.84 | 0 (0) | 0 |
| HDL112 (E4 ND) | k4 + singles co-pos | 528 / 1806 / 904 | 19/2397/822 | 19 | 5.68 | 1.84 | 0 | 0 |
| 6450 (E3 Aab+) | default (k4) | 175 / 1343 / 1077 | 89/1568/938 | 89 | 0.51 | 1.85 | 3 | 0 |
| 6450 (E3 Aab+) | k0 (no floor) | 180 / 1379 / 1169 | 89/1704/935 | 89 | 1.13 | 1.90 | 3 | 0 |
| 6450 (E3 Aab+) | k3 | 180 / 1375 / 1157 | 89/1688/935 | 89 | 0.95 | 1.89 | 3 | 0 |
| 6450 (E3 Aab+) | k4 + neck split | 174 / 1347 / 1077 | 89/1568/941 | 89 | 0.51 | 1.85 | 0 (3) | 0 |
| 6450 (E3 Aab+) | k4 + singles co-pos | 175 / 1343 / 1042 | 89/1533/938 | 89 | 0.00 | 1.84 | 3 | 35 |
| 6533 (E4 T1D) | default (k4) | 301 / 1977 / 1916 | 117/3120/957 | 117 | 5.72 | 1.74 | 22 | 0 |
| 6533 (E4 T1D) | k0 (no floor) | 434 / 4928 / 5176 | 113/9477/948 | 113 | 38.26 | 2.74 | 22 | 0 |
| 6533 (E4 T1D) | k3 | 304 / 2515 / 2682 | 115/4432/954 | 115 | 13.47 | 1.89 | 21 | 0 |
| 6533 (E4 T1D) | k4 + neck split | 298 / 2007 / 1925 | 118/3120/992 | 118 | 5.72 | 1.74 | 6 (37) | 0 |
| 6533 (E4 T1D) | k4 + singles co-pos | 301 / 1977 / 1336 | 117/2540/957 | 117 | 0.00 | 1.72 | 22 | 580 |

| image | tissue mm2 | QuPath wall s (12 combos) | pixel trace s | geometry s per combo | max RSS GB (16 workers) |
|---|---|---|---|---|---|
| 6476 (E1 ND) | 196 | 90 | 45 | 0.4-1.6 | 97 |
| HDL112 (E4 ND) | 111 | 52 | 23 | 0.2-1.1 | 85 |
| 6450 (E3 Aab+) | 69 | 38 | 16 | 0.2-1.2 | 68 |
| 6533 (E4 T1D) | 101 | 58 | 29 | 0.4-2.2 | 84 |

**Why FLOOR_K = 4 (provisional).** The floor binds only where tau x anchor falls below the latent noise; it never
changes the low-clamp images (6476, HDL112: identical at k 0/3/4). Isolated INS-only cell-size blobs (20-80 um^2,
>= 10 px from other hormone signal; all 24 images, `/scratch/io60/work/isletC2/py/speck2.py`) in the three most clamped images vs the low-clamp
ND images of the cohort (1.4-7.8 /mm^2): k 3 leaves 6461 (E2 ND) 24.6, 6521 (E1 Aab+) 18.9 and 6533 (E4 T1D) 15.4;
k 4 gives 9.0, 9.4 and 6.5, and 6533 matches its same-panel ND partner HDL112 (detector: 5.72 vs 5.68 INS-only
singles/mm^2); k 5 gives 4.3, 5.3, 2.9. Over the 24 images the effective INS threshold is 3.0-7.4 latent sigma
above mu at k 0; k 4 raises the floor of that range to 4. The benchmark decides between 3, 4 and 5.

**Tau sweeps** (k 4): tau_CD99 2-4 changes the object count by < 2% on HDL112 (it moves objects between routes);
tau_H 4 -> 7 takes HDL112 singles from 1583 to 482 and Islets from 563 to 498: tau_H still dominates.

**Neck split** cut 37 necks on 6533 (22 -> 6 multi-core objects; Islets 301 -> 298), 3 on 6450, 1 on 6476, none
on HDL112. **Singles co-positivity** removes all INS-only singles where INS is heavily clamped (6533: 580 singles,
30%; 6450: 35) and nothing where it is not (6476, HDL112).

**Benchmark v2** (6450panc E3 Aab+, 6533 E4 T1D, HDL112 E4 ND; `/scratch/io60/benchmark_v2`):

| image | stratum | tiles per category | tiles | categories filled | empty pools | short pools | seconds |
|---|---|---|---|---|---|---|---|
| 6450panc | E3_Aab+ | 2 | 38 | 19/20 | pos_legacy_only | - | 173 |
| 6533 | E4_T1D | 2 | 38 | 19/20 | pos_legacy_only | - | 236 |
| HDL112 | E4_ND | 3 | 58 | 20/20 | - | neg_exo_cd99 (1/3, pool 1) | 316 |

- Candidates: 10 sources per image (v3 default k4, k0, k3, neck, copos, tau_CD99 2 and 4, tau_H 4 and 7, legacy
  Islets2) -> 3,746 / 11,045 / 4,241 groups. `hard_near_threshold` needs the tau variants: with only k/neck/copos
  sources HDL112 had no near-threshold candidate (the floor never binds there).
- Coverage (`private/coverage.tsv`): E4_ND at quota in 19/20 categories; `neg_exo_cd99` SHORT (1/3: HDL112 has one
  exocrine window with a CD99+ single candidate; its CD99 anchor is the cohort's highest, 864). E3_Aab+ and E4_T1D
  are PARTIAL at 2/3 until 6505 and 6580 are tiled; `pos_legacy_only` pools are empty there (every legacy islet is
  found by some v3 setting). 7 strata are NOT_RUN.
- Leak check: 134 tiles, 0 problems (names, donor IDs, pixel size 0.5078, palette, OME elements, mtimes, key).
- Full cohort at the defaults: 20 categories x 39 tiles (at least one per image and category) ~ 750 tiles of 400 um.

## 7. Known limitations

- **Clamping is not removed, only floored.** The half-max anchor still tracks the INS clamp fraction (Spearman
  -0.96); the floor keeps the threshold at >= k latent sigma. On the most clamped images INS positivity is therefore
  defined by the censored-normal model, on the others by the anchor.
- **Hormone noise still dominates the small objects.** Most singles and clusters are hormone-only (route 2); tau_H
  moves singles 3x between 4 and 7. Singles remain the least reliable class; the benchmark's `pos_single`,
  `neg_speckle` and `hard_near_threshold` categories target exactly this.
- **Hormone-positive bridges are not guarded** (only CD99-only necks can be cut, and only with `NECK_SPLIT`).
- **The occupied-tissue mask depends on DAPI** having an empty-space peak; on 9/24 images DAPI is unimodal and every
  tissue pixel counts as background.
- **6457**: its Tissue was imported into the `.qpdata` from `results/Geojson/Tissue/...` on 2026-09-28 and exported with
  the others, so it runs like every image. Do NOT use `--tissue-geojson` with the old file: the polygons are the same,
  but each ring starts at a different vertex, so its hash would not match script 1's (review 2026-09-29).
- **Images still to be annotated** (6579's pLN; 6347 and 6580) need `export_annotations.sh` + `islet_anchors.py
  --skip-existing` after drawing: the anchor check refuses script 1 until then.
- **pLN is disjoint from Tissue** in the four images that have one; the detector and bands exclude pLN either way.

## 8. Report only (not corrected)

- **Panel 2 (E3) SST speckle is aliased to disease.** Isolated SST-only cell-size blobs: E3 T1D median 4.0/mm^2
  (n = 5, 2.5-5.4), E3 Aab+ 7.7 (n = 2, 2.2-13.1), against 0.7 (0.2-1.8) in ND and 2.5 (0.7-7.9) in T1D of
  E1/E2/E4. E3 has no ND donor, so the panel-2 excess cannot be separated from disease.
- **CD99 anchors are lower in T1D.** Median CD99 anchor T1D/ND = 0.66 over the cohort (262 vs 398 counts; geometric
  mean 0.57); within group 0.49 (E1), 0.74 (E2), 0.51 (E4). With tau_CD99 relative to the anchor, the absolute CD99
  threshold is correspondingly lower in T1D images.

## 9. Decisions for you

1. **Annotation budget.** The default design (quota 3 per category and stratum, at least 1 per image) gives about
   750 tiles of 400 um for the cohort; `QUOTA_PER_STRATUM=2` about 550. Everything else has a provisional default
   (tau_CD99 3, tau_H 5, k 4, no neck split, no singles co-positivity) that the annotated benchmark will score
   against the sweep variants above.
