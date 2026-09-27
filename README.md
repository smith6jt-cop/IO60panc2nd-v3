# IO60panc2nd QuPath scripts

Scripts for the IO60panc2nd QuPath 0.7 project (PhenoCycler pancreas, nPOD donors).
They live inside the QuPath project folder because each stage loads
`scripts/lib/common.groovy` through `PROJECT_BASE_DIR`.

| Path | What |
|---|---|
| `0_…8_*.groovy`, `lib/` | Current pipeline stages (0–4 and 8 updated for the v3 re-run). |
| `v2/` | The stage versions that produced the v2 data. |
| `v1/`, `WorkingScripts2026April23/` | Earlier history. |
| `export_helpers/` | Export utilities (cells/islets GeoJSON, density maps). |
| `v3/` | The v3 re-run: storage/backup helpers (Phase 0), driver, islet method. |

The v3 plan (reuse, data management, normalization, islet redefinition, compute placement)
is summarised in `v3/README.md` once the driver lands. Nothing here touches data by itself;
every v3 helper logs to `results/v3/logs/` and writes manifests to `results/v3/manifests/`.

## v3 benchmark scoring and confound tests

`v3/islet_benchmark_score.py` and `v3/confound_tests.py` are pure-Python analysis
modules (see `v3/requirements-analysis.txt`) for the v3 endocrine-object
redefinition (Drotar et al. 2025): endocrine objects are classed by mean
caliper diameter (EndocrineSingle <= 10 um, EndocrineCluster 10-35 um,
Islet >= 35 um), recomputed from geometry for both detector output
(`pred/<code>.geojson`) and annotator ground truth (`truth/<code>.geojson`)
described by a benchmark `key.tsv` and `donors.tsv`. `islet_benchmark_score.py`
matches predictions to truth per tile by IoU, reports precision/recall/F1
per stratum (stain_group x disease_status x size_class x tile_kind) with a
donor-block bootstrap, and ranks detector parameter sweeps by macro-F1.
`confound_tests.py` checks whether donor-level outcomes are confounded by
stain batch rather than disease biology (a within-disease-stratum
permutation test, per-batch sensitivity of a disease contrast, and a
MixedLM with a donor random intercept). Run `pytest v3/tests` to exercise
both modules against synthetic fixtures (no project data files required).
