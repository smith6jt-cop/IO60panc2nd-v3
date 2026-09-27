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
