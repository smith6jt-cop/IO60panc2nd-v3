# scripts/v3: the v3 re-run execution layer

This directory holds the local driver for the v3 QuPath passes, the steps it runs, and the HiPerGator
templates. The islet method itself is kept separately in `qupath/*_v3.groovy`, `islet_anchors.py` and
`islet_benchmark_*`. The storage and backup helpers (`lib.sh`, `archive_v2.sh`, `export_annotations.*`,
`globus_*`, `verify_share_backup.sh`) come from Phase 0.

| File | Purpose |
|---|---|
| `run_v3.sh` | Resumable scheduler. It handles locks, lanes, the memory budget, the free-space and busy guards, and logs. |
| `v3ctl.py` | Every decision that involves JSON: plan, instantiate a step, validate its sentinel, check drift, convert to Parquet, write the run log, write the report. |
| `steps.tsv` | The step table: step → script(s), passes, save policy, lane, heap, threads, timeout, dependencies, overrides, postcondition, post hook. |
| `footer_template.groovy` | Code appended to every step script. It checks the postcondition, saves when the policy allows, and writes the provenance/sentinel. |
| `preflight.groovy` | Gate for passes I and C. It checks channels, pixel size 0.5078, exactly one Tissue, and that every pLN is at least 99% inside Tissue, and it removes stale `Other: Gal3Norm`. |
| `exports.groovy` | Read-only exports in four modes: `cells` (full measurement TSV), `context`, `annotations` and `islets` (GeoJSON). |
| `qupath/export_cells_geojson_v3.groovy` | Fixed cell-polygon exporter: a correct duplicate check, and the output is stream-verified against every detection UUID. |
| `import_reviewed_islets.groovy` | Pass R (HiPerGator): replaces the islet-type annotations with the reviewed GeoJSON. |
| `hpg/` | HiPerGator bundle (see the HiPerGator section below). It cannot be run on this workstation. |
| `tests/` | Fixture builder, unit tests (`test_v3ctl.py`), end-to-end suite (`test_run_v3_fixture.sh`) and deliberately broken steps. |

## Passes and steps

```
I  preflight -> 1 islets -> 4 bands                          (light lane)
C  preflight -> 0 InstanSeg (gpu lane) -> 2 measurements (cpu) -> C_cells (TSV -> Parquet), C_geojson
J  3 signed distance -> 8 parenting -> J_context (TSV -> Parquet), J_annotations (TSV), J_islets (GeoJSON)
X  R_export: reviewed islets -> results/v3/reviewed_islets/<tag>.geojson (after gate 1, for HiPerGator)
R  R_import: HiPerGator only; runs before pass J there
```

Scripts are chosen through `steps.tsv`:
- `a|b` lists candidates. `--script-set auto` (the default) uses the v3 copy when it exists and falls back to the original otherwise.
- `orig` or `v3` forces one set.
- `--script STEP=PATH` substitutes a single script, for testing.
- The plan prints which file each step will use. Check it before a real run: once the islet-method files appear, `auto` picks them up.

Overrides such as `N_THREADS=24` for step 2 rewrite the declaration in the script copy. Each override must match exactly once, or the step is refused.

## Usage

```bash
run_v3.sh --pass I --dry-run                                     # the plan for all 24 images; changes nothing
run_v3.sh --pass I,C --images 6374,6533 --lanes-cpu 1            # real run on two images
run_v3.sh --pass J --images 6374 --force 8                       # re-run 8 and everything downstream of it
run_v3.sh --accept-external --images 6374 --note "gate-1 edits"  # after editing an image in the GUI
run_v3.sh --project /scratch/io60/fixture/project/project.qpproj --pass I,C,J --script-set orig   # fixture
```

- **Lanes:** `--lanes-gpu/-cpu/-light` (default 1 each) are pools of flock slots shared by every `run_v3.sh` on the host.
- **Memory budget:** `--mem-budget` (default 400 GB) caps the sum of `-Xmx` of the steps running at once. It uses `-Xmx` rather than actual use, because the JVM grows toward `-Xmx`: step 2 reached 69 GB resident on the 2048² fixture.
- **Timeouts and thread counts:** per step, from the table. Threads are set through `ForkJoinPool.common.parallelism`, and step 2 has its own `N_THREADS`.
- **Defaults are starting points:** lane counts, heaps and timeouts are for the pilot to benchmark (1, 2 and 4 lanes).

## How a step runs

1. **Plan under the image lock** (`/scratch/io60/locks/<tag>.lock`). The lock is held for the whole chain of that image. QuPath inherits the lock descriptor, so even an orphaned QuPath keeps its image locked.
2. **Guards:**
   - Free space: the data filesystem needs at least 2 × data.qpdata + 50 GB, because QuPath keeps the old file as `.bkp` while it writes, and `/` needs at least 300 GB.
   - A QuPath GUI appearing mid-run stops new steps.
3. **Slot and memory:** wait for a lane slot and a memory reservation.
4. **`v3ctl.py prepare`:**
   - It creates a private worker project under `/scratch/io60/work/driver/<run>/<tag>/<step>.<ns>/`: a copy of `project.qpproj`, with `data/ classifiers/ scripts/ resources/` symlinked.
   - It copies the step script, applies the overrides and appends the footer with a baked-in configuration (attempt id, sequence number, save policy, postcondition, qpdata stat before).
   - It moves the step's previous records to `history/`.
5. **Run QuPath:** `JAVA_TOOL_OPTIONS=-Xmx<heap>g … timeout … QuPath script -p <worker project> -i <image> <script>`. Nothing uses `-s`.
6. **The footer** runs only if the body finished. It:
   - checks the postcondition;
   - for `footer-saves`, calls `getProjectEntry().saveImageData(...)` only when the check passed (and not when the step reported `changed: false`);
   - for `self-saves` (script 0), requires that data.qpdata was written during this JVM;
   - for `never`, requires that it was not written;
   - writes `provenance/<tag>/<step>.json` atomically.
7. **`v3ctl.py finish`:**
   - The step passes only if the sentinel exists, belongs to this attempt, says `pass`, the script hash matches, and data.qpdata still matches what the footer recorded.
   - It copies the script to `provenance/<tag>/<step>.groovy`, checks for drift and appends a row to `run_log.csv` (under flock).
8. **`v3ctl.py post`** (C_cells, J_context):
   - DuckDB reads the TSV with explicit types: text columns as VARCHAR, everything else as DOUBLE (`NaN` → NaN, empty → NULL).
   - It writes zstd Parquet and verifies that the Parquet rows, the TSV rows and the detection count all agree, and that Object ID is non-null and unique.
   - For J_context it also checks that the context has the same Object IDs as the cells Parquet.
   - It then deletes the TSV, except for `0080870_6533` and `0081366_HDL112`, whose TSV moves to `cells_tsv/`.

**Why a sentinel and not the exit code.** It was checked on QuPath 0.7 / Groovy 5:
- A top-level `return` exits 0.
- A throw exits 1.
- A trailing `return` followed by the footer compiles, but the footer never runs.
- A body that ends in method or class definitions still runs the footer.

In every failing case the sentinel is missing, and the driver records FAIL and stops that image's chain.

**Footer rule.** No `import` lines, fully-qualified names, and every identifier, closure parameters included, starts with `__v3_`. Groovy rejects re-declaring a name the step already declared.

## Resumability, custody, drift

Every record carries a per-image sequence number, `seq`.

**When a step is done.** A step counts as done when all of the following hold, and it is skipped:
- its `<step>.json` says pass;
- its post hook passed for the same attempt;
- its recorded outputs still exist with the same size;
- every `after` dependency is done and older (a re-run upstream step makes it stale);
- every `needs` dependency is done.

**What re-runs.** A step that re-runs (failed, missing or `--force`d) takes everything downstream of it in the selected passes along with it. A step whose prerequisite is outside the selected passes and not done is **BLOCKED**, and so is the rest of that image's chain.

**Custody.** Before anything runs, data.qpdata must match the `after` stat of the image's highest-seq record: size exactly, mtime within 1 s so that Globus copies still match. If it does not match, the file was changed outside run_v3 (GUI review, a restore, an interrupted save), and the image is **blocked** rather than re-run. Re-running could overwrite review edits. `--accept-external --note …` records the current file and unblocks it.

**Drift.** On reload, QuPath 0.7 takes an image's server metadata (channel names, pixel size) from the `project.qpproj` entry. `saveImageData()` then rewrites only `data/<id>/server.json`, and headless runs without `-s` never sync `project.qpproj`, so a metadata change made by a step silently reverts. After each step the driver compares both the worker project entry and a freshly written `server.json` with the project entry. On any difference it writes `<step>.drift.json` and blocks the image. Nothing is merged. Fix the entry in the GUI, delete the file, and re-run (or pass `--allow-drift`). The run report also lists worker-versus-project entry differences.

**Busy check.** The driver refuses to start while a QuPath process (GUI or headless) or a Senior rebuild runs. It uses `assert_project_idle` from `lib.sh`, narrowed to real processes by name, because that helper also matches shell wrappers that merely mention QuPath. `--allow-busy` overrides it and forces a single CPU lane. Ctrl-C/TERM lets running steps finish and starts no new ones; QuPath runs in `timeout`'s own process group.

## Output layout (under `--out`; default `results/v3` for the master project, `<project dir>/results/v3` otherwise)

```
provenance/<tag>/<step>.json (+ .post.json, .groovy, .drift.json), external.<seq>.json, history/
cells/<tag>.parquet   cells_tsv/<tag>.tsv (6533, HDL112)   cell_geojson/<tag>.geojson
context/<tag>.parquet   annotations/<tag>.tsv   islets/<tag>.geojson   reviewed_islets/<tag>.geojson
logs/run_v3/<tag>/<step>.log (appended per attempt)   logs/run_v3/<run>.driver.log   run_log.csv
```

Worker projects and time files stay in `/scratch/io60/work/driver/<run>/`. TSVs are written to `/scratch/io60/tmp/run_v3/` and deleted after conversion.

## Safety rules

- The driver never opens the target `project.qpproj` in QuPath; every step gets a private copy. It never uses `-s`, and nothing auto-merges into the project.
- Test against the fixture (`tests/`), never the master project. On the master, use only `--dry-run` until you run the pilot deliberately.
- `/home` has little headroom. After step 2 the .qpdata files grow several-fold, and the guard stops the run before a save could fill the disk.
- Before any real run:
  - check which scripts the plan prints (the islet-method v3 files are picked up automatically by `auto`);
  - fix 6450panc's pixel size, which is 0.5077663810243286 µm in its project entry, so preflight fails it by design.

## Tests

```bash
/scratch/io60/venv-v3/bin/python -m pytest -q scripts/v3/tests/test_v3ctl.py     # 19 unit tests, no QuPath
scripts/v3/tests/make_fixture.sh --force                                         # 2048² crop of 6450panc on /scratch
TEST_ALLOW_BUSY=1 scripts/v3/tests/test_run_v3_fixture.sh                         # 25 end-to-end checks, ~3 min
```

## HiPerGator (`hpg/`; templates, not runnable here)

The bundle relies on identical absolute paths: `/home/smith6jt/IO60panc2nd` holds the copied `project.qpproj`, `classifiers/` and `scripts/`, and `data`, `results` and `Images/<file>` are symlinks into `/blue/maigan/smith6jt/IO60panc2nd_v3`. With that layout, `run_v3.sh`, `steps.tsv` and the footer run unchanged there.

| File | Purpose |
|---|---|
| `setup_hpg.sh` | Checks and repairs that layout, the portable QuPath, the InstanSeg model, `~/.djl.ai` and the venv. |
| `c0_gpu.sbatch` | Step 0 only, one image per task. It first runs the DJL/CUDA fast-fail check (`djl_cuda_check.groovy`, tested locally), logs `nvidia-smi dmon`, can stage to `$TMPDIR` with a per-image symlink repoint restored on exit, then runs `run_v3.sh --steps 0`. |
| `c2_cpu.sbatch` | Step 2 array. |
| `j_cpu.sbatch` | `--pass R,J`; it forces `R_import` when the reviewed GeoJSON's sha256 differs from the imported one. |
| `hpg_submit.sh` | Prints `showQos`, `slurmInfo` and `squeue -A` for both accounts. Capacity comes from `sacctmgr` GrpTRES minus running `tres-alloc`; it fills clive, then maigan (at most half, at most 1 GPU), then maigan-b (CPU only), and splits the image list into one array per account. Every decision is logged. |

**Provenance must travel with the data.** Sync `results/v3/provenance/` both ways with `data/`, or step 0 on HiPerGator is blocked because the preflight record is missing.

**Known limitations / TODOs:**
- The `showQos`/`slurmInfo` parsing, the partition names per GPU type and the Python module name need a sample of the real output (marked `TODO(HiPerGator)`).
- shellcheck was not available on this workstation; the scripts pass only `bash -n`.
- Pass J does not check a gate-1 sign-off. The review file format is not fixed yet.
- Input-content changes are not tracked in general: a step re-runs when a step it depends on re-runs, not when an input file changes. The exception is the reviewed-islet hash in `j_cpu.sbatch`.
