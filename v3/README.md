# scripts/v3: the v3 re-run execution layer

This directory holds the local driver for the v3 QuPath passes, the steps it runs, and the HiPerGator
templates. The islet method itself is kept separately in `qupath/{1,3,4,8}_*_v3.groovy`, `islet_anchors.py` and
`islet_benchmark_*`. The storage and backup helpers (`lib.sh`, `archive_v2.sh`, `export_annotations.*`,
`globus_*`, `verify_share_backup.sh`) come from Phase 0.

| File | Purpose |
|---|---|
| `run_v3.sh` | Resumable scheduler. It handles locks, lanes, the memory budget, deadlines, the free-space and busy guards, and logs. |
| `v3ctl.py` | Every decision that involves JSON: plan, resolve images, instantiate a step, validate its sentinel, check drift, convert to Parquet, merge provenance across sites, write the run log and the report. |
| `steps.tsv` | The step table: step → script(s), passes, save policy, lane, heap, threads, timeout, dependencies, what an external edit does to it, overrides, postcondition, post hook. |
| `footer_template.groovy` | Code appended to every step script. It checks the postcondition, re-checks that data.qpdata is the file the step started from, saves when the policy allows and the deadline is far enough, and writes the provenance/sentinel. |
| `preflight.groovy` | Gate for passes I and C. It checks channels, pixel size 0.5078, exactly one Tissue, and that every pLN is disjoint from Tissue (≤ 1% inside) or contained in it (≥ 99%); a partial overlap fails. It removes stale `Other: Gal3Norm`. |
| `qupath/0_cell_detection_panc_v3.groovy` | Step 0: InstanSeg over every **Tissue and pLN** region (v3 copy of `scripts/0_cell_detection_panc.groovy`). |
| `qupath/7_tiled_pixel_classifiers_v3.groovy` | Script 7 with the tile area = union of Tissue and pLN (not a run_v3 step; run it from the GUI as before). |
| `exports.groovy` | Read-only exports in four modes: `cells` (measurement TSV without signed distances), `context`, `annotations` and `islets` (GeoJSON). |
| `qupath/export_cells_geojson_v3.groovy` | Fixed cell-polygon exporter: a correct duplicate check, and the output is stream-verified against every detection UUID. |
| `import_reviewed_islets.groovy` | Pass R (HiPerGator): replaces the islet-type annotations with the reviewed GeoJSON. |
| `hpg/` | HiPerGator bundle (see the HiPerGator section below). It cannot be run on this workstation. |
| `tests/` | Fixture builder (4 variants), unit tests (`test_v3ctl.py`), end-to-end suite (`test_run_v3_fixture.sh`), mocked-SLURM test of the HiPerGator bundle (`test_hpg_submit.sh`) and deliberately broken steps. |

## Passes and steps

```
I  preflight -> 1 islets -> 4 bands                          (light lane)
C  preflight -> 0 InstanSeg over Tissue + pLN (gpu lane) -> 2 measurements (cpu) -> C_cells (TSV -> Parquet), C_geojson
J  3 signed distance -> 8 parenting -> J_context (TSV -> Parquet), J_annotations (TSV), J_islets (GeoJSON)
X  R_export: reviewed islets -> results/v3/reviewed_islets/<tag>.geojson (after gate 1, for HiPerGator)
R  R_import: HiPerGator only; runs before pass J there
```

Scripts are chosen through `steps.tsv`:
- `a|b` lists candidates. `--script-set v3` (**the default**) requires the v3 copy and fails the plan if it is missing; it never falls back.
- `orig` forces the originals (the fixture tests use it for 1/3/4/8); `auto` takes the v3 copy when it exists.
- Step 0 has one candidate, the v3 copy, whatever the set.
- `--script STEP=PATH` substitutes a single script, for testing.
- The plan prints which file each step will use.

Overrides such as `N_THREADS=24` for step 2 rewrite the declaration in the script copy. Each override must match exactly once, or the step is refused.

### Lymph nodes (pLN) and step 0

Every script that selects a region handles both `Tissue` and `pLN` annotations:
- A pLN **disjoint** from Tissue (6533, HDL112, HDL115, 6551) is segmented as a region of its own.
  - Disjoint means sharing at most 10 µm² with Tissue.
  - Its detections become its children, and step 8 parents them to pLN.
- A pLN **contained** in Tissue (≥ 99% of its area) is not selected: the Tissue pass already segments it.
- A **partial** overlap is refused, because it would segment the shared area twice.
- Two selected regions sharing more than 10 µm² are refused and reported.
- The postcondition `regions_detected(classes=Tissue|pLN)` requires every Tissue and every pLN annotation to hold detection centroids.
- Script 0 keeps the original's other rules:
  - explicit selection, never `selectAnnotations()`;
  - `removeDetections()` first;
  - children of every region parked and restored;
  - the same InstanSeg parameters and 13 channels;
  - the result is verified, the duplicate check runs, and it throws rather than save on failure.

## Usage

```bash
run_v3.sh --pass I --dry-run                                        # the plan for all 24 images; changes nothing
run_v3.sh --pass I,C --images id:1,id:37 --lanes-cpu 1              # real run on two images (exact: id:N, name or tag)
run_v3.sh --pass J --images 6374_Scan1.er.qptiff_-_resolution_1 --force 8   # re-run 8 and everything downstream of it
run_v3.sh --accept-external --images id:1 --note "gate-1 edits"     # after editing an image in the GUI
run_v3.sh --accept-external --images id:1 --note "deleted cells" --rerun-from 2
run_v3.sh --merge-provenance --images-file pilot.txt                # adopt the records HiPerGator wrote
run_v3.sh --project /scratch/io60/fixture/project/project.qpproj --pass I,C,J --script-set orig   # fixture
```

- **Images:** `--images` takes a comma list of exact image names, tags, or `id:<entry id>`. `--images-file` takes one per line; blank lines and `#` comments are ignored.
  - Every token must match exactly one image, and no image may be named twice.
  - Substrings are refused, and the error lists the images whose name contains the token.
  - The resolved list is printed and written to the driver log.
- **Force:** `--force` and `--rerun-from` must name steps inside the selected passes or steps; otherwise the plan is refused.
- **Lanes:** `--lanes-gpu/-cpu/-light` (default 1 each) are pools of flock slots shared by every `run_v3.sh` on the host.
- **Memory budget:** `--mem-budget` (default **300 GB**) caps the sum of the reservations of running steps.
  - A reservation is the step's `-Xmx`, because the JVM grows toward it (step 2 reached 69 GB resident on the 2048² fixture).
  - A step with a Parquet post hook reserves at least the DuckDB limit (`--duckdb-mem-gb`, default 48) + 8 GB.
- **Deadlines:** no save starts later than `--save-margin-min` (default 10) before the step's deadline.
  - The deadline is its timeout or, under SLURM, the job's end (`V3_DEADLINE_EPOCH`), whichever comes first.
  - A step killed mid-save leaves data.qpdata half-written.
- **Timeouts and thread counts:** per step, from the table. Threads are set through `ForkJoinPool.common.parallelism`, and step 2 has its own `N_THREADS`.
- **Defaults are starting points:** lane counts, heaps and timeouts are for the pilot to benchmark (1, 2 and 4 lanes).

## How a step runs

1. **Plan under the image lock** (`/scratch/io60/locks/<tag>.lock`, plus a mkdir lock in `V3_SHARED_LOCKDIR` on a cluster).
   - The lock is held for the whole chain of that image, and QuPath inherits the flock descriptor.
   - A re-plan that fails, or that does not describe exactly this image, is a PLAN ERROR, never "done".
2. **Guards:**
   - Free space: the data filesystem needs at least 2 × data.qpdata + 50 GB, and `/` needs at least 300 GB.
   - A QuPath GUI appearing mid-run stops new steps.
   - A SLURM deadline closer than the save margin + 5 min stops new steps.
3. **Slot and memory:** wait for a lane slot and a memory reservation.
4. **`v3ctl.py prepare`:**
   - It re-checks custody; the file may have changed since the plan.
   - It creates a private worker project under `/scratch/io60/work/driver/<run>/<tag>/<step>.<ns>/`.
   - It copies the step script, applies the overrides and appends the footer with a baked-in configuration. The configuration holds the attempt id, sequence number, save policy, postcondition, qpdata stat at the start, deadline and save margin, and the scripts repo commit and dirty flag.
   - It writes the same configuration as `step_config.json` (`{CFG}`, read by step 0).
   - It moves the step's previous records to `history/`.
   - It prints `KEY=VALUE` lines that run_v3.sh **parses** (never `eval`s), and a NUL-separated argument file.
5. **Run QuPath:** `JAVA_TOOL_OPTIONS=-Xmx<heap>g … timeout … QuPath script -p <worker project> -i <image> <script>`. Nothing uses `-s`.
6. **The footer** runs only if the body finished. It:
   - checks the postcondition;
   - for `footer-saves`, **re-stats data.qpdata** and refuses to save (throws, **no sentinel**) if it differs from the stat taken when the step began;
   - skips the save if the deadline is within the margin, which is a failing sentinel;
   - otherwise saves when the check passed and the step reported a change. A save that throws means no sentinel;
   - refuses (no sentinel) if data.qpdata changed during a `never` step, or during a `footer-saves` step that did not save;
   - for `self-saves` (step 0, which guards its own save the same way), requires that data.qpdata was written during this JVM;
   - writes `provenance/<tag>/<site>/<step>.json` atomically. The file records the site, the scripts git commit, dirty flag and whether the step script itself is committed, and the `.bkp` stat.
7. **`v3ctl.py finish`:**
   - The step passes only if the sentinel exists, belongs to this attempt, says `pass`, the script hash matches, and data.qpdata still matches what the footer recorded.
   - With no sentinel, it quotes the footer's or step 0's refusal from the log.
   - It names a leftover `data.qpdata.bkp`. QuPath renames the old file to `.bkp` before writing and deletes it after a complete write; `readImageData()` silently loads the `.bkp` if data.qpdata is unreadable.
   - It writes `<step>.finish.json` (its verdict), copies the script, checks for drift and appends a row to `run_log.csv`.
8. **`v3ctl.py post`** (C_cells, J_context):
   - DuckDB reads the TSV with explicit types: text columns as VARCHAR, everything else as DOUBLE.
   - It writes zstd Parquet and verifies that the Parquet rows, the TSV rows and the detection count all agree, and that Object ID is non-null and unique.
   - For C_cells, `Cell: DAPI: Mean` must be non-null and non-NaN on every row.
   - For J_context, the context must have the same Object IDs as the cells Parquet.
   - It then deletes the TSV, except for `0080870_6533` and `0081366_HDL112`, whose TSV moves to `cells_tsv/`.
   - The cells export leaves out `Signed distance…` columns (they belong to the context export), so its column set is the same whether it runs before or after step 3.

**Why a sentinel and not the exit code.** QuPath 0.7 / Groovy 5 behaves as follows:
- A top-level `return` exits 0.
- A throw exits 1.
- A trailing `return` followed by the footer compiles, but the footer never runs.
- A body that ends in method or class definitions still runs the footer.

In every failing case the sentinel is missing, and the driver records FAIL and stops that image's chain.

**Footer rule.** No `import` lines, fully-qualified names, and every identifier, closure parameters included, starts with `__v3_`.

## Resumability, custody, sites, drift

Every record carries a per-image sequence number `seq` and the data.qpdata stat before and after it.

**When a step is done.** A step counts as done, and is skipped, when all of the following hold:
- its sentinel says pass **and** its `finish` record passed for the same attempt. An orphaned QuPath's late sentinel does not count;
- its post hook passed for the same attempt, and its recorded outputs still exist with the same size;
- its **script is still the one the plan resolves**: same source path and sha256. A different script set, an edited script or a `--script` override makes it stale;
- **no external edit was accepted after it**, if its `on_external` is `stale` (preflight, C_cells, C_geojson, 3, 8, J_*, R_export: they read the edited image);
- it is not named by an external edit's `--rerun-from` (or downstream of one);
- every `after` dependency is done and older.

`needs` dependencies only gate whether a step may **start**. A stale preflight never makes a reviewed step 1 re-run.

**What re-runs.** A step that re-runs (failed, missing, stale or `--force`d) takes everything downstream of it in the selected passes along with it. A step whose prerequisite is outside the selected passes and not done is **BLOCKED**, and so is the rest of that image's chain. Step 3 and J_context need step 2, so pass J cannot start while step 2 is unfinished (e.g. in another job).

**Custody.** Before anything runs, the records, sorted by seq, must form **one chain**: each record starts from the file the previous one left, and a failed record did not change the file. The chain must also end at the current data.qpdata: size exactly, mtime within 1 s so that Globus copies still match.
- If the file changed outside run_v3 (GUI review, a restore, an interrupted save), the image is **blocked** rather than re-run, because re-running could overwrite review edits.
- `--accept-external --note … [--rerun-from STEP]` records the current file under the image lock. Afterwards every `on_external=stale` step, and STEP with its downstream, re-run. At gate 1: accept, then `--pass X` re-exports the reviewed islets.
- Two records with the same seq (the image was processed on two sites from the same state), or a broken chain, block the image as `custody_conflict`. Resolve it by hand: keep one side's file and records.

**Sites.** Each machine writes only `provenance/<tag>/<site>/`. The site is `V3_SITE`: the workstation defaults to its short host name, and every HiPerGator node uses `hpg`.
- Syncing `results/v3/provenance/` both ways therefore never overwrites a record.
- Another site's records count only after `run_v3.sh --merge-provenance --images …` validated the combined chain against the current data.qpdata and mirrored them into `<site>/merged/<other site>/`. Until then the image is blocked as `unmerged`.
- An in-flight marker of another site blocks it too, unless that attempt wrote its finish record there.
- Records in the old flat layout are refused (`legacy_layout`).

**Drift.** On reload, QuPath 0.7 takes an image's server metadata (channel names, pixel size) from the `project.qpproj` entry.
- `saveImageData()` rewrites only `data/<id>/server.json`, and headless runs without `-s` never sync `project.qpproj`, so a metadata change made by a step silently reverts.
- After each step the driver compares both the worker project entry and a freshly written `server.json` with the project entry.
- On any difference it writes `<step>.drift.json` and blocks the image. Nothing is merged.
- Fix the entry in the GUI, delete the file, and re-run (or pass `--allow-drift`).

**Busy check.** The driver refuses to start while a QuPath process (GUI or headless) or a Senior rebuild runs.
- It uses `assert_project_idle` from `lib.sh`, narrowed to real processes by name.
- `--allow-busy` overrides it and forces a single CPU lane.
- Ctrl-C/TERM lets running steps finish and starts no new ones.

## Output layout (under `--out`; default `results/v3` for the master project, `<project dir>/results/v3` otherwise)

```
provenance/<tag>/<site>/<step>.json (+ .finish.json, .post.json, .groovy, .drift.json), external.<seq>.json,
                        history/, merged/<other site>/
cells/<tag>.parquet   cells_tsv/<tag>.tsv (6533, HDL112)   cell_geojson/<tag>.geojson
context/<tag>.parquet   annotations/<tag>.tsv   islets/<tag>.geojson   reviewed_islets/<tag>.geojson
logs/run_v3/<tag>/<step>.log (appended per attempt)   logs/run_v3/<run>.driver.log   run_log.csv
```

Worker projects and time files stay in `/scratch/io60/work/driver/<run>/`. TSVs are written to `/scratch/io60/tmp/run_v3/` and deleted after conversion.

## Safety rules

- The driver never opens the target `project.qpproj` in QuPath; every step gets a private copy. It never uses `-s`, and nothing auto-merges into the project.
- Test against the fixtures (`tests/`), never the master project. On the master, use only `--dry-run` until you run the pilot deliberately.
- `/home` has little headroom. After step 2 the .qpdata files grow several-fold, and the guard stops the run before a save could fill the disk.
- Before any real run:
  - check which scripts the plan prints;
  - fix 6450panc's pixel size, which is 0.5077663810243286 µm in its project entry, so preflight fails it by design.

## Tests

```bash
/scratch/io60/venv-v3/bin/python -m pytest -q scripts/v3/tests/test_v3ctl.py     # 52 unit tests, no QuPath
scripts/v3/tests/test_hpg_submit.sh                                               # 19 checks, mocked SLURM
TEST_ALLOW_BUSY=1 scripts/v3/tests/test_run_v3_fixture.sh                         # end-to-end on 4 fixtures, ~10 min
```

`make_fixture.sh [--force] [--variant pln|pln-partial|pln-sliver]` builds the fixtures on /scratch from a 2048² crop of 6450panc:
- the default variant has a pLN inside Tissue;
- `pln` has a disjoint pLN and a contained one;
- `pln-partial` fails preflight;
- `pln-sliver` passes preflight with a warning and step 0 refuses it.

## HiPerGator (`hpg/`; templates, not runnable here)

The bundle relies on identical absolute paths:
- `/home/smith6jt/IO60panc2nd` holds the copied `project.qpproj`, `classifiers/` and `scripts/`.
- `data`, `results` and `Images/<file>` are symlinks into `/blue/maigan/smith6jt/IO60panc2nd_v3`.
- With that layout, `run_v3.sh`, `steps.tsv` and the footer run unchanged there.

`hpg_env.sh` sets:
- `V3_SITE=hpg`;
- `V3_SHARED_LOCKDIR=$BLUE/locks`. flock does not work across nodes, so image locks also use mkdir there; a lock whose owner job has left squeue is broken.

| File | Purpose |
|---|---|
| `setup_hpg.sh` | Checks and repairs that layout, the portable QuPath, the InstanSeg model, `~/.djl.ai` and the venv. |
| `c0_gpu.sbatch` | Step 0 only, one image per task. It first checks that the job's time left covers the step timeout + margin, then runs the DJL/CUDA fast-fail check (`djl_cuda_check.groovy`), logs `nvidia-smi dmon`, can stage to `$TMPDIR`, then runs `run_v3.sh --steps 0`. |
| `c2_cpu.sbatch` | Step 2 array (time budget checked the same way). |
| `j_cpu.sbatch` | `--pass R,J`; it forces `R_import` when the reviewed GeoJSON's sha256 differs from the imported one (`v3ctl.py show-record`). |
| `hpg_submit.sh` | Normalizes the images file (exact names, blank lines dropped), refuses images already in one of your queued/running v3 jobs (registry `$BLUE/work/submit/jobs/`), and refuses images whose plan would not run (unmerged, blocked, nothing to do); `--drop-unready` submits the rest. Capacity comes from `sacctmgr` GrpTRES minus the group's RUNNING **and PENDING** jobs. It fills clive, then maigan (at most half, at most 1 GPU), then maigan-b (CPU only, at most `--maigan-b-max` 2 tasks), and sets `--time` = the sum of the job's step timeouts + margin. |

**Provenance travels with the data, but is adopted explicitly.** The round trip is:
1. Sync `data/<id>/` and `results/v3/provenance/` to HiPerGator.
2. There, run `run_v3.sh --merge-provenance --images-file …`, then `hpg_submit.sh`.
3. Sync both back.
4. On the workstation, run `--merge-provenance` again.

A merge refuses when the data.qpdata on disk is not the file the other site's last record left.

**Known limitations / TODOs:**
- The `showQos`/`slurmInfo` parsing, the partition names per GPU type and the Python module name need a sample of the real output (marked `TODO(HiPerGator)`).
- shellcheck was not available on this workstation; the scripts pass `bash -n` and the mocked test.
- Pass J does not check a gate-1 sign-off. The review file format is not fixed yet.
- The scripts copied to HiPerGator are not a git checkout unless you sync `.git` too, so records made there hold `repo_git.commit = null`. The source sha256 is always recorded.
