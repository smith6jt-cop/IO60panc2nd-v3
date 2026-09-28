#!/usr/bin/env bash
# test_run_v3_fixture.sh -- end-to-end tests of run_v3.sh on the /scratch fixture projects.
#
# Never touches the real project: every run passes --project /scratch/io60/fixture/<project>/project.qpproj
# (outputs land in <project>/results/v3). Uses the GPU for step 0 (2048x2048 px). Provenance site: ws.
#   1  full run of passes I, C, J from a freshly built fixture        -> every step PASS, image DONE; step 0 is
#                                                                         the v3 copy (Tissue + contained pLN)
#   2  second run                                                      -> nothing runs (resumability)
#   3  a step body ending in method/class definitions                  -> footer still runs, PASS
#   4  a step that leaves with a top-level `return`                    -> no sentinel, FAILED, chain stops
#   5  a step that throws                                              -> no sentinel, FAILED, qpdata untouched
#   6  a step that fails its postcondition                             -> footer does not save, FAILED
#   7  plain re-run                                                    -> failed + downstream steps re-run, DONE
#   8  data.qpdata changed outside run_v3 (touch)                      -> BLOCKED; --accept-external unblocks and
#                                                                         re-runs exactly the steps that read the
#                                                                         image (on_external=stale); the cells export
#                                                                         keeps its column set after step 3 ran
#   9  image lock held by someone else                                 -> LOCKED, nothing runs
#  10  free-space guard not met                                        -> STOP (exit 3), nothing runs
#  11  a step that changes the server metadata (on a copy)             -> DRIFT recorded, image blocked
#  12  --images substring / --force outside the selection              -> planning refused (exit 2)
#  13  data.qpdata written by something else during a step (copy)      -> footer refuses to save, NO sentinel,
#                                                                         next plan BLOCKED; accept-external
#                                                                         --rerun-from 0 records the invalidation
#  14  a save within --save-margin-min of the deadline (copy)          -> skipped, FAILED, qpdata untouched
#  15  v3 script 0's own guard: config says the file was different     -> REFUSED at start, qpdata untouched
#  16  DISJOINT pLN fixture (+ a contained pLN): passes I, C, J        -> pLN segmented as its own region, its
#                                                                         cells parented to pLN (v3 step 8)
#  17  pLN partly inside Tissue                                        -> preflight FAILS (partial)
#  18  pLN 5 px into Tissue (<1% of it, ~1700 um^2 shared)            -> preflight FAILS (partial, by step 0's
#                                                                         10 um^2 rule), before step 0 runs
#  19  two provenance sites on one image (copy)                        -> blocked until --merge-provenance,
#                                                                         then each site continues the chain
# Usage: test_run_v3_fixture.sh [--keep-fixture]   (default rebuilds the fixtures first)
set -uo pipefail
V3DIR=$(cd "$(dirname "$0")/.." && pwd)
R=$V3DIR/run_v3.sh
S=$V3DIR/tests/steps
PY=/scratch/io60/venv-v3/bin/python
QUPATH=/home/smith6jt/QuPath/bin/QuPath
FROOT=/scratch/io60/fixture
FIXDIR=$FROOT/project
FIX=$FIXDIR/project.qpproj
OUT=$FIXDIR/results/v3
TAG=fixture_6450.ome.tif
QP=$FIXDIR/data/1/data.qpdata
WORK=/scratch/io60/work/driver
LOGS=/scratch/io60/work/driver/tests/$(date +%Y%m%dT%H%M%S)
V3_8=/home/smith6jt/IO60panc2nd/scripts/v3/qupath/8_hierarchy_parenting_v3.groovy
export V3_SITE=ws
unset V3_SHARED_LOCKDIR V3_DEADLINE_EPOCH
mkdir -p "$LOGS"
rm -rf -- "$WORK"/t-*          # run dirs of earlier suite runs (the ids below are fixed)
FAILS=0
# TEST_ALLOW_BUSY=1: other sessions' headless QuPath runs (on their own private projects) would make
# run_v3.sh refuse to start; the fixture is private, so the tests may pass --allow-busy.
EXTRA=()
[[ ${TEST_ALLOW_BUSY:-0} == 1 ]] && EXTRA=(--allow-busy)

check() { local d=$1; shift; if "$@"; then printf '  PASS  %s\n' "$d"; else printf '  FAIL  %s\n' "$d"; FAILS=$((FAILS + 1)); fi; }
show() { grep -E ' (start|PASS|FAIL|BLOCKED|LOCKED|STOP|DRIFT|post )| DONE| FAILED|refus|images:|matches|not among' "$1" | grep -v '^== ' | sed 's/^/      /' | head -40; }
runp() {   # runp <project file> <id> <args...> ; sets RC
    local proj=$1 id=$2; shift 2
    "$R" --project "$proj" --script-set orig --run-id "$id" "${EXTRA[@]}" "$@" > "$LOGS/$id.out" 2>&1
    RC=$?
    printf '\n[%s] exit %s: run_v3.sh %s\n' "$id" "$RC" "$*"
    show "$LOGS/$id.out"
}
run() { runp "$FIX" "$@"; }
status_is() { [[ $(cat "$WORK/$1/status/$TAG" 2>/dev/null) == "$2"* ]]; }
rows_in() {  # rows_in <out dir> <run id> [status] -> number of run_log rows
    "$PY" - "$1/run_log.csv" "$2" "${3:-}" <<'EOF'
import csv, sys
rows = [r for r in csv.DictReader(open(sys.argv[1])) if r['run_id'] == sys.argv[2] and (not sys.argv[3] or r['status'] == sys.argv[3])]
print(len(rows))
EOF
}
rows() { rows_in "$OUT" "$@"; }
ran_steps() {  # ran_steps <out dir> <run id> -> comma list of the steps that ran (QuPath attempts)
    "$PY" -c "import csv, sys; print(','.join(r['step'] for r in csv.DictReader(open(sys.argv[1])) if r['run_id'] == sys.argv[2] and r['action'] == 'run'))" "$1/run_log.csv" "$2"
}
jq_py() {   # jq_py <json file> <python expression over d> -> prints its value
    "$PY" -c "import json, sys; d = json.load(open(sys.argv[1])); print($2)" "$1"
}
cells_columns() { "$PY" -c "import pyarrow.parquet as pq, sys; print('|'.join(pq.read_schema(sys.argv[1]).names))" "$1"; }
qpstat() { stat -c '%s %.9Y' "$1"; }
copy_fixture() {   # copy_fixture <src project dir> <dest dir>: a private copy without results/
    rm -rf -- "$2"; cp -a "$1" "$2"; rm -rf -- "$2/results"
}

if [[ ${1:-} != --keep-fixture ]]; then
    echo '== rebuilding the fixtures'
    for v in '' pln pln-partial pln-sliver; do
        "$V3DIR/tests/make_fixture.sh" --force ${v:+--variant "$v"} > "$LOGS/make_fixture${v:+_$v}.out" 2>&1 \
            || { cat "$LOGS/make_fixture${v:+_$v}.out"; exit 1; }
        grep -E 'annotations|pLN' "$LOGS/make_fixture${v:+_$v}.out" | sed 's/^/   /'
    done
    rm -rf -- "$OUT"
fi
PD=$OUT/provenance/$TAG/$V3_SITE

echo '== 1. full run: passes I, C, J'
run t-full --pass I,C,J
check 'exit 0 and image DONE' test "$RC" = 0 -a -n "$(status_is t-full DONE && echo y)"
check '12 steps + 2 post hooks passed' test "$(rows t-full PASS)" = 14
check 'step 0 ran the v3 copy (Tissue + pLN regions)' test "$(jq_py "$PD/0.json" "d['script']['source_path'].endswith('v3/qupath/0_cell_detection_panc_v3.groovy')")" = True
check 'step 0 postcondition: Tissue and the contained pLN both hold detections' \
    test "$(jq_py "$PD/0.json" "sorted((r['class'], r['detections'] > 0) for r in d['postcondition']['details']['regions'])")" = "[('Tissue', True), ('pLN', True)]"
check 'the contained pLN was segmented with Tissue, not on its own' \
    test "$(jq_py "$PD/0.json" "(len(d['step_result']['details']['selected']), len(d['step_result']['details']['contained']))")" = '(1, 1)'
check 'records carry site, scripts git commit, finish verdict' \
    test "$(jq_py "$PD/1.json" "(d['site'], len(d['script']['repo_git']['commit'] or ''))")" = "('ws', 40)" -a \
         "$(jq_py "$PD/1.finish.json" "d['status']")" = PASS
check 'step 8 postcondition forbids ROOT and Gal3Norm parents' test "$(jq_py "$PD/8.json" "d['postcondition']['details']['forbidden_parented']")" = 0
cells_columns "$OUT/cells/$TAG.parquet" > "$LOGS/cells_columns.t-full"

echo '== 2. resumability'
run t-resume --pass I,C,J
check 'exit 0, nothing ran' test "$RC" = 0 -a "$(rows t-resume)" = 0
check 'driver says 0 images with work' grep -q '0 with work to do' "$LOGS/t-resume.out"

echo '== 3. footer after a body that ends with method and class definitions'
run t-methods --steps preflight --force preflight --script "preflight=$S/ends_with_methods.groovy"
check 'PASS' test "$RC" = 0 -a "$(rows t-methods PASS)" = 1
run t-methods-restore --steps preflight --force preflight          # put the real preflight record back
check 'real preflight re-run passes, saves nothing (changed=false)' test "$RC" = 0

echo '== 4. early top-level return in step 4'
before=$(qpstat "$QP")
run t-return --pass I,J --force 4 --script "4=$S/early_return.groovy"
check 'exit 1, FAILED at 4' test "$RC" = 1 -a -n "$(status_is t-return 'FAILED at 4' && echo y)"
check 'no sentinel reported' grep -q 'no sentinel' "$LOGS/t-return.out"
check 'chain stopped: only step 4 ran' test "$(rows t-return)" = 1
check 'data.qpdata untouched' test "$(qpstat "$QP")" = "$before"

echo '== 5. step 3 throws'
run t-throw --steps 3 --force 3 --script "3=$S/throws.groovy"
check 'exit 1, FAILED at 3, no sentinel' test "$RC" = 1 -a -n "$(status_is t-throw 'FAILED at 3' && echo y)"
check 'data.qpdata untouched' test "$(qpstat "$QP")" = "$before"

echo '== 6. step 1 fails its postcondition'
run t-postcond --steps 1 --force 1 --script "1=$S/removes_islets.groovy"
check 'exit 1, footer reports islets_present failure' grep -q 'FAIL -- footer: islet-class annotations' "$LOGS/t-postcond.out"
check 'footer did not save' test "$(qpstat "$QP")" = "$before"

echo '== 7. recovery: plain re-run'
run t-recover --pass I,C,J
check 'exit 0 and DONE' test "$RC" = 0 -a -n "$(status_is t-recover DONE && echo y)"
check 're-ran 1, 4 and downstream 3, 8, J_* (not 0, 2)' test "$(ran_steps "$OUT" t-recover)" = '1,4,3,8,J_context,J_annotations,J_islets'

echo '== 8. data.qpdata changed outside run_v3'
touch "$QP"
run t-custody --pass J
check 'exit 1, BLOCKED (custody_broken), logged' test "$RC" = 1 -a -n "$(status_is t-custody 'BLOCKED (custody_broken)' && echo y)" -a "$(rows t-custody BLOCKED)" = 1
"$R" --project "$FIX" --accept-external --images "$TAG" --note 'test: touched on purpose' > "$LOGS/t-accept.out" 2>&1
check 'accept-external recorded' grep -q 'accepted the current data.qpdata' "$LOGS/t-accept.out"
run t-after-accept --pass I,C,J
check 'exit 0; exactly the steps that read the image re-ran (not 1, 4, 0, 2)' \
    test "$RC" = 0 -a "$(ran_steps "$OUT" t-after-accept)" = 'preflight,C_cells,C_geojson,3,8,J_context,J_annotations,J_islets'
cells_columns "$OUT/cells/$TAG.parquet" > "$LOGS/cells_columns.t-after-accept"
check 'cells export after step 3 ran: same columns as before it ran (no signed distances)' \
    cmp -s "$LOGS/cells_columns.t-full" "$LOGS/cells_columns.t-after-accept"
run t-after-accept2 --pass I,C,J
check 'then nothing to run' test "$RC" = 0 -a "$(rows t-after-accept2)" = 0

echo '== 9. image lock held elsewhere'
# one process holds the lock (with `flock file sleep`, killing flock leaves sleep holding the inherited fd)
"$PY" -c 'import fcntl, sys, time; f = open(sys.argv[1], "a"); fcntl.flock(f, fcntl.LOCK_EX); time.sleep(120)' \
    "/scratch/io60/locks/$TAG.lock" &
holder=$!
sleep 1
run t-lock --pass J --force J_islets
check 'exit 1, LOCKED, nothing ran (one LOCKED row)' test "$RC" = 1 -a -n "$(status_is t-lock LOCKED && echo y)" -a "$(rows t-lock LOCKED)" = 1 -a "$(rows t-lock)" = 1
"$R" --project "$FIX" --accept-external --images "$TAG" --note 'must wait for the lock' > "$LOGS/t-lock-accept.out" 2>&1
check 'accept-external refuses while the image lock is held' test $? -ne 0 -a -n "$(grep -l 'image lock' "$LOGS/t-lock-accept.out")"
kill "$holder" 2>/dev/null; wait "$holder" 2>/dev/null

echo '== 10. free-space guard'
run t-disk --pass J --force J_islets --min-root-free-gb 999999
check 'exit 3 with a DISK LOW stop, nothing ran (one STOPPED row)' test "$RC" = 3 -a "$(rows t-disk STOPPED)" = 1 -a "$(rows t-disk)" = 1
check 'message names the floor' grep -q 'DISK LOW' "$LOGS/t-disk.out"

echo '== 11. server-metadata drift (on a copy of the fixture)'
copy_fixture "$FIXDIR" "$FROOT/project_drift"
DP=$FROOT/project_drift/project.qpproj
runp "$DP" t-drift --steps preflight --script "preflight=$S/renames_channel.groovy"
check 'exit 1, DRIFT at preflight' test "$RC" = 1 -a -n "$(status_is t-drift 'DRIFT at preflight' && echo y)"
dj=$FROOT/project_drift/results/v3/provenance/$TAG/$V3_SITE/preflight.drift.json
check 'drift record names server.json and the renamed channel' grep -q '_renamed' "$dj"
runp "$DP" t-drift2 --pass I
check 'next run blocks the image until the drift file is resolved' grep -q 'BLOCKED: unresolved' "$LOGS/t-drift2.out"

echo '== 12. image selection is exact; --force must be selected'
run t-substring --pass J --images fixture --dry-run
check '--images substring: planning refused (exit 2), names the candidates' test "$RC" = 2 -a -n "$(grep -l 'matches no image' "$LOGS/t-substring.out")"
run t-byid --pass J --images id:1 --dry-run
check '--images id:1 resolves the image' test "$RC" = 0 -a -n "$(grep -l "image  $TAG" "$LOGS/t-byid.out")"
run t-force-out --pass I --force 8 --dry-run
check '--force outside the selected passes: refused (exit 2)' test "$RC" = 2 -a -n "$(grep -l 'not among the selected steps' "$LOGS/t-force-out.out")"

echo '== 13. something else writes data.qpdata during a step (on a copy, provenance included)'
rm -rf -- "$FROOT/project_guard"; cp -a "$FIXDIR" "$FROOT/project_guard"   # cp -a keeps mtimes: custody holds
GP=$FROOT/project_guard/project.qpproj
GOUT=$FROOT/project_guard/results/v3
GQP=$FROOT/project_guard/data/1/data.qpdata
runp "$GP" t-guard --steps preflight --force preflight --script "preflight=$S/touches_qpdata.groovy"
check 'exit 1, FAILED at preflight' test "$RC" = 1 -a -n "$(status_is t-guard 'FAILED at preflight' && echo y)"
check 'no sentinel; finish quotes the footer refusal' grep -q 'no sentinel.*REFUSED: step preflight: not saving: data.qpdata changed since the step began' "$LOGS/t-guard.out"
check 'no preflight sentinel on disk' test ! -e "$GOUT/provenance/$TAG/$V3_SITE/preflight.json"
runp "$GP" t-guard2 --steps preflight
check 'the next plan blocks the image (custody broken)' test "$RC" = 1 -a -n "$(status_is t-guard2 'BLOCKED (custody_broken)' && echo y)"
"$R" --project "$GP" --accept-external --images "$TAG" --note 'test: accepted, cells edited' --rerun-from 0 > "$LOGS/t-guard-accept.out" 2>&1
check 'accept-external --rerun-from 0 names what is now stale' grep -q 'Now stale: .*--rerun-from 0' "$LOGS/t-guard-accept.out"
runp "$GP" t-guard-plan --pass C --dry-run
check 'dry run: 0 re-runs (declared by the external record)' grep -qE 'RUN +0 .*rerun-from 0' "$LOGS/t-guard-plan.out"

echo '== 14. the footer never starts a save within --save-margin-min of the deadline (on the copy)'
before=$(qpstat "$GQP")
runp "$GP" t-deadline --steps preflight --force preflight --script "preflight=$S/changes_annotation.groovy" --save-margin-min 100000
check 'exit 1, footer: save skipped' test "$RC" = 1 -a -n "$(grep -l 'save skipped' "$LOGS/t-deadline.out")"
check 'data.qpdata untouched' test "$(qpstat "$GQP")" = "$before"
check 'the failing sentinel records the deadline' test "$(jq_py "$GOUT/provenance/$TAG/$V3_SITE/preflight.json" "(d['status'], d['save']['saved'], d['save']['deadline_epoch_ms'] > 0)")" = "('fail', False, True)"

echo '== 15. v3 script 0 refuses when data.qpdata is not the file the driver recorded (on the copy)'
cfg=$LOGS/fake_step_config.json
"$PY" -c "import json, os, sys; s = os.stat(sys.argv[1]); json.dump({'attempt_id': 'fake', 'qpdata_before': {'exists': True, 'size': s.st_size + 1, 'mtime_ns': s.st_mtime_ns}}, open(sys.argv[2], 'w'))" "$GQP" "$cfg"
before=$(qpstat "$GQP")
JAVA_TOOL_OPTIONS=-Xmx8g "$QUPATH" script -p "$GP" -i "$TAG" -a "$cfg" "$V3DIR/qupath/0_cell_detection_panc_v3.groovy" > "$LOGS/t-s0guard.out" 2>&1
rc=$?
check 'exit non-zero with [detect] REFUSED at start' test "$rc" -ne 0 -a -n "$(grep -l '\[detect\] REFUSED at start' "$LOGS/t-s0guard.out")"
check 'data.qpdata untouched' test "$(qpstat "$GQP")" = "$before"

echo '== 16. disjoint pLN (+ a contained pLN): passes I, C, J with the v3 step 8'
PLN=$FROOT/project_pln/project.qpproj
POUT=$FROOT/project_pln/results/v3
PPD=$POUT/provenance/$TAG/$V3_SITE
rm -rf -- "$POUT"
runp "$PLN" t-pln --pass I,C,J --script "8=$V3_8"
check 'exit 0 and DONE, 14 PASS rows' test "$RC" = 0 -a -n "$(status_is t-pln DONE && echo y)" -a "$(rows_in "$POUT" t-pln PASS)" = 14
check 'step 0 segmented Tissue and the disjoint pLN as regions; the contained pLN with Tissue' \
    test "$(jq_py "$PPD/0.json" "(sorted(s.split(' (')[0] for s in d['step_result']['details']['selected']), [s.split(' (')[0] for s in d['step_result']['details']['contained']])")" = \
         "(['Tissue', \"pLN 'pLN_disjoint'\"], [\"pLN 'pLN_contained'\"])"
check 'every region holds detections (postcondition), none duplicated' \
    test "$(jq_py "$PPD/0.json" "(d['postcondition']['pass'], d['postcondition']['details']['empty_regions'], d['step_result']['details']['duplicate_pairs'])")" = '(True, 0, 0)'
check 'step 8 parented cells to pLN, none to ROOT' \
    test "$(jq_py "$PPD/8.json" "(d['postcondition']['details']['parent_class_counts'].get('pLN', 0) > 0, d['postcondition']['details']['parent_class_counts'].get('ROOT', 0))")" = '(True, 0)'
res=$("$PY" - "$POUT/context/$TAG.parquet" <<'EOF'
import duckdb, sys
px = 0.5078
n, n_pln = duckdb.connect().execute(f"""SELECT count(*), count(*) FILTER (WHERE Parent = 'pLN_disjoint') FROM read_parquet('{sys.argv[1]}')
    WHERE "Centroid X µm" BETWEEN {1350 * px} AND {2032 * px} AND "Centroid Y µm" BETWEEN {200 * px} AND {1500 * px}""").fetchone()
print(n > 0 and n == n_pln, n, n_pln)
EOF
)
check "every cell inside the disjoint pLN is parented to it ($res)" test "${res%% *}" = True

echo '== 17. pLN partly inside Tissue'
runp "$FROOT/project_pln_partial/project.qpproj" t-pln-partial --pass C
PP=$FROOT/project_pln_partial/results/v3/provenance/$TAG/$V3_SITE/preflight.json
check 'exit 1, FAILED at preflight' test "$RC" = 1 -a -n "$(status_is t-pln-partial 'FAILED at preflight' && echo y)"
check 'preflight says partial' test "$(jq_py "$PP" "(d['step_result']['details']['checks']['pln_vs_tissue']['ok'], sorted(d['step_result']['details']['checks']['pln_vs_tissue']['by_kind'].items()))")" = "(False, [('contained', 1), ('partial', 1)])"

echo '== 18. pLN 5 px into Tissue: preflight fails by step 0 rule, step 0 never runs'
SQ=$FROOT/project_pln_sliver/data/1/data.qpdata
runp "$FROOT/project_pln_sliver/project.qpproj" t-pln-sliver --pass C
SP=$FROOT/project_pln_sliver/results/v3/provenance/$TAG/$V3_SITE
check 'exit 1, FAILED at preflight' test "$RC" = 1 -a -n "$(status_is t-pln-sliver 'FAILED at preflight' && echo y)"
check 'preflight says partial and explains the shared area' test "$(jq_py "$SP/preflight.json" "(sorted(d['step_result']['details']['checks']['pln_vs_tissue']['by_kind'].items()), any('trim it off the Tissue' in w for w in d['step_result']['details']['warnings']))")" = "([('contained', 1), ('partial', 1)], True)"
check 'step 0 never ran (no sentinel)' test ! -e "$SP/0.json"

echo '== 19. two provenance sites on one image (on a copy)'
copy_fixture "$FIXDIR" "$FROOT/project_sites"
SPJ=$FROOT/project_sites/project.qpproj
V3_SITE=siteA runp "$SPJ" t-siteA-I --pass I
check 'site A: pass I DONE' test "$RC" = 0
V3_SITE=siteB runp "$SPJ" t-siteB-C --pass C
check 'site B: blocked until the records of site A are merged' test "$RC" = 1 -a -n "$(grep -l 'BLOCKED: records from site(s) siteA are not merged' "$LOGS/t-siteB-C.out")"
V3_SITE=siteB "$R" --project "$SPJ" --merge-provenance --images "$TAG" > "$LOGS/t-merge-B.out" 2>&1
check 'site B merges site A (mirror + custody chain ok)' grep -q 'merged (siteA: +' "$LOGS/t-merge-B.out"
V3_SITE=siteB runp "$SPJ" t-siteB-C2 --pass C
check "site B: pass C runs 0, 2, exports; preflight done by site A" test "$RC" = 0 -a "$(ran_steps "$FROOT/project_sites/results/v3" t-siteB-C2)" = '0,2,C_cells,C_geojson'
V3_SITE=siteA runp "$SPJ" t-siteA-J --pass J
check 'site A: blocked until site B is merged' test "$RC" = 1 -a -n "$(grep -l 'not merged' "$LOGS/t-siteA-J.out")"
V3_SITE=siteA "$R" --project "$SPJ" --merge-provenance --images "$TAG" > "$LOGS/t-merge-A.out" 2>&1
V3_SITE=siteA runp "$SPJ" t-siteA-J2 --pass J
check 'site A: after merging, pass J runs on the file site B left' test "$RC" = 0 -a -n "$(status_is t-siteA-J2 DONE && echo y)"

echo
echo "== per-step timings of the full run (t-full):"
"$PY" - "$OUT/run_log.csv" <<'EOF'
import csv, sys
for r in csv.DictReader(open(sys.argv[1])):
    if r['run_id'] == 't-full':
        rss = f"{int(r['peak_rss_kb']) / 1e6:.1f} GB" if r['peak_rss_kb'] else '-'
        print(f"   {r['step']:<14} {r['status']:<5} {float(r['seconds']):>6.1f} s  peak RSS {rss:>8}  {r['message'][:70]}")
EOF
echo "logs: $LOGS"
((FAILS == 0)) && echo 'ALL FIXTURE TESTS PASSED' || echo "$FAILS FIXTURE TEST(S) FAILED"
exit $((FAILS > 0))
