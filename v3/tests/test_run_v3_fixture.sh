#!/usr/bin/env bash
# test_run_v3_fixture.sh -- end-to-end tests of run_v3.sh on the /scratch fixture project.
#
# Never touches the real project: every run passes --project /scratch/io60/fixture/project/project.qpproj
# (outputs land in /scratch/io60/fixture/project/results/v3). Uses the GPU once (step 0 on 2048x2048 px).
#   1  full run of passes I, C, J from a freshly built fixture        -> every step PASS, image DONE
#   2  second run                                                      -> nothing runs (resumability)
#   3  a step body ending in method/class definitions                  -> footer still runs, PASS
#   4  a step that leaves with a top-level `return`                    -> no sentinel, FAILED, chain stops
#   5  a step that throws                                              -> no sentinel, FAILED, qpdata untouched
#   6  a step that fails its postcondition                             -> footer does not save, FAILED
#   7  plain re-run                                                    -> failed + downstream steps re-run, DONE
#   8  data.qpdata changed outside run_v3 (touch)                      -> BLOCKED; --accept-external unblocks
#   9  image lock held by someone else                                 -> LOCKED, nothing runs
#  10  free-space guard not met                                        -> STOP (exit 3), nothing runs
#  11  a step that changes the server metadata (on a copy)             -> DRIFT recorded, image blocked
# Usage: test_run_v3_fixture.sh [--keep-fixture]   (default rebuilds the fixture first)
set -uo pipefail
V3DIR=$(cd "$(dirname "$0")/.." && pwd)
R=$V3DIR/run_v3.sh
S=$V3DIR/tests/steps
PY=/scratch/io60/venv-v3/bin/python
FIXDIR=/scratch/io60/fixture/project
FIX=$FIXDIR/project.qpproj
OUT=$FIXDIR/results/v3
TAG=fixture_6450.ome.tif
QP=$FIXDIR/data/1/data.qpdata
WORK=/scratch/io60/work/driver
LOGS=/scratch/io60/work/driver/tests/$(date +%Y%m%dT%H%M%S)
mkdir -p "$LOGS"
rm -rf -- "$WORK"/t-*          # run dirs of earlier suite runs (the ids below are fixed)
FAILS=0
# TEST_ALLOW_BUSY=1: other sessions' headless QuPath runs (on their own private projects) would make
# run_v3.sh refuse to start; the fixture is private, so the tests may pass --allow-busy.
EXTRA=()
[[ ${TEST_ALLOW_BUSY:-0} == 1 ]] && EXTRA=(--allow-busy)

check() { local d=$1; shift; if "$@"; then printf '  PASS  %s\n' "$d"; else printf '  FAIL  %s\n' "$d"; FAILS=$((FAILS + 1)); fi; }
run() {   # run <id> <args...> ; sets RC
    local id=$1; shift
    "$R" --project "$FIX" --script-set orig --run-id "$id" "${EXTRA[@]}" "$@" > "$LOGS/$id.out" 2>&1
    RC=$?
    printf '\n[%s] exit %s: run_v3.sh %s\n' "$id" "$RC" "$*"
    grep -E ' (start|PASS|FAIL|BLOCKED|LOCKED|STOP|DRIFT|post )| DONE| FAILED|refus|images:' "$LOGS/$id.out" | grep -v '^== ' | sed 's/^/      /' | head -40
}
status_is() { [[ $(cat "$WORK/$1/status/$TAG" 2>/dev/null) == "$2"* ]]; }
rows() {  # rows <run id> [status] -> number of run_log rows
    "$PY" - "$OUT/run_log.csv" "$1" "${2:-}" <<'EOF'
import csv, sys
rows = [r for r in csv.DictReader(open(sys.argv[1])) if r['run_id'] == sys.argv[2] and (not sys.argv[3] or r['status'] == sys.argv[3])]
print(len(rows))
EOF
}
qpstat() { stat -c '%s %.9Y' "$QP"; }

if [[ ${1:-} != --keep-fixture ]]; then
    echo '== rebuilding the fixture'
    "$V3DIR/tests/make_fixture.sh" --force > "$LOGS/make_fixture.out" 2>&1 || { cat "$LOGS/make_fixture.out"; exit 1; }
    tail -3 "$LOGS/make_fixture.out"
    rm -rf -- "$OUT"
fi

echo '== 1. full run: passes I, C, J'
run t-full --pass I,C,J
check 'exit 0 and image DONE' test "$RC" = 0 -a -n "$(status_is t-full DONE && echo y)"
check '12 steps + 2 post hooks passed' test "$(rows t-full PASS)" = 14

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
before=$(qpstat)
run t-return --pass I,J --force 4 --script "4=$S/early_return.groovy"
check 'exit 1, FAILED at 4' test "$RC" = 1 -a -n "$(status_is t-return 'FAILED at 4' && echo y)"
check 'no sentinel reported' grep -q 'no sentinel' "$LOGS/t-return.out"
check 'chain stopped: only step 4 ran' test "$(rows t-return)" = 1
check 'data.qpdata untouched' test "$(qpstat)" = "$before"

echo '== 5. step 3 throws'
run t-throw --steps 3 --force 3 --script "3=$S/throws.groovy"
check 'exit 1, FAILED at 3, no sentinel' test "$RC" = 1 -a -n "$(status_is t-throw 'FAILED at 3' && echo y)"
check 'data.qpdata untouched' test "$(qpstat)" = "$before"

echo '== 6. step 1 fails its postcondition'
run t-postcond --steps 1 --force 1 --script "1=$S/removes_islets.groovy"
check 'exit 1, footer reports islets_present failure' grep -q 'FAIL -- footer: islet-class annotations' "$LOGS/t-postcond.out"
check 'footer did not save' test "$(qpstat)" = "$before"

echo '== 7. recovery: plain re-run'
run t-recover --pass I,C,J
check 'exit 0 and DONE' test "$RC" = 0 -a -n "$(status_is t-recover DONE && echo y)"
check 're-ran 1, 4 and downstream 3, 8, J_* (not 0, 2)' \
    test "$("$PY" -c "import csv; print(','.join(r['step'] for r in csv.DictReader(open('$OUT/run_log.csv')) if r['run_id']=='t-recover' and r['action']=='run'))")" = '1,4,3,8,J_context,J_annotations,J_islets'

echo '== 8. data.qpdata changed outside run_v3'
touch "$QP"
run t-custody --pass J
check 'exit 1, BLOCKED (custody_broken), logged' test "$RC" = 1 -a -n "$(status_is t-custody 'BLOCKED (custody_broken)' && echo y)" -a "$(rows t-custody BLOCKED)" = 1
"$R" --project "$FIX" --accept-external --images "$TAG" --note 'test: touched on purpose' > "$LOGS/t-accept.out" 2>&1
check 'accept-external recorded' grep -q 'accepted the current data.qpdata' "$LOGS/t-accept.out"
run t-after-accept --pass I,C,J
check 'exit 0, nothing to run' test "$RC" = 0 -a "$(rows t-after-accept)" = 0

echo '== 9. image lock held elsewhere'
# one process holds the lock (with `flock file sleep`, killing flock leaves sleep holding the inherited fd)
"$PY" -c 'import fcntl, sys, time; f = open(sys.argv[1], "a"); fcntl.flock(f, fcntl.LOCK_EX); time.sleep(120)' \
    "/scratch/io60/locks/$TAG.lock" &
holder=$!
sleep 1
run t-lock --pass J --force J_islets
check 'exit 1, LOCKED, nothing ran (one LOCKED row)' test "$RC" = 1 -a -n "$(status_is t-lock LOCKED && echo y)" -a "$(rows t-lock LOCKED)" = 1 -a "$(rows t-lock)" = 1
kill "$holder" 2>/dev/null; wait "$holder" 2>/dev/null

echo '== 10. free-space guard'
run t-disk --pass J --force J_islets --min-root-free-gb 999999
check 'exit 3 with a DISK LOW stop, nothing ran (one STOPPED row)' test "$RC" = 3 -a "$(rows t-disk STOPPED)" = 1 -a "$(rows t-disk)" = 1
check 'message names the floor' grep -q 'DISK LOW' "$LOGS/t-disk.out"

echo '== 11. server-metadata drift (on a copy of the fixture)'
rm -rf -- /scratch/io60/fixture/project_drift
cp -a "$FIXDIR" /scratch/io60/fixture/project_drift
rm -rf -- /scratch/io60/fixture/project_drift/results
"$R" --project /scratch/io60/fixture/project_drift/project.qpproj --script-set orig --run-id t-drift "${EXTRA[@]}" \
    --steps preflight --script "preflight=$S/renames_channel.groovy" > "$LOGS/t-drift.out" 2>&1
RC=$?
grep -E 'DRIFT|PASS|FAIL' "$LOGS/t-drift.out" | sed 's/^/      /' | head -5
check 'exit 1, DRIFT at preflight' test "$RC" = 1 -a -n "$(status_is t-drift 'DRIFT at preflight' && echo y)"
dj=/scratch/io60/fixture/project_drift/results/v3/provenance/$TAG/preflight.drift.json
check 'drift record names server.json and the renamed channel' grep -q '_renamed' "$dj"
"$R" --project /scratch/io60/fixture/project_drift/project.qpproj --script-set orig --run-id t-drift2 "${EXTRA[@]}" \
    --pass I > "$LOGS/t-drift2.out" 2>&1
check 'next run blocks the image until the drift file is resolved' grep -q 'BLOCKED: unresolved' "$LOGS/t-drift2.out"

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
