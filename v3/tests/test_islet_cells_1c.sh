#!/usr/bin/env bash
# test_islet_cells_1c.sh -- tests of v3 step 1c (qupath/1c_islet_objects_to_cells_v3.groovy) on PRIVATE copies of
# the run_v3 fixture (/scratch/io60/fixture/project, built by make_fixture.sh; its cells come from the fixture
# suite's pass C). Never touches the real project or the shared fixture: every QuPath run gets its own project copy
# under $WORK with data/1 COPIED. Each result is re-checked independently by check_islet_cells_1c.py (shapely).
#
#   1  script 1 (v3 detector), then 1c with the defaults, saved          -> every check passes
#   2  1c with SINGLE_RULE=caliper and COVER_MIN=0.3 on the same script-1 state -> every check passes
#   3  1c with APPLY=false                                              -> saved image keeps the pixel objects
#   4  1c on its own output                                             -> refused, data.qpdata untouched
#   5  synthetic scenarios (islet_cells_1c_synthetic.groovy)            -> the known answers
#
# The fixture's anchors CSV (ANCHORS) predates the current fixture annotations, so script 1 runs with
# ANCHOR_CHECK=warn here; 1c itself does not use the anchors. Usage: test_islet_cells_1c.sh [--keep]
set -uo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
V3=$(cd "$HERE/.." && pwd)
S1=$V3/qupath/1_panc_objects_v3.groovy
S1C=$V3/qupath/1c_islet_objects_to_cells_v3.groovy
EXPORT=$HERE/islet_cells_1c_export.groovy
SYNTH=$HERE/islet_cells_1c_synthetic.groovy
CHECK=$HERE/check_islet_cells_1c.py
PY=/scratch/io60/venv-v3/bin/python
QUPATH=/home/smith6jt/QuPath/bin/QuPath
FIXSRC=/scratch/io60/fixture/project
ANCHORS=${ANCHORS:-/scratch/io60/work/isletC2/fixture/fixture_anchors.csv}
WORK=${WORK:-/scratch/io60/work/islet1c/tests}
IMG=fixture_6450.ome.tif
export JAVA_TOOL_OPTIONS="${JAVA_TOOL_OPTIONS_1C:--Xmx16g -Djava.util.concurrent.ForkJoinPool.common.parallelism=16}"
KEEP=0; [[ ${1:-} == --keep ]] && KEEP=1
case $WORK in /home/smith6jt/IO60panc2nd*|/scratch/io60/fixture*) echo "WORK must be private" >&2; exit 2 ;; esac
rm -rf -- "$WORK"; mkdir -p "$WORK/logs"
FAILS=0
ok()   { printf '  PASS  %s\n' "$1"; }
bad()  { printf '  FAIL  %s\n' "$1"; FAILS=$((FAILS + 1)); }
check() { local d=$1; shift; if "$@"; then ok "$d"; else bad "$d"; fi; }

copy_fixture() {   # copy_fixture <dir>: private project with data/1 COPIED (never a symlink)
    local d=$1
    mkdir -p "$d/data"
    cp -p "$FIXSRC/project.qpproj" "$d/"
    cp -p "$FIXSRC/data/metadata.json" "$d/data/" 2>/dev/null || true
    cp -rp "$FIXSRC/data/1" "$d/data/1"
    cp -rp "$FIXSRC/classifiers" "$d/"
    ln -s /home/smith6jt/IO60panc2nd/scripts "$d/scripts"   # read-only: lib/common.groovy
}
qp() {   # qp <log> <project dir> [--save] <script> [KEY=VALUE ...]  -> QuPath's exit code
    local log=$1 proj=$2; shift 2
    local save=(); if [[ $1 == --save ]]; then save=(-s); shift; fi
    local script=$1; shift
    local a=(); for x in "$@"; do a+=(-a "$x"); done
    choom -n 1000 -- "$QUPATH" script -p "$proj/project.qpproj" -i "$IMG" "${save[@]}" "${a[@]}" "$script" > "$log" 2>&1
}
stamp() { stat -c '%s %Y' "$1"; }
run_check() {   # run_check <name> <before> <after> <run> [checker args]
    local name=$1 b=$2 a=$3 r=$4; shift 4
    "$PY" "$CHECK" --before "$b" --after "$a" --run "$r" --json-out "$r.check.json" "$@" > "$WORK/logs/$name.check.log" 2>&1
    local rc=$?
    sed 's/^/      /' "$WORK/logs/$name.check.log" | grep -E 'FAIL|check1c'
    return $rc
}

echo '== 0. private fixture copy + script 1 (v3 detector)'
P=$WORK/p_main
copy_fixture "$P"
qp "$WORK/logs/s1.log" "$P" --save "$S1" ANCHORS_CSV="$ANCHORS" ANCHOR_CHECK=warn SAVE_CLASSIFIERS=false
check 'script 1 ran and saved' grep -q '^\[01v3\] fixture_6450.ome.tif: \[Islet' "$WORK/logs/s1.log"
grep -E '^\[01v3\] fixture_6450.ome.tif: \[' "$WORK/logs/s1.log" | sed 's/^/      /'
cp -p "$P/data/1/data.qpdata" "$WORK/after_s1.qpdata"
qp "$WORK/logs/export_before.log" "$P" "$EXPORT" OUT="$WORK/before"
check 'fixture has cells (run the fixture suite pass C first otherwise)' \
    "$PY" -c "import json,sys; sys.exit(0 if json.load(open('$WORK/before.info.json'))['n_cells'] > 0 else 1)"

echo '== 1. 1c with the defaults (COVER_MIN 0.5, SINGLE_RULE cells, TOUCH_GAP_PX 1)'
qp "$WORK/logs/1c_t1.log" "$P" --save "$S1C" WRITE_GEOJSON="$WORK/t1/run.geojson"
rc=$?
check '1c exit 0' test $rc -eq 0
grep -E '^\[1c\] (fixture|  dropped|  old class|verification)' "$WORK/logs/1c_t1.log" | sed 's/^/      /'
qp "$WORK/logs/export_t1.log" "$P" "$EXPORT" OUT="$WORK/t1/after"
check 'independent check (defaults)' run_check t1 "$WORK/before" "$WORK/t1/after" "$WORK/t1/run"
cp -p "$P/data/1/data.qpdata" "$WORK/after_1c.qpdata"

echo '== 2. SINGLE_RULE=caliper, COVER_MIN=0.3 (same script-1 state)'
cp -p "$WORK/after_s1.qpdata" "$P/data/1/data.qpdata"
qp "$WORK/logs/1c_t2.log" "$P" --save "$S1C" WRITE_GEOJSON="$WORK/t2/run.geojson" SINGLE_RULE=caliper COVER_MIN=0.3
check '1c exit 0 (caliper, 0.3)' test $? -eq 0
grep -E '^\[1c\] fixture' "$WORK/logs/1c_t2.log" | sed 's/^/      /'
qp "$WORK/logs/export_t2.log" "$P" "$EXPORT" OUT="$WORK/t2/after"
check 'independent check (caliper, 0.3)' run_check t2 "$WORK/before" "$WORK/t2/after" "$WORK/t2/run" --single-rule caliper --cover-min 0.3

echo '== 3. APPLY=false'
cp -p "$WORK/after_s1.qpdata" "$P/data/1/data.qpdata"
qp "$WORK/logs/1c_t3.log" "$P" --save "$S1C" WRITE_GEOJSON="$WORK/t3/run.geojson" APPLY=false
check '1c exit 0 (APPLY=false)' test $? -eq 0
qp "$WORK/logs/export_t3.log" "$P" "$EXPORT" OUT="$WORK/t3/after"
check 'independent check (not applied: pixel objects kept)' run_check t3 "$WORK/before" "$WORK/t3/after" "$WORK/t3/run" --not-applied

echo '== 4. 1c on its own output is refused'
cp -p "$WORK/after_1c.qpdata" "$P/data/1/data.qpdata"
s0=$(stamp "$P/data/1/data.qpdata")
qp "$WORK/logs/1c_t4.log" "$P" --save "$S1C"
rc=$?
check '1c refuses (exit != 0)' test $rc -ne 0
check 'refusal names the reason' grep -q 'already carry n_cells' "$WORK/logs/1c_t4.log"
check 'data.qpdata untouched' test "$(stamp "$P/data/1/data.qpdata")" == "$s0"

echo '== 5. synthetic scenarios'
S=$WORK/p_synth
copy_fixture "$S"
mkdir -p "$WORK/t5"
qp "$WORK/logs/synth.log" "$S" --save "$SYNTH" PRIVATE_PROJECT="$S" EXPECT="$WORK/t5/expected.json"
check 'synthetic image built' grep -q '^\[synthetic1c\]' "$WORK/logs/synth.log"
qp "$WORK/logs/export_t5b.log" "$S" "$EXPORT" OUT="$WORK/t5/before" CELLS=all
qp "$WORK/logs/1c_t5.log" "$S" --save "$S1C" WRITE_GEOJSON="$WORK/t5/run.geojson"
check '1c exit 0 (synthetic)' test $? -eq 0
grep -E '^\[1c\] (fixture|  dropped|  old class)' "$WORK/logs/1c_t5.log" | sed 's/^/      /'
qp "$WORK/logs/export_t5a.log" "$S" "$EXPORT" OUT="$WORK/t5/after" CELLS=all
check 'independent check + expected answers (synthetic)' run_check t5 "$WORK/t5/before" "$WORK/t5/after" "$WORK/t5/run" --expect "$WORK/t5/expected.json"

echo
if ((FAILS)); then echo "test_islet_cells_1c: $FAILS FAILED (logs in $WORK/logs)"; exit 1; fi
echo "test_islet_cells_1c: all passed (logs in $WORK/logs)"
((KEEP)) || rm -rf -- "$WORK/p_main" "$WORK/p_synth" "$WORK"/*.qpdata
exit 0
