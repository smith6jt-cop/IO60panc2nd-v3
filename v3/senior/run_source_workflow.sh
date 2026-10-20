#!/usr/bin/env bash
# run_source_workflow.sh -- the Senior source workflow on the v3 cohort, in the v3 worktree:
#
#     QuPath cohort manifest -> ingest -> geometry QC -> REDSEA -> expression
#
# It is the declared workflow, run as the repository says to run it (from core/phenocycler, which owns
# the only config.ini), through senior_v3.sh so the worktree's own code and data are used. Nothing here
# types a cell: the pipeline stops after `expression`.
#
# Run it detached, because it outlives a terminal command:
#     setsid nohup scripts/v3/senior/run_source_workflow.sh --wait-for-structures > /dev/null 2>&1 &
#
#   --redsea-workers N      default 5: one REDSEA worker held 19 GB on 6539 (550k cells, 59 channels); the
#                           largest images are four times that size, and QuPath may be open beside it
#   --geometry-workers N    default 2 (config.ini's own)
#   --wait-for-structures   start only when no local structure driver (run_v3.sh --pass S) is running
#   --force                 passed to the pipeline: rerun stages whose output is already present
#
# A donor is done when its partition has content; a rerun skips what is there.
set -euo pipefail
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
WT=/home/smith6jt/IO60panc2nd/Islet-Explorer-Senior-v3
LOGS=/home/smith6jt/IO60panc2nd/results/v3/logs
REDSEA=5 GEOMETRY=2 WAIT=0
EXTRA=()
while (($#)); do
    case $1 in
        --redsea-workers) REDSEA=$2; shift ;;
        --geometry-workers) GEOMETRY=$2; shift ;;
        --wait-for-structures) WAIT=1 ;;
        --force) EXTRA+=(--force) ;;
        *) echo "unknown option: $1" >&2; exit 2 ;;
    esac
    shift
done
mkdir -p "$LOGS"
LOG=$LOGS/senior_source_workflow_$(date +%Y%m%dT%H%M%S).log
say() { echo "$(date -Is) [run_source_workflow.sh] $*" | tee -a "$LOG"; }

[[ -s $WT/core/phenocycler/data/qupath_manifest.json ]] || { say "no manifest in the worktree: create it first"; exit 2; }
python3 - "$WT/core/phenocycler/data/qupath_manifest.json" <<'EOF' | tee -a "$LOG"
import json, sys
images = json.load(open(sys.argv[1]))["images"]
donors = sorted(str(i["donor_id"]) for i in images)
assert "6457" not in donors and not any("6457" in i["image_id"] for i in images), "6457 is not processed"
print(f"manifest: {len(images)} images, donors {' '.join(donors)}")
EOF

if ((WAIT)); then
    say "waiting for the local structure drivers to finish"
    while pgrep -f "run_v3[.]sh --pass S" > /dev/null; do sleep 30; done
fi
say "starting: --redsea-workers $REDSEA --geometry-workers $GEOMETRY ${EXTRA[*]:-} (log: $LOG)"
set +e
choom -n 1000 -- "$HERE/senior_v3.sh" --cd core/phenocycler python -m phenocycler.pipeline run --config config.ini \
    --pipelined --through expression --redsea-workers "$REDSEA" --geometry-workers "$GEOMETRY" "${EXTRA[@]}" >> "$LOG" 2>&1
rc=$?
set -e
say "finished with exit code $rc"
for stage in cells geometry_qc redsea expression; do
    say "$stage: $(find "$WT/core/phenocycler/data/$stage" -mindepth 1 -maxdepth 1 -type d -name 'donor_id=*' 2>/dev/null | wc -l) donor partitions"
done
exit $rc
