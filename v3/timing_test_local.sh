#!/usr/bin/env bash
# timing_test_local.sh -- compute pilot, workstation half: pass C on one image with per-step timing.
#   1. preflight; 2. snapshot the post-preflight data/<id> + provenance for the HiPerGator half;
#   3. step 0 (InstanSeg, local GPU) with nvidia-smi dmon; 4. steps 2, C_cells, C_geojson (CPU).
# Usage: timing_test_local.sh <entry id>
set -euo pipefail
source "$(dirname "$0")/lib.sh"
ID=${1:?entry id}
cd "$PROJ"
LOG=$V3/logs/timing_test/id$ID; mkdir -p "$LOG"
ST=$SCRATCH/hpg_stage; mkdir -p "$ST/data" "$ST/results/v3/provenance"
run() { local name=$1; shift
        /usr/bin/time -v -o "$LOG/$name.time" bash scripts/v3/run_v3.sh "$@" --images "id:$ID" > "$LOG/$name.out" 2>&1
        log "timing_test id:$ID $name: exit $? ($(grep -E 'Elapsed \(wall' "$LOG/$name.time" | awk '{print $NF}'))"; }
run preflight --steps preflight
TAG=$(python3 -c "import json,sys; e=[e for e in json.load(open('project.qpproj'))['images'] if e['entryID']==$ID][0]; print(e['imageName'].replace(' - resolution #','_-_resolution_').replace(' ','_'))")
rm -rf "$ST/data/$ID"; cp -a "data/$ID" "$ST/data/$ID"
rm -rf "$ST/results/v3/provenance/$TAG"; cp -a "$V3/provenance/$TAG" "$ST/results/v3/provenance/"
log "timing_test id:$ID: post-preflight snapshot in $ST (tag $TAG)"
nvidia-smi dmon -s um -d 5 -o T > "$LOG/gpu_dmon.log" 2>&1 & DMON=$!
run step0 --steps 0 || true
kill "$DMON" 2>/dev/null || true
run step2_exports --steps 2,C_cells,C_geojson || true
log "timing_test id:$ID done; logs in $LOG"
