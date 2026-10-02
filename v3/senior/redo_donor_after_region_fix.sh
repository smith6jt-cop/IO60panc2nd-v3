#!/usr/bin/env bash
# redo_donor_after_region_fix.sh -- the Senior source workflow for ONE donor whose image was redone from cell
# detection (its Tissue or pLN annotation was corrected in QuPath; run_v3.sh --accept-external --rerun-from 0, then
# --steps preflight,0,2,C_cells,C_geojson,7).
#
# The declared workflow is run as it is, from the one checkout and the one config.ini (run_source_workflow.sh):
#
#   - `ingest` is cohort-wide and refuses a non-empty cells root, so the whole cells root is moved aside and every
#     donor is ingested again from its own cells export (24 of them unchanged);
#   - geometry QC, REDSEA and expression skip a donor whose partition has content, so ONLY this donor's three
#     partitions are moved aside -- the others are left exactly as they are;
#   - a stage also writes a per-donor audit file (data/audit/donors/<name>/donor_id=<donor>) and refuses to
#     overwrite one, so this donor's audit files are moved aside with its partitions (2026-10-02: the expression
#     stage wrote 6347's new partition and then stopped on yesterday's audit file);
#   - the cohort manifest describes each image's files (counts, sizes, hashes; nothing validates them), and this
#     image's files changed, so it is created again from the same specification, beside the run, and put in
#     place when the run has finished.
#
# Nothing is deleted: what is moved aside goes to results/v3/senior_superseded/<stamp>_<donor>_region_fix/.
#
#   redo_donor_after_region_fix.sh --donor 6347 --tag <image tag> --cells-log <the cells pass's log> [--detections N]
#
#   --cells-log    waits until that log shows `C_geojson: PASS` for the tag (the new cells and outlines are exported);
#                  a FAIL / BLOCKED / STOPPED line for the tag stops this script before anything is moved
#   --detections   the detection count the cells pass reported; the new manifest must describe that many outlines
set -euo pipefail
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
WT=/home/smith6jt/IO60panc2nd/Islet-Explorer-Senior-v3
DATA=$WT/core/phenocycler/data
V3=/home/smith6jt/IO60panc2nd/results/v3
DONOR="" TAG="" CELLS_LOG="" DETECTIONS=""
while (($#)); do
    case $1 in
        --donor) DONOR=$2; shift ;;
        --tag) TAG=$2; shift ;;
        --cells-log) CELLS_LOG=$2; shift ;;
        --detections) DETECTIONS=$2; shift ;;
        *) echo "unknown option: $1" >&2; exit 2 ;;
    esac
    shift
done
[[ $DONOR =~ ^[0-9A-Za-z]+$ && -n $TAG && -s $CELLS_LOG ]] || { echo "need --donor, --tag and an existing --cells-log" >&2; exit 2; }
[[ $DONOR != 6457 && $TAG != *6457* ]] || { echo "6457 is not processed" >&2; exit 2; }
say() { echo "$(date +%H:%M:%S) [redo_donor] $*"; }
die() { say "STOPPED: $*"; exit 1; }

# ---- 1. the new cells and outlines are exported ----------------------------------------------------------------
say "waiting for the cells pass to export $TAG (log $CELLS_LOG)"
until grep -q "$TAG: C_geojson: PASS" "$CELLS_LOG"; do
    grep -E "$TAG: .*(FAIL|BLOCKED|STOPPED)" "$CELLS_LOG" && die "the cells pass did not get through"
    grep -q "worker project entries differing" "$CELLS_LOG" && ! grep -q "$TAG: C_geojson: PASS" "$CELLS_LOG" && die "the cells pass ended without exporting the outlines"
    sleep 15
done
grep -q "$TAG: C_cells: PASS" "$CELLS_LOG" || die "C_geojson passed but C_cells did not"
cells=$V3/export/cells/$TAG.parquet outlines=$V3/export/cell_geojson/$TAG.geojson
[[ -s $cells && -s $outlines ]] || die "exports not on disk: $cells / $outlines"
say "exports on disk: $(stat -c '%s bytes %y' "$cells" | cut -c1-40) / $(stat -c '%s bytes %y' "$outlines" | cut -c1-40)"

# ---- 2. move aside what is rebuilt ---------------------------------------------------------------------------
for stage in cells geometry_qc redsea expression; do
    [[ -d $DATA/$stage/donor_id=$DONOR ]] || die "no $stage partition for donor $DONOR in $DATA -- nothing to redo from"
done
SUP=$V3/senior_superseded/$(date +%Y%m%dT%H%M%S)_${DONOR}_region_fix
mkdir -p "$SUP/geometry_qc" "$SUP/redsea" "$SUP/expression" "$SUP/audit"
touch "$SUP/.started"            # anything of another donor's later stages newer than this was written by this run
cp -p "$DATA/qupath_manifest.json" "$SUP/qupath_manifest.json"
for f in ingest_rejects.parquet panel_availability.json; do [[ -e $DATA/$f ]] && cp -p "$DATA/$f" "$SUP/$f"; done
[[ -e $DATA/audit/expression_availability.parquet ]] && cp -p "$DATA/audit/expression_availability.parquet" "$SUP/audit/"
mv -- "$DATA/cells" "$SUP/cells"
for stage in geometry_qc redsea expression; do
    mv -- "$DATA/$stage/donor_id=$DONOR" "$SUP/$stage/donor_id=$DONOR"
done
for a in "$DATA"/audit/donors/*/donor_id="$DONOR"; do
    [[ -d $a ]] || continue
    name=$(basename "$(dirname "$a")")
    mkdir -p "$SUP/audit/donors/$name"
    mv -- "$a" "$SUP/audit/donors/$name/donor_id=$DONOR"
done
say "moved aside to $SUP: the cells root (all donors) and donor $DONOR's geometry_qc / redsea / expression partitions and audit files"

# ---- 3. the manifest, created again beside the run -----------------------------------------------------------
NEW=$DATA/qupath_manifest.json.new
rm -f -- "$NEW"
( "$HERE/senior_v3.sh" --cd core/phenocycler python -m phenocycler.pipeline manifest create \
      --spec data/cohort-spec-v3.json --out data/qupath_manifest.json.new > "$SUP/manifest_create.log" 2>&1 ) &
manifest_pid=$!

# ---- 4. the declared workflow --------------------------------------------------------------------------------
say "source workflow: ingest (all donors), then geometry QC / REDSEA / expression for donor $DONOR"
"$HERE/run_source_workflow.sh" --redsea-workers 2 --geometry-workers 1 || die "the source workflow failed (results/v3/logs/senior_source_workflow_*.log)"

# ---- 5. checks ----------------------------------------------------------------------------------------------
for stage in cells geometry_qc redsea expression; do
    n=$(find "$DATA/$stage" -mindepth 1 -maxdepth 1 -type d -name 'donor_id=*' | wc -l)
    [[ $n -eq 25 ]] || die "$stage holds $n donor partitions, not 25"
    [[ -n $(find "$DATA/$stage/donor_id=$DONOR" -type f -name '*.parquet' -newer "$cells" | head -1) ]] || die "$stage/donor_id=$DONOR was not rebuilt after the new cells export"
done
for stage in geometry_qc redsea expression; do   # the other donors' partitions must be the ones that were there
    for d in "$DATA/$stage"/donor_id=*; do
        [[ $d == */donor_id=$DONOR ]] && continue
        [[ -z $(find "$d" -type f -newer "$SUP/.started" | head -1) ]] || die "$d changed during this run; only donor $DONOR's later stages may"
    done
done
"$HERE/senior_v3.sh" --cd core/phenocycler python - "$DONOR" <<'EOF'
import sys
import pyarrow.parquet as pq
from pathlib import Path
donor = sys.argv[1]
for stage in ("cells", "geometry_qc", "redsea", "expression"):
    files = sorted(Path("data", stage, f"donor_id={donor}").glob("*.parquet"))
    rows = sum(pq.ParquetFile(f).metadata.num_rows for f in files)
    print(f"  {stage:<12} donor {donor}: {rows:,} rows in {len(files)} file(s)")
EOF

wait "$manifest_pid" || die "the manifest could not be created ($SUP/manifest_create.log)"
python3 - "$NEW" "$TAG" "${DETECTIONS:-0}" <<'EOF'
import json, sys
m = json.load(open(sys.argv[1]))
assert len(m["images"]) == 25, f"{len(m['images'])} images in the new manifest"
assert not any("6457" in i["image_id"] or str(i["donor_id"]) == "6457" for i in m["images"]), "6457 is not processed"
tag = sys.argv[2]
hit = [i for i in m["images"] if i["cell_geojson"]["path"].endswith(tag + ".geojson")]
assert len(hit) == 1, f"{len(hit)} manifest images for {tag}"
n = hit[0]["cell_geometry"]["feature_count"]
want = int(sys.argv[3])
assert want == 0 or n == want, f"the manifest describes {n:,} outlines, the cells pass reported {want:,}"
print(f"  manifest: 25 images; {tag}: {n:,} cell outlines")
EOF
mv -- "$NEW" "$DATA/qupath_manifest.json"
say "manifest replaced (the earlier one is in $SUP)"
say "done"
