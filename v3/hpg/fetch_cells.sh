#!/usr/bin/env bash
# fetch_cells.sh -- bring ONE image's finished cells pass back from HiPerGator and adopt its records.
# Run on the workstation; uses the maintainer's ControlMaster (`ssh hpg`).
#
# HiPerGator only ever runs the cells pass (preflight, 0, 2, C_cells, C_geojson). When an image is
# done there, its data file and exports come down ONCE and every later step is local, so no data
# file is ever written at two sites. Nothing is fetched unless the hpg records say every step
# passed; nothing local is overwritten unless the local file is still the stripped file that was
# sent up.
#
# Usage: fetch_cells.sh <entry id> [--out DIR] [--check] [--profile pln]
#
# --profile pln fetches an image of the IO60pLN project (~/IO60pLN, /blue/.../IO60pLN_v3). Its local data
# file is NOT a stripped file: it is the project's own file with the July cells, so it is replaced only
# when it is byte-for-byte the size of its copy in /data/IO60panc2nd_archive/pln_qpdata_2026-10-01 (nothing
# that is not backed up is overwritten) and no QuPath window is open.
set -euo pipefail
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
# shellcheck source=../lib.sh
source "$HERE/../lib.sh"          # PROJ V3
ID=${1:?entry id}; shift || true
OUT="" CHECK=0 PROFILE=""
while (($#)); do
    case $1 in
        --out) OUT=$2; shift ;;
        --check) CHECK=1 ;;
        --profile) PROFILE=$2; shift ;;
        *) echo "unknown option: $1" >&2; exit 2 ;;
    esac
    shift
done
BLUE=/blue/maigan/smith6jt/IO60panc2nd_v3
VARIANT=()
case $PROFILE in
    "") ;;
    pln)
        PROJ=/home/smith6jt/IO60pLN
        BLUE=/blue/maigan/smith6jt/IO60pLN_v3
        PLN_BACKUP=/data/IO60panc2nd_archive/pln_qpdata_2026-10-01
        VARIANT=(--project "$PROJ/project.qpproj" --script "preflight=$HERE/../pln/preflight_pln.groovy"
                 --script "0=$HERE/../pln/0_cell_detection_pln.groovy") ;;
    *) echo "--profile must be pln" >&2; exit 2 ;;
esac
: "${OUT:=$PROJ/results/v3/export}"
REL=${OUT#"$PROJ"/results/}                         # e.g. v3/export

TAG=$(python3 -c 'import json, re, sys
e = next(e for e in json.load(open(sys.argv[1]))["images"] if str(e["entryID"]) == sys.argv[2])
print(re.sub(r"[^A-Za-z0-9._-]+", "_", e["imageName"]))' "$PROJ/project.qpproj" "$ID")
RP=$BLUE/results/$REL

# ---- 1. the hpg records must say every step of the cells pass passed ----
remote_check='
p=$1/provenance/$2/hpg
for s in preflight 0 2 C_cells C_geojson; do
  f=$p/$s.finish.json
  if [ -s "$f" ]; then python3 -c "import json,sys; print(sys.argv[2], json.load(open(sys.argv[1])).get(\"status\"))" "$f" $s; else echo "$s missing"; fi
done
f=$p/C_cells.post.json
if [ -s "$f" ]; then python3 -c "import json,sys; d=json.load(open(sys.argv[1])); print(\"C_cells.post\", d.get(\"status\"))" "$f"; else echo "C_cells.post missing"; fi
stat -c "qpdata %s" $3/data/$4/data.qpdata 2>/dev/null || echo "qpdata-missing"
for f in $1/cells/$2.parquet $1/cell_geojson/$2.geojson; do if [ -s "$f" ]; then echo "file $(stat -c %s "$f") $f"; fi; done
exit 0
'
remote_state=$(ssh hpg bash -s -- "$RP" "$TAG" "$BLUE" "$ID" <<< "$remote_check")
echo "$remote_state" | sed 's/^/  hpg: /'
n_pass=$(grep -c -i -E '^(preflight|0|2|C_cells|C_geojson|C_cells\.post) pass$' <<< "$remote_state" || true)
n_file=$(grep -c '^file ' <<< "$remote_state" || true)
if ((n_pass != 6 || n_file != 2)) || ! grep -q '^qpdata ' <<< "$remote_state"; then
    echo "NOT READY: entry $ID ($TAG): $n_pass of 6 records pass, $n_file of 2 exports present" >&2
    exit 3
fi
((CHECK)) && { echo "ready: entry $ID ($TAG)"; exit 0; }

# ---- 2. nothing local may be lost: the local file is the stripped one, or (pln) it is backed up ----
local_size=$(stat -c %s "$PROJ/data/$ID/data.qpdata")
if [[ $PROFILE == pln ]]; then
    backup_size=$(stat -c %s "$PLN_BACKUP/data/$ID/data.qpdata" 2> /dev/null || echo -1)
    hpg_size=$(sed -n 's/^qpdata //p' <<< "$remote_state")
    if ((local_size == hpg_size)); then
        echo "entry $ID ($TAG): the local file already has HiPerGator's size; nothing to fetch"
        exit 0
    fi
    if ((local_size != backup_size)); then
        echo "REFUSED: local data/$ID/data.qpdata ($local_size bytes) differs from its backup ($backup_size bytes) -- it changed after the backup; resolve by hand" >&2
        exit 4
    fi
    if pgrep -f 'QuPath/bin/QuPath$' > /dev/null; then
        echo "REFUSED for now: a QuPath window is open; the IO60pLN data file is not replaced under it" >&2
        exit 5
    fi
elif ((local_size > 50 * 1024 * 1024)); then
    echo "REFUSED: local data/$ID/data.qpdata is ${local_size} bytes -- not the stripped file; resolve by hand" >&2
    exit 4
fi
mkdir -p "$OUT/cells" "$OUT/cell_geojson" "$OUT/provenance/$TAG" "$OUT/logs/hpg"
log "fetching entry $ID ($TAG) from HiPerGator"
rsync -a --partial hpg:"$RP/cells/$TAG.parquet" "$OUT/cells/"
rsync -a --partial hpg:"$RP/cell_geojson/$TAG.geojson" "$OUT/cell_geojson/"
rsync -a hpg:"$RP/provenance/$TAG/hpg" "$OUT/provenance/$TAG/"
rsync -a hpg:"$RP/logs/run_v3/$TAG" "$OUT/logs/hpg/" 2> /dev/null || true
# the data file last, to a temporary name, then moved into place in one step
rsync -a --partial hpg:"$BLUE/data/$ID/data.qpdata" "$PROJ/data/$ID/data.qpdata.from_hpg"
rsync -a hpg:"$BLUE/data/$ID/server.json" hpg:"$BLUE/data/$ID/summary.json" "$PROJ/data/$ID/"
mv -f "$PROJ/data/$ID/data.qpdata.from_hpg" "$PROJ/data/$ID/data.qpdata"

# ---- 3. adopt the hpg records: the driver validates the chain against the file now on disk ----
"$HERE/../run_v3.sh" --merge-provenance --from hpg --images "id:$ID" --out "$OUT" "${VARIANT[@]}"
log "fetched and merged entry $ID ($TAG): $(stat -c %s "$PROJ/data/$ID/data.qpdata") bytes"
