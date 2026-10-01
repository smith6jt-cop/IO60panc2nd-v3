#!/usr/bin/env bash
# islets_all.sh -- islet objects for the MASTER project, after the maintainer chose the detector (Gate 1).
#
# The same four steps the 3-image comparison ran on private copies (islet_comparison.sh), now on the
# project itself, through the same driver and footer:
#
#   detect    step 1 with the chosen detector, step 4 (names), R_export (the pixel objects as GeoJSON)
#   resolve   resolve_islets_to_cells.py: each pixel object -> the outline of its whole cells
#             (a cell with a nucleus and more than half its area inside); the pixel objects are kept
#             on disk (pixel_objects/<tag>.geojson) but are NOT put back into the project
#   import    R_import (the cell-based objects replace the pixel objects), then step 4 again, so the
#             names and numbers belong to the objects that are there
#
# After this the maintainer reviews the 25 images in QuPath; each reviewed image is then accepted with
#   run_v3.sh --accept-external --images id:<n> --note "islets reviewed" --out <export>
# and pass J runs on it.
#
# Usage: islets_all.sh <detect|resolve|import|all> --detector A|B --images "id:5,id:33,..." [--lanes N]
#   A = scripts/v3/qupath/1_panc_objects_v3.groovy (per-image thresholds; needs islet_anchors.csv)
#   B = scripts/1_panc_objects.groovy (the maintainer's Islets2 script, unchanged)
set -euo pipefail
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
# shellcheck source=lib.sh
source "$HERE/lib.sh"            # PROJ V3
STAGE=${1:?stage}; shift
DETECTOR="" IMAGES="" LANES=2      # step 1 peaks at 55-98 GB per image: two at a time beside an open QuPath
while (($#)); do
    case $1 in
        --detector) DETECTOR=$2; shift ;;
        --images) IMAGES=$2; shift ;;
        --lanes) LANES=$2; shift ;;
        *) echo "unknown option: $1" >&2; exit 2 ;;
    esac
    shift
done
[[ $DETECTOR == A || $DETECTOR == B ]] || { echo "--detector A|B is required (the maintainer's choice; there is no default)" >&2; exit 2; }
[[ -n $IMAGES ]] || { echo "--images is required" >&2; exit 2; }
[[ $IMAGES != *6457* ]] || { echo "6457 is not processed" >&2; exit 2; }
EXPORT=$V3/export
PYS=/home/smith6jt/IO60panc2nd/Islet-Explorer-Senior/.venv/bin/python

variant=()
[[ $DETECTOR == B ]] && variant=(--script "1=$PROJ/scripts/1_panc_objects.groovy")
driver() {   # every call names the detector's script: a later step is refused when its prerequisite ran another one
    "$HERE/run_v3.sh" --out "$EXPORT" --images "$IMAGES" --lanes-light "$LANES" --lanes-cpu 1 --mem-budget 300 \
        --heap preflight=24 --heap 4=40 --heap R_export=40 --heap R_import=40 "${variant[@]}" "$@"
}
tags() { "$PYS" - "$PROJ/project.qpproj" "$IMAGES" <<'EOF'
import json, re, sys
want = {x.split(":", 1)[1] for x in sys.argv[2].split(",") if x.startswith("id:")}
for e in json.load(open(sys.argv[1]))["images"]:
    if str(e["entryID"]) in want:
        print(re.sub(r"[^A-Za-z0-9._-]+", "_", e["imageName"]))
EOF
}

detect() {
    if [[ $DETECTOR == A ]]; then
        driver --steps preflight,1,4,R_export --override "1:ANCHORS_CSV='$EXPORT/islet_anchors.csv'" --override "1:SAVE_CLASSIFIERS=false"
    else
        driver --steps preflight,1,4,R_export
    fi
}
resolve() {
    mkdir -p "$EXPORT/pixel_objects" "$EXPORT/resolved_islets"
    local tag
    while read -r tag; do
        # R_export wrote the pixel objects where R_import reads; the original is kept aside once
        [[ -s $EXPORT/pixel_objects/$tag.geojson ]] || cp -p "$EXPORT/reviewed_islets/$tag.geojson" "$EXPORT/pixel_objects/$tag.geojson"
        "$PYS" "$HERE/resolve_islets_to_cells.py" --objects "$EXPORT/pixel_objects/$tag.geojson" \
            --cells "$EXPORT/cell_geojson/$tag.geojson" --out "$EXPORT/resolved_islets/$tag.geojson" \
            --members "$EXPORT/resolved_islets/$tag.members.tsv" --summary "$EXPORT/resolved_islets/$tag.summary.json" \
            > "$EXPORT/resolved_islets/$tag.log"
        cp -p "$EXPORT/resolved_islets/$tag.geojson" "$EXPORT/reviewed_islets/$tag.geojson"
        log "$tag: $("$PYS" -c 'import json, sys; s = json.load(open(sys.argv[1])); print(s["pixel_objects"], "pixel objects ->", s["resolved_objects"], "cell-based objects,", s["member_cells"], "member cells")' "$EXPORT/resolved_islets/$tag.summary.json")"
    done < <(tags)
}
import_objects() {
    driver --steps R_import --force R_import
    driver --steps 4 --force 4
}

case $STAGE in
    detect) detect ;;
    resolve) resolve ;;
    import) import_objects ;;
    all) detect; resolve; import_objects ;;
    *) echo "unknown stage: $STAGE" >&2; exit 2 ;;
esac
