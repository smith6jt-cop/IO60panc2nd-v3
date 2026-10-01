#!/usr/bin/env bash
# islet_comparison.sh -- the 3-image islet comparison the maintainer asked for (2026-09-30):
# both detectors, both resolved to whole cells, side by side in one QuPath project to look at.
#
#   A  the CD99 + hormone threshold detector (qupath/1_panc_objects_v3.groovy; per-image thresholds,
#      background measured over Tissue as drawn -- no DAPI mask)
#   B  the maintainer's own scripts/1_panc_objects.groovy (Islets2 classifier, size-adaptive rules)
#
# Each runs on a PRIVATE copy of the image's data file (the master project is never written), through
# the same driver and footer as everything else. The pixel objects are then resolved to cells by
# resolve_islets_to_cells.py (Senior's majority-of-area rule) and re-imported; the pixel objects are
# kept beside them as a reference class so the two can be toggled.
#
# Usage: islet_comparison.sh <stage> [variant] --ids "5 33 37"
#   stages: prepare | detect A|B | resolve A|B | import A|B | assemble | table
set -euo pipefail
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
# shellcheck source=lib.sh
source "$HERE/lib.sh"            # PROJ V3
STAGE=${1:?stage}; shift
VARIANT=""
[[ ${1:-} == A || ${1:-} == B ]] && { VARIANT=$1; shift; }
IDS=""
while (($#)); do
    case $1 in
        --ids) IDS=$2; shift ;;
        *) echo "unknown option: $1" >&2; exit 2 ;;
    esac
    shift
done
[[ -n $IDS ]] || { echo "--ids \"<entry ids>\" is required" >&2; exit 2; }
REV=$V3/review/islet_comparison
EXPORT=$V3/export                                  # the new cells pass: cell_geojson/, islet_anchors.csv
PYS=/home/smith6jt/IO60panc2nd/Islet-Explorer-Senior/.venv/bin/python
REF_CLASS="Pixel object (before cells)"
ids_csv=$(sed 's/[^ ]\+/id:&/g; s/ \+/,/g' <<< "$IDS")

tag_of() { python3 -c 'import json, re, sys
e = next(e for e in json.load(open(sys.argv[1]))["images"] if str(e["entryID"]) == sys.argv[2])
print(re.sub(r"[^A-Za-z0-9._-]+", "_", e["imageName"]))' "$PROJ/project.qpproj" "$1"; }

run_driver() {   # run_driver <variant> <extra run_v3.sh args...>
    local v=$1; shift
    "$HERE/run_v3.sh" --project "$REV/$v/project.qpproj" --out "$REV/$v/out" --images "$ids_csv" \
        --lanes-light 1 --lanes-cpu 1 --mem-budget 300 "$@"
}

case $STAGE in
prepare)
    for v in A B; do
        mkdir -p "$REV/$v/data"
        python3 - "$PROJ/project.qpproj" "$REV/$v/project.qpproj" $IDS <<'EOF'
import json, sys
src, dst, ids = sys.argv[1], sys.argv[2], set(sys.argv[3:])
p = json.load(open(src))
p["images"] = [e for e in p["images"] if str(e["entryID"]) in ids]
assert len(p["images"]) == len(ids), "an entry id is not in the project"
json.dump(p, open(dst, "w"), indent=2)
EOF
        for d in classifiers scripts resources; do [[ -e $PROJ/$d ]] && ln -sfn "$PROJ/$d" "$REV/$v/$d"; done
        for i in $IDS; do
            [[ -s $REV/$v/data/$i/data.qpdata ]] && { log "$v: data/$i already copied"; continue; }
            mkdir -p "$REV/$v/data/$i"
            cp -p "$PROJ/data/$i/data.qpdata" "$PROJ/data/$i/server.json" "$PROJ/data/$i/summary.json" "$REV/$v/data/$i/"
            log "$v: copied data/$i ($(stat -c %s "$REV/$v/data/$i/data.qpdata") bytes)"
        done
    done ;;
detect)
    [[ -n $VARIANT ]] || { echo "detect needs A or B" >&2; exit 2; }
    if [[ $VARIANT == A ]]; then
        run_driver A --steps preflight,1,4,R_export \
            --override "1:ANCHORS_CSV='$EXPORT/islet_anchors.csv'" --override "1:SAVE_CLASSIFIERS=false"
    else
        run_driver B --steps preflight,1,4,R_export --script "1=$PROJ/scripts/1_panc_objects.groovy"
    fi ;;
resolve)
    [[ -n $VARIANT ]] || { echo "resolve needs A or B" >&2; exit 2; }
    out=$REV/$VARIANT/out
    mkdir -p "$out/pixel_objects" "$out/resolved"
    for i in $IDS; do
        tag=$(tag_of "$i")
        # R_export wrote the pixel objects where R_import reads; keep the original aside once
        [[ -s $out/pixel_objects/$tag.geojson ]] || cp -p "$out/reviewed_islets/$tag.geojson" "$out/pixel_objects/$tag.geojson"
        "$PYS" "$HERE/resolve_islets_to_cells.py" --objects "$out/pixel_objects/$tag.geojson" \
            --cells "$EXPORT/cell_geojson/$tag.geojson" --out "$out/resolved/$tag.geojson" \
            --members "$out/resolved/$tag.members.tsv" --summary "$out/resolved/$tag.summary.json" > "$out/resolved/$tag.log"
        # what gets imported: the resolved objects + the pixel objects as a reference class
        "$PYS" - "$out/pixel_objects/$tag.geojson" "$out/resolved/$tag.geojson" "$out/reviewed_islets/$tag.geojson" "$REF_CLASS" <<'EOF'
import json, sys, uuid
pix, res, dst, ref = sys.argv[1:5]
P = json.load(open(pix)); P = P["features"] if isinstance(P, dict) else P
R = json.load(open(res))["features"]
out = list(R)
for f in P:
    cls = ((f.get("properties") or {}).get("classification") or {}).get("name")
    if cls in ("Tissue", "pLN"):
        continue
    g = dict(f)
    g["id"] = str(uuid.uuid5(uuid.NAMESPACE_URL, "pixel-object:" + str(f.get("id"))))
    g["properties"] = {"objectType": "annotation", "classification": {"name": ref, "color": [150, 150, 150]},
                       "name": "pixel " + str(cls)}
    out.append(g)
json.dump({"type": "FeatureCollection", "features": out}, open(dst, "w"))
print(f"{dst}: {len(R)} resolved + {len(out) - len(R)} reference objects")
EOF
    done ;;
import)
    [[ -n $VARIANT ]] || { echo "import needs A or B" >&2; exit 2; }
    run_driver "$VARIANT" --steps R_import,4 --force 4 ;;
assemble)
    mkdir -p "$REV/review/data"
    for d in classifiers scripts resources; do [[ -e $PROJ/$d ]] && ln -sfn "$PROJ/$d" "$REV/review/$d"; done
    python3 - "$PROJ/project.qpproj" "$REV" $IDS <<'EOF'
import copy, json, os, shutil, sys
src, rev, ids = sys.argv[1], sys.argv[2], sys.argv[3:]
p = json.load(open(src)); entries = {str(e["entryID"]): e for e in p["images"]}
label = {"A": "A threshold detector", "B": "B Islets2"}
images = []
for i in ids:
    for k, v in enumerate(("A", "B")):
        e = copy.deepcopy(entries[i]); new_id = (k + 1) * 1000 + int(i)
        short = e["imageName"].split("_Scan")[0]
        e["entryID"] = new_id; e["imageName"] = f"{short} | {label[v]}"; e.pop("randomizedName", None)
        dst = os.path.join(rev, "review", "data", str(new_id)); os.makedirs(dst, exist_ok=True)
        for f in ("data.qpdata", "server.json", "summary.json"):
            shutil.copy2(os.path.join(rev, v, "data", i, f), os.path.join(dst, f))
        images.append(e)
p["images"] = images; p["lastID"] = max(e["entryID"] for e in images)
json.dump(p, open(os.path.join(rev, "review", "project.qpproj"), "w"), indent=2)
print("review project:", os.path.join(rev, "review", "project.qpproj"), [e["imageName"] for e in images])
EOF
    ;;
table)
    "$PYS" - "$REV" "$PROJ/project.qpproj" $IDS <<'EOF'
import json, os, re, sys
rev, proj, ids = sys.argv[1], sys.argv[2], sys.argv[3:]
P = {str(e["entryID"]): e for e in json.load(open(proj))["images"]}
rows = []
for i in ids:
    tag = re.sub(r"[^A-Za-z0-9._-]+", "_", P[i]["imageName"]); short = P[i]["imageName"].split("_Scan")[0]
    for v in ("A", "B"):
        f = os.path.join(rev, v, "out", "resolved", tag + ".summary.json")
        if not os.path.exists(f):
            continue
        s = json.load(open(f)); c = s["resolved_by_cell_count"]; d = s["resolved_by_class_drotar"]
        rows.append([short, v, s["pixel_objects"], s["resolved_objects"], sum(s["dropped_no_member_by_pixel_class"].values()),
                     d.get("Islet", 0), d.get("EndocrineCluster", 0), d.get("EndocrineSingle", 0),
                     c.get("1 cell", 0), c.get("2-3", 0), c.get("4-10", 0), c.get(">10", 0),
                     s["member_cells"], s["majority_cells_without_nucleus"]])
head = ["image", "detector", "pixel objects", "cell-based objects", "dropped (no whole cell)", "Islet >=35um", "Cluster", "Single <=10um",
        "1 cell", "2-3 cells", "4-10 cells", ">10 cells", "member cells", "majority cells w/o nucleus"]
out = os.path.join(rev, "comparison_table.tsv")
with open(out, "w") as h:
    h.write("\t".join(head) + "\n")
    for r in rows:
        h.write("\t".join(map(str, r)) + "\n")
print("\t".join(head)); [print("\t".join(map(str, r))) for r in rows]; print("->", out)
EOF
    ;;
*) echo "unknown stage: $STAGE" >&2; exit 2 ;;
esac
