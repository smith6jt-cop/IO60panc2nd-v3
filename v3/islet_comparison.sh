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
# THE MEMBER RULE (2026-10-01). Area alone was rejected by the maintainer on these three images ("Too many
# false positives handling segmentations around endocrine-positive cells"). `resolve` therefore takes
# --rule: `area` (the rejected rule, kept in out/resolved/ as the record) or one of islet_cell_test.py's
# rules, written to out/resolved_<rule>/. --alt-rules adds other rules to the same image as reference
# classes ("Alt: <rule>"), so they can be compared in place. --values redsea reads the source workflow's
# REDSEA-corrected values from <--redsea-root>/redsea/donor_id=<donor>/data_0.parquet.
#
# Usage: islet_comparison.sh <stage> [variant] --ids "5 33 37" [--rule R] [--alt-rules "R2 R3"]
#                            [--values raw|redsea --redsea-root DIR] [--review-dir NAME]
#   stages: prepare | detect A|B | resolve A|B | import A|B | refresh A|B | assemble | table
set -euo pipefail
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
# shellcheck source=lib.sh
source "$HERE/lib.sh"            # PROJ V3
STAGE=${1:?stage}; shift
VARIANT=""
[[ ${1:-} == A || ${1:-} == B ]] && { VARIANT=$1; shift; }
IDS="" RULE="" ALT_RULES="" VALUES=raw REDSEA_ROOT="" REVIEW_DIR=review
while (($#)); do
    case $1 in
        --ids) IDS=$2; shift ;;
        --rule) RULE=$2; shift ;;
        --alt-rules) ALT_RULES=$2; shift ;;
        --values) VALUES=$2; shift ;;
        --redsea-root) REDSEA_ROOT=$2; shift ;;
        --review-dir) REVIEW_DIR=$2; shift ;;
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

SPEC=/home/smith6jt/IO60panc2nd/Islet-Explorer-Senior-v3/core/phenocycler/data/cohort-spec-v3.json
donor_of() { python3 -c 'import json, sys
name = next(e["imageName"] for e in json.load(open(sys.argv[1]))["images"] if str(e["entryID"]) == sys.argv[2])
print(next(str(r["donor_id"]) for r in json.load(open(sys.argv[3]))["images"] if r["image_id"] == name))' "$PROJ/project.qpproj" "$1" "$SPEC"; }
resolved_dir() { [[ $1 == area ]] && echo resolved || echo "resolved_$1"; }

tag_of() { python3 -c 'import json, re, sys
e = next(e for e in json.load(open(sys.argv[1]))["images"] if str(e["entryID"]) == sys.argv[2])
print(re.sub(r"[^A-Za-z0-9._-]+", "_", e["imageName"]))' "$PROJ/project.qpproj" "$1"; }

run_driver() {   # run_driver <variant> <extra run_v3.sh args...>
    local v=$1; shift
    local -a variant=()
    # B's step 1 is the maintainer's script on EVERY call: a later step is refused when the plan names
    # a different script than the one its prerequisite ran
    [[ $v == B ]] && variant=(--script "1=$PROJ/scripts/1_panc_objects.groovy")
    # heaps sized to what these steps were measured to use (R_export 10 GB, step 4 <= 19 GB, preflight
    # <= 16 GB peak RSS), so the comparison fits beside the cells pass inside the host's 300 GB budget
    "$HERE/run_v3.sh" --project "$REV/$v/project.qpproj" --out "$REV/$v/out" --images "$ids_csv" \
        --lanes-light 1 --lanes-cpu 1 --mem-budget 300 \
        --heap preflight=24 --heap 4=32 --heap R_export=32 --heap R_import=32 "${variant[@]}" "$@"
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
        run_driver B --steps preflight,1,4,R_export
    fi ;;
resolve)
    [[ -n $VARIANT ]] || { echo "resolve needs A or B" >&2; exit 2; }
    [[ -n $RULE ]] || { echo "resolve needs --rule (area, or a rule of islet_cell_test.py): the maintainer's choice, no default" >&2; exit 2; }
    [[ $VALUES == raw || -n $REDSEA_ROOT ]] || { echo "--values redsea needs --redsea-root" >&2; exit 2; }
    out=$REV/$VARIANT/out
    mkdir -p "$out/pixel_objects"
    for i in $IDS; do
        tag=$(tag_of "$i")
        # R_export wrote the pixel objects where R_import reads; keep the original aside once
        [[ -s $out/pixel_objects/$tag.geojson ]] || cp -p "$out/reviewed_islets/$tag.geojson" "$out/pixel_objects/$tag.geojson"
        for r in $RULE $ALT_RULES; do
            d=$out/$(resolved_dir "$r"); mkdir -p "$d"
            test_args=()
            if [[ $r != area ]]; then
                values=(--values "$VALUES")
                [[ $VALUES == redsea ]] && values+=(--redsea "$REDSEA_ROOT/redsea/donor_id=$(donor_of "$i")/data_0.parquet")
                "$PYS" "$HERE/islet_cell_test.py" --rule "$r" --cells "$EXPORT/cells/$tag.parquet" --objects "$out/pixel_objects/$tag.geojson" \
                    --structures "$EXPORT/structures" --tag "$tag" "${values[@]}" \
                    --out "$d/$tag.cell_test.parquet" --summary "$d/$tag.cell_test.json" > "$d/$tag.cell_test.log"
                test_args=(--cell-test "$d/$tag.cell_test.parquet")
            elif [[ -s $d/$tag.summary.json ]]; then
                continue                                   # the rejected area-only objects are the record: never rebuilt
            fi
            "$PYS" "$HERE/resolve_islets_to_cells.py" --objects "$out/pixel_objects/$tag.geojson" \
                --cells "$EXPORT/cell_geojson/$tag.geojson" --out "$d/$tag.geojson" "${test_args[@]}" \
                --members "$d/$tag.members.tsv" --summary "$d/$tag.summary.json" > "$d/$tag.log"
            log "$VARIANT $tag, rule $r: $("$PYS" -c 'import json, sys; s = json.load(open(sys.argv[1])); print(s["pixel_objects"], "pixel objects ->", s["resolved_objects"], "cell-based objects,", s["member_cells"], "member cells")' "$d/$tag.summary.json")"
        done
        # what gets imported: the chosen rule's objects, the pixel objects as a reference class, and each
        # alternative rule's objects as its own reference class
        "$PYS" - "$out" "$tag" "$REF_CLASS" "$(resolved_dir "$RULE")" $(for r in $ALT_RULES; do echo "$r=$(resolved_dir "$r")"; done) <<'EOF'
import json, sys, uuid
out, tag, ref, chosen = sys.argv[1:5]
alts = [a.split("=", 1) for a in sys.argv[5:]]
ALT_COLOR = {"top-marker": [255, 0, 255], "cd99-level": [255, 255, 0], "cd99-or-hormone-top": [0, 255, 0], "area": [255, 255, 255]}
P = json.load(open(f"{out}/pixel_objects/{tag}.geojson")); P = P["features"] if isinstance(P, dict) else P
R = json.load(open(f"{out}/{chosen}/{tag}.geojson"))["features"]
feats = list(R)
for f in P:
    cls = ((f.get("properties") or {}).get("classification") or {}).get("name")
    if cls in ("Tissue", "pLN"):
        continue
    g = dict(f)
    g["id"] = str(uuid.uuid5(uuid.NAMESPACE_URL, "pixel-object:" + str(f.get("id"))))
    g["properties"] = {"objectType": "annotation", "classification": {"name": ref, "color": [150, 150, 150]},
                       "name": "pixel " + str(cls)}
    feats.append(g)
n_alt = {}
for rule, sub in alts:
    for f in json.load(open(f"{out}/{sub}/{tag}.geojson"))["features"]:
        g = dict(f)
        g["id"] = str(uuid.uuid5(uuid.NAMESPACE_URL, f"alt:{rule}:" + str(f.get("id"))))
        m = f["properties"]["measurements"]
        g["properties"] = {"objectType": "annotation", "classification": {"name": f"Alt: {rule}", "color": ALT_COLOR.get(rule, [255, 255, 255])},
                           "name": f"{rule}: {m['n_cells']} cells", "measurements": {"n_cells": m["n_cells"]}}
        feats.append(g)
        n_alt[rule] = n_alt.get(rule, 0) + 1
dst = f"{out}/reviewed_islets/{tag}.geojson"
json.dump({"type": "FeatureCollection", "features": feats}, open(dst, "w"))
print(f"{dst}: {len(R)} objects of the chosen rule + {len(feats) - len(R) - sum(n_alt.values())} pixel objects + alternatives {n_alt}")
EOF
    done ;;
import)
    [[ -n $VARIANT ]] || { echo "import needs A or B" >&2; exit 2; }
    # two runs, in this order: the driver follows the step table, where 4 precedes R_import, and the
    # names and numbers must be given to the objects AFTER they are replaced by the cell-based ones
    run_driver "$VARIANT" --steps R_import --force R_import
    run_driver "$VARIANT" --steps 4 --force 4 ;;
refresh)
    # the private copies were made right after segmentation; once the master's cells carry their
    # measurements (same detections, step 2 only adds values) the copies take that file, so the cells
    # in the review project can be read like those of any other image. The driver is told (the data
    # file changed outside it), then the cell-based objects are imported again.
    [[ -n $VARIANT ]] || { echo "refresh needs A or B" >&2; exit 2; }
    for i in $IDS; do
        cp -p "$PROJ/data/$i/data.qpdata" "$PROJ/data/$i/summary.json" "$REV/$VARIANT/data/$i/"
        log "$VARIANT: data/$i replaced by the master's measured file ($(stat -c %s "$REV/$VARIANT/data/$i/data.qpdata") bytes)"
    done
    run_driver "$VARIANT" --accept-external --note "cells replaced by the master project's measured cells (same detections)" || true
    run_driver "$VARIANT" --steps R_import --force R_import
    run_driver "$VARIANT" --steps 4 --force 4 ;;
assemble)
    mkdir -p "$REV/$REVIEW_DIR/data"
    # its OWN class list (never a link: QuPath writes classes.json when the project is saved, and a
    # link would let a review session rewrite the master project's classes); scripts are read-only links
    [[ -L $REV/$REVIEW_DIR/classifiers ]] && unlink "$REV/$REVIEW_DIR/classifiers"
    mkdir -p "$REV/$REVIEW_DIR/classifiers"
    ln -sfn "$PROJ/scripts" "$REV/$REVIEW_DIR/scripts"
    python3 - "$PROJ/project.qpproj" "$REV" "$REF_CLASS" "$REVIEW_DIR" "$ALT_RULES" $IDS <<'EOF'
import copy, json, os, shutil, sys
src, rev, ref, review, alt_rules, ids = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5].split(), sys.argv[6:]
classes = json.load(open(os.path.join(os.path.dirname(src), "classifiers", "classes.json")))
rgb = lambda r, g, b: (255 << 24 | r << 16 | g << 8 | b) - (1 << 32)        # QuPath's packed ARGB, signed
extra = [(ref, -6908266)] + [(f"Alt: {r}", {"top-marker": rgb(255, 0, 255), "cd99-level": rgb(255, 255, 0),
                                              "cd99-or-hormone-top": rgb(0, 255, 0), "area": rgb(255, 255, 255)}.get(r, rgb(255, 255, 255)))
                             for r in alt_rules]
for name, color in extra:
    if not any(c["name"] == name for c in classes["pathClasses"]):
        classes["pathClasses"].append({"name": name, "color": color})
json.dump(classes, open(os.path.join(rev, review, "classifiers", "classes.json"), "w"), indent=2)
p = json.load(open(src)); entries = {str(e["entryID"]): e for e in p["images"]}
label = {"A": "A threshold detector", "B": "B Islets2"}
images = []
for i in ids:
    for k, v in enumerate(("A", "B")):
        e = copy.deepcopy(entries[i]); new_id = (k + 1) * 1000 + int(i)
        short = e["imageName"].split("_Scan")[0]
        e["entryID"] = new_id; e["imageName"] = f"{short} | {label[v]}"; e.pop("randomizedName", None)
        dst = os.path.join(rev, review, "data", str(new_id)); os.makedirs(dst, exist_ok=True)
        for f in ("data.qpdata", "server.json", "summary.json"):
            shutil.copy2(os.path.join(rev, v, "data", i, f), os.path.join(dst, f))
        images.append(e)
p["images"] = images; p["lastID"] = max(e["entryID"] for e in images)
p["uri"] = "file:" + os.path.join(rev, review, "project.qpproj")
json.dump(p, open(os.path.join(rev, review, "project.qpproj"), "w"), indent=2)
print("review project:", os.path.join(rev, review, "project.qpproj"), [e["imageName"] for e in images])
EOF
    ;;
table)
    "$PYS" - "$REV" "$PROJ/project.qpproj" "${RULE:-area} $ALT_RULES" $IDS <<'EOF'
import json, os, re, sys
rev, proj, rules, ids = sys.argv[1], sys.argv[2], sys.argv[3].split(), sys.argv[4:]
P = {str(e["entryID"]): e for e in json.load(open(proj))["images"]}
rows = []
for i in ids:
    tag = re.sub(r"[^A-Za-z0-9._-]+", "_", P[i]["imageName"]); short = P[i]["imageName"].split("_Scan")[0]
    for v in ("A", "B"):
        for rule in rules:
            f = os.path.join(rev, v, "out", "resolved" if rule == "area" else "resolved_" + rule, tag + ".summary.json")
            if not os.path.exists(f):
                continue
            s = json.load(open(f)); c = s["resolved_by_cell_count"]; d = s["resolved_by_class_drotar"]
            rows.append([short, v, rule, s["pixel_objects"], s["resolved_objects"], sum(s["dropped_no_member_by_pixel_class"].values()),
                         d.get("Islet", 0), d.get("EndocrineCluster", 0), d.get("EndocrineSingle", 0),
                         c.get("1 cell", 0), c.get("2-3", 0), c.get("4-10", 0), c.get(">10", 0),
                         s["member_cells"], s.get("majority_cells_failing_the_cell_test", 0), s["majority_cells_without_nucleus"]])
head = ["image", "detector", "member rule", "pixel objects", "cell-based objects", "dropped (no member)", "Islet >=35um", "Cluster", "Single <=10um",
        "1 cell", "2-3 cells", "4-10 cells", ">10 cells", "member cells", "majority cells failing the cell test", "majority cells w/o nucleus"]
out = os.path.join(rev, "comparison_table.tsv" if rules == ["area"] else "comparison_table_rules.tsv")
with open(out, "w") as h:
    h.write("\t".join(head) + "\n")
    for r in rows:
        h.write("\t".join(map(str, r)) + "\n")
print("\t".join(head)); [print("\t".join(map(str, r))) for r in rows]; print("->", out)
EOF
    ;;
*) echo "unknown stage: $STAGE" >&2; exit 2 ;;
esac
