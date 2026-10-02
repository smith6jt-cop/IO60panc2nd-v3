#!/usr/bin/env bash
# islet_classifiers.sh -- the per-image islet classifier, from the fit to the maintainer's script 1.
#
# The maintainer, 2026-10-01: "All I wanted was the pixel classifier to be more accurate and not filter out so many
# cells while adapting to batch effects that changed the classifier for each image ... Then the rest of the workflow
# works without major changes." So the ONLY thing that changes is which classifier script 1 loads:
# qupath/1_panc_objects_perimage.groovy is scripts/1_panc_objects.groovy with that one line changed.
#
#   fit      islet_classifier_fit.py per image -> <project>/classifiers/pixel_classifiers/Islets_<tag>.json
#            and IsletHormone_<tag>.json (the hormone's own outline, for Hormone Fragments) -- new files per image;
#            none of the maintainer's classifiers is touched -- and a report in $V3/islet_classifier/fit/<tag>.json.
#            Needs the image's cells export, its structure objects and the source workflow's REDSEA partition.
#   markers  IsletHormone_<tag>.json again from the thresholds of an existing fit report (no fit, no image read)
#   trial    on a PRIVATE copy of each image's data file (the master project is not written): first the maintainer's
#            script 1 as it is (Islets2), exported to trial/out_islets2/reviewed_islets/ for comparison; then script 1
#            with the image's own classifier, step 1b (outlines erode to the cell outlines; Hormone Fragments --
#            maintainer, 2026-10-02), step 4 (names), exported to trial/out/reviewed_islets/. The outlines as script 1
#            left them are exported to trial/out/script1_outlines/ before 1b runs (step 1_outlines). The copies end
#            holding the new islets, so trial/project.qpproj can be opened to look at them.
#   check    QuPath's objects from the trial: the same cells as the fit expected? and the counts beside Islets2's
#   montage  pictures of both sets of objects on the image, for looking without QuPath
#   table    the checks (or, before a trial, the fit reports) of the given images as one table
#            check, montage and table read the trial copies' objects; with --master they read the MASTER project's
#            (the exports `apply` wrote under $V3/export) and write beside the trial's: fit/<tag>.master.check.json,
#            pictures_master/<donor>.png. Nothing is written to the project by any of the three.
#   apply    script 1 with the image's own classifier, then 1b and 4, on the master project. Only after the maintainer
#            has looked.
#   applied  what `apply` left on the master project, per image, read from the driver's own records
#            ($V3/export/provenance): script 1's outlines, the Islets and Hormone Fragments 1b made of them, what it
#            had to repair or skip, and the classes the export then found. Reads only.
#
# The test set is every non-diabetic / control donor (maintainer, 2026-10-01: "use all non-diabetic/controls for
# testing classifier since the endocrine differences with disease interfere"):  --ids "$ND_IDS"
#
# Usage: islet_classifiers.sh <fit|markers|trial|check|montage|table|apply|applied> --ids "<entry ids>" | --nd [--master] [--jobs N] [--lanes N]
set -euo pipefail
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
# shellcheck source=lib.sh
source "$HERE/lib.sh"            # PROJ V3
STAGE=${1:?stage}; shift
IDS="" JOBS=4 LANES=2 MASTER=0
# disease.status == ND in the maintainer's ~/IO60panc2nd/donor_metadata_panc.xlsx (all nine controls):
# HDL112 (38), HDL115 (17), 6374 (1), 6442 (8), 6461 (9), 6476 (11), 6482 (4), 6493 (6), 6539 (5)
ND_IDS="38 17 1 8 9 11 4 6 5"
while (($#)); do
    case $1 in
        --ids) IDS=$2; shift ;;
        --nd) IDS=$ND_IDS ;;
        --master) MASTER=1 ;;
        --jobs) JOBS=$2; shift ;;
        --lanes) LANES=$2; shift ;;
        *) echo "unknown option: $1" >&2; exit 2 ;;
    esac
    shift
done
[[ -n $IDS ]] || { echo "--ids \"<entry ids>\" or --nd is required" >&2; exit 2; }
EXPORT=$V3/export
OUT=$V3/islet_classifier
TRIAL=$OUT/trial
# whose objects check / montage / table read: the trial copies', or (--master) the master project's as `apply` exported them
OBJ=$TRIAL/out SUF="" PICS=$OUT/pictures
((MASTER)) && { OBJ=$EXPORT; SUF=.master; PICS=$OUT/pictures_master; }
PYS=/home/smith6jt/IO60panc2nd/Islet-Explorer-Senior/.venv/bin/python
SENIOR_DATA=/home/smith6jt/IO60panc2nd/Islet-Explorer-Senior-v3/core/phenocycler/data
SPEC=$SENIOR_DATA/cohort-spec-v3.json
STEP1=$HERE/qupath/1_panc_objects_perimage.groovy
ids_csv=$(sed 's/[^ ]\+/id:&/g; s/ \+/,/g' <<< "$IDS")

# id -> "tag<TAB>image file<TAB>donor"; 6457 is not processed and is refused by name
describe() { python3 - "$PROJ/project.qpproj" "$SPEC" "$1" <<'EOF'
import json, re, sys
from urllib.parse import unquote, urlparse
entry = next((e for e in json.load(open(sys.argv[1]))["images"] if str(e["entryID"]) == sys.argv[3]), None)
assert entry is not None, f"entry id {sys.argv[3]} is not in the project"
name = entry["imageName"]
assert "6457" not in name, "6457 is not processed"
donor = next(str(r["donor_id"]) for r in json.load(open(sys.argv[2]))["images"] if r["image_id"] == name)
builder = entry["serverBuilder"]
uri = builder.get("uri") or builder.get("builder", {}).get("uri")
print(re.sub(r"[^A-Za-z0-9._-]+", "_", name), unquote(urlparse(uri).path), donor, sep="\t")
EOF
}

case $STAGE in
fit)
    mkdir -p "$OUT/fit" "$OUT/logs"
    running=0
    for i in $IDS; do
        IFS=$'\t' read -r tag image donor < <(describe "$i")
        redsea=$SENIOR_DATA/redsea/donor_id=$donor/data_0.parquet
        [[ -s $redsea ]] || { log "$tag: no REDSEA partition for donor $donor yet ($redsea) -- the source workflow has not reached it; skipped"; continue; }
        [[ -s $EXPORT/cells/$tag.parquet ]] || die "$tag: no cells export"
        log "fit $tag (donor $donor)"
        rm -f "$OUT/fit/$tag.check.json" "$OUT/fit/$tag.master.check.json"   # a check belongs to the classifier it checked
        choom -n 1000 -- nice -n 5 "$PYS" "$HERE/islet_classifier_fit.py" --tag "$tag" --image "$image" --cells "$EXPORT/cells/$tag.parquet" \
            --redsea "$redsea" --structures "$EXPORT/structures" \
            --classifier "$PROJ/classifiers/pixel_classifiers/Islets_$tag.json" --report "$OUT/fit/$tag.json" \
            --mask "$OUT/fit/$tag.mask.npy" --objects "$OUT/fit/$tag.objects.npy" > "$OUT/logs/fit_$tag.log" 2>&1 &
        ((++running % JOBS == 0)) && wait
    done
    wait
    for i in $IDS; do
        IFS=$'\t' read -r tag _ _ < <(describe "$i")
        [[ -s $OUT/fit/$tag.json && -s $PROJ/classifiers/pixel_classifiers/Islets_$tag.json ]] && log "$tag: classifier written" \
            || log "$tag: NO classifier -- see $OUT/logs/fit_$tag.log: $(tail -n 1 "$OUT/logs/fit_$tag.log" 2> /dev/null)"
    done ;;
markers)
    for i in $IDS; do
        IFS=$'\t' read -r tag _ _ < <(describe "$i")
        [[ -s $OUT/fit/$tag.json ]] || die "$tag: no fit report -- run fit first"
        "$PYS" "$HERE/islet_classifier_fit.py" --tag "$tag" --classifier "$PROJ/classifiers/pixel_classifiers/Islets_$tag.json" \
            --report "$OUT/fit/$tag.json" --marker-only
    done ;;
trial)
    mkdir -p "$TRIAL/data"
    python3 - "$PROJ/project.qpproj" "$TRIAL/project.qpproj" $IDS <<'EOF'
import json, os, sys
src, dst, ids = sys.argv[1], sys.argv[2], set(sys.argv[3:])
p = json.load(open(src))
have = {str(e["entryID"]) for e in json.load(open(dst))["images"]} if os.path.exists(dst) else set()
p["images"] = [e for e in p["images"] if str(e["entryID"]) in ids | have]
assert ids <= {str(e["entryID"]) for e in p["images"]}, "an entry id is not in the project"
json.dump(p, open(dst, "w"), indent=2)
EOF
    for d in classifiers scripts resources; do [[ -e $PROJ/$d ]] && ln -sfn "$PROJ/$d" "$TRIAL/$d"; done
    for i in $IDS; do
        IFS=$'\t' read -r tag _ _ < <(describe "$i")
        [[ -s $PROJ/classifiers/pixel_classifiers/Islets_$tag.json ]] || die "$tag: no classifier -- run fit first"
        [[ -s $PROJ/classifiers/pixel_classifiers/IsletHormone_$tag.json ]] || die "$tag: no IsletHormone_ classifier -- run markers (or fit)"
        [[ -s $TRIAL/data/$i/data.qpdata ]] && { log "trial: data/$i already copied"; continue; }
        mkdir -p "$TRIAL/data/$i"
        cp -p "$PROJ/data/$i/data.qpdata" "$PROJ/data/$i/server.json" "$PROJ/data/$i/summary.json" "$TRIAL/data/$i/"
        log "trial: copied data/$i ($(stat -c %s "$TRIAL/data/$i/data.qpdata") bytes)"
    done
    # the comparison first: the maintainer's script 1 exactly as it is, for the images that do not have it yet
    todo=""
    for i in $IDS; do
        IFS=$'\t' read -r tag _ _ < <(describe "$i")
        [[ -s $TRIAL/out_islets2/reviewed_islets/$tag.geojson ]] || todo+=" id:$i"
    done
    if [[ -n $todo ]]; then
        "$HERE/run_v3.sh" --project "$TRIAL/project.qpproj" --out "$TRIAL/out_islets2" --images "$(tr ' ' ',' <<< "${todo# }")" \
            --steps preflight,1,4,R_export --script "1=$PROJ/scripts/1_panc_objects.groovy" --force 1 --lanes-light "$LANES" --lanes-cpu 1 \
            --mem-budget 300 --heap preflight=24 --heap 4=32 --heap R_export=32
    fi
    # script 1; its outlines exported as it drew them (1_outlines); then they erode to the cells (1b); names (4); export
    "$HERE/run_v3.sh" --project "$TRIAL/project.qpproj" --out "$TRIAL/out" --images "$ids_csv" \
        --steps preflight,1,1_outlines,1b,4,R_export --script "1=$STEP1" --force 1 --lanes-light "$LANES" --lanes-cpu 1 \
        --mem-budget 300 --heap preflight=24 --heap 1_outlines=32 --heap 4=32 --heap R_export=32 ;;
check)
    running=0
    for i in $IDS; do
        IFS=$'\t' read -r tag _ donor < <(describe "$i")
        baseline=()
        [[ -s $TRIAL/out_islets2/reviewed_islets/$tag.geojson ]] && baseline=(--baseline "$TRIAL/out_islets2/reviewed_islets/$tag.geojson")
        before=()
        [[ -s $OBJ/reviewed_islets/$tag.geojson ]] || { log "$tag: no objects in $OBJ/reviewed_islets -- skipped"; continue; }
        [[ -s $OBJ/script1_outlines/$tag.geojson ]] && before=(--script1-objects "$OBJ/script1_outlines/$tag.geojson" \
            --cell-outlines "$EXPORT/cell_geojson/$tag.geojson")
        "$PYS" "$HERE/islet_classifier_check.py" --tag "$tag" --fit "$OUT/fit/$tag.json" --fit-objects "$OUT/fit/$tag.objects.npy" \
            --qupath-objects "$OBJ/reviewed_islets/$tag.geojson" "${before[@]}" "${baseline[@]}" --cells "$EXPORT/cells/$tag.parquet" \
            --redsea "$SENIOR_DATA/redsea/donor_id=$donor/data_0.parquet" --structures "$EXPORT/structures" \
            --out "$OUT/fit/$tag$SUF.check.json" > "$OUT/logs/check${SUF}_$tag.log" 2>&1 &
        ((++running % JOBS == 0)) && wait
    done
    wait
    for i in $IDS; do
        IFS=$'\t' read -r tag _ _ < <(describe "$i")
        [[ -s $OUT/fit/$tag$SUF.check.json ]] && log "$tag: checked" || log "$tag: check FAILED -- $OUT/logs/check${SUF}_$tag.log"
    done ;;
montage)
    mkdir -p "$PICS"
    running=0
    for i in $IDS; do
        IFS=$'\t' read -r tag image donor < <(describe "$i")
        baseline=()
        [[ -s $TRIAL/out_islets2/reviewed_islets/$tag.geojson ]] && baseline=(--baseline "$TRIAL/out_islets2/reviewed_islets/$tag.geojson")
        before=()
        [[ -s $OBJ/reviewed_islets/$tag.geojson ]] || { log "$tag: no objects in $OBJ/reviewed_islets -- skipped"; continue; }
        [[ -s $OBJ/script1_outlines/$tag.geojson ]] && before=(--script1 "$OBJ/script1_outlines/$tag.geojson")
        "$PYS" "$HERE/islet_classifier_montage.py" --tag "$tag" --image "$image" --objects "$OBJ/reviewed_islets/$tag.geojson" \
            "${before[@]}" --cell-outlines "$EXPORT/cell_geojson/$tag.geojson" \
            "${baseline[@]}" --cells "$EXPORT/cells/$tag.parquet" --redsea "$SENIOR_DATA/redsea/donor_id=$donor/data_0.parquet" \
            --structures "$EXPORT/structures" --anchors "$EXPORT/islet_anchors.csv" --label "$donor" \
            --out "$PICS/$donor.png" > "$OUT/logs/montage${SUF}_$tag.log" 2>&1 &
        ((++running % JOBS == 0)) && wait
    done
    wait
    ls -la "$PICS" ;;
apply)
    for i in $IDS; do
        IFS=$'\t' read -r tag _ _ < <(describe "$i")
        [[ -s $PROJ/classifiers/pixel_classifiers/Islets_$tag.json ]] || die "$tag: no classifier -- run fit first"
    done
    for i in $IDS; do
        IFS=$'\t' read -r tag _ _ < <(describe "$i")
        [[ -s $PROJ/classifiers/pixel_classifiers/IsletHormone_$tag.json ]] || die "$tag: no IsletHormone_ classifier -- run markers (or fit)"
    done
    # --out is the v3 export tree: its records are the ones the master data files were last written under
    "$HERE/run_v3.sh" --out "$EXPORT" --images "$ids_csv" --steps preflight,1,1_outlines,1b,4,R_export --script "1=$STEP1" --lanes-light "$LANES" --lanes-cpu 1 \
        --mem-budget 300 --heap preflight=24 --heap 1_outlines=32 --heap 4=32 --heap R_export=32 ;;
applied)
    "$PYS" - "$EXPORT/provenance" $(for i in $IDS; do describe "$i" | cut -f1,3 | tr '\t' '='; done) <<'EOF'
import json, sys
from pathlib import Path
prov = Path(sys.argv[1])
cols = ["donor", "script 1 outlines", "Islets", "one whole cell", "two or more", "outlines split", "pieces without a whole cell dropped",
        "whole cells", "Hormone Fragments", "no cell and no hormone (removed)", "overlays repaired", "openings part by part", "openings skipped",
        "exported classes", "1b status", "export status"]
print(*cols, sep="\t")
def latest(tag, step):          # the newest record of a step on this image, whichever site wrote it
    found = sorted(prov.glob(f"{tag}/*/{step}.json"), key=lambda f: f.stat().st_mtime)
    return json.loads(found[-1].read_text()) if found else None
for item in sys.argv[2:]:
    tag, donor = item.rsplit("=", 1)
    b, e = latest(tag, "1b"), latest(tag, "R_export")
    if b is None:
        print(donor, "1b has not run", sep="\t"); continue
    d = (b.get("step_result") or {}).get("details") or {}
    by = ((e or {}).get("step_result") or {}).get("details", {}).get("by_class", {})
    print(donor, d.get("islets_from_script_1"), d.get("islets"), d.get("islets_one_whole_cell"), d.get("islets_two_or_more"), d.get("outlines_split"),
          d.get("pieces_without_a_whole_cell_dropped"), d.get("whole_cells"), d.get("hormone_fragments"), d.get("outlines_without_cell_or_hormone_removed"),
          d.get("overlays_repaired"), d.get("openings_done_part_by_part"), d.get("openings_skipped"),
          ", ".join(f"{k} {v}" for k, v in sorted(by.items())), b.get("status"), (e or {}).get("status", "not exported"), sep="\t")
EOF
    ;;
table)
    "$PYS" - "$OUT/fit" "$SUF" $(for i in $IDS; do describe "$i" | cut -f1,3 | tr '\t' '='; done) <<'EOF'
import json, sys
from pathlib import Path
fit, suffix = Path(sys.argv[1]), sys.argv[2]
cols = ["donor", "objects", "CD99", "INS", "GCG", "SST", "islets", "1 cell", "2-3", "4-10", ">10", "hormone fragments", "cells inside", "neighbours inside", "%",
        "enclosed non-endocrine inside", "endocrine outside", "%", "confirmed endocrine outside", "%", "of confirmed"]
print(*cols, sep="\t")
num = lambda v: "not used" if v is None else f"{v:.0f}"
for item in sys.argv[3:]:
    tag, donor = item.rsplit("=", 1)
    check, report = fit / f"{tag}{suffix}.check.json", fit / f"{tag}.json"
    if not report.exists():
        print(donor, "no classifier", sep="\t"); continue
    r = json.loads(report.read_text()); t = r["thresholds"]
    thr = [num(t[k]) for k in ("CD99", "INS", "GCG", "SST")]
    rows = []
    if check.exists():
        c = json.loads(check.read_text())
        if "islets2" in c:
            rows.append(("Islets2 (script 1 as it was)", [""] * 4, c["islets2"], ""))
        if "script_1" in c:
            rows.append(("this image's own classifier, script 1", thr, c["script_1"], ""))
            rows.append(("... then eroded to the cell outlines (1b)", thr, c["qupath"], c["fragments"]["hormone_fragments"]))
        else:
            rows.append(("this image's own classifier", thr, c["qupath"], ""))
    else:
        rows.append(("this image's own classifier (expected)", thr, r["after_script_1"], ""))
    for name, th, s, frag in rows:
        by = s.get("objects_by_cells", {})
        print(donor, name, *th, s["objects_with_cells"], by.get("1", ""), by.get("2-3", ""), by.get("4-10", ""), by.get(">10", ""), frag,
              s["cells_inside"], s["inside_other_at_the_edge"], s["inside_other_at_the_edge_percent"],
              s["inside_other_enclosed"], s["endocrine_outside"], s["endocrine_outside_percent"], s["confirmed_outside"], s["confirmed_outside_percent"],
              r["cells"]["confirmed"], sep="\t")
EOF
    ;;
*) echo "unknown stage: $STAGE" >&2; exit 2 ;;
esac
