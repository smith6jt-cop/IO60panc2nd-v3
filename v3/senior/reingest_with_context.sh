#!/usr/bin/env bash
# reingest_with_context.sh -- after pass J: give the Senior source workflow each cell's place.
#
# Pass J (signed distance, parenting, exports) writes, per image, the context table: every cell's Parent (Islet_N /
# Cluster_N / Single_N / Annotation (pLN) / Annotation (Tissue)), its centroid and its signed distance to each
# annotation class. The source workflow ran before the islets existed, on the cells table alone; this runs the
# declared workflow again with the (cells, context) pair:
#
#   1. the cohort specification names each image's context table (make_cohort_spec.py --context) and the manifest
#      is created from it -- BEFORE the run, because ingest takes the context paths from the manifest;
#   2. `ingest` is cohort-wide and refuses a non-empty cells root: the cells root is moved aside and every donor is
#      ingested again, now joined to its context;
#   3. geometry QC and REDSEA read cell outlines and pixels, not the Parent: their partitions stay and are skipped;
#   4. `expression` carries cell_region and the islet numbers, so the expression root and its per-donor audit files
#      are moved aside and it is rebuilt for all donors (seconds per donor).
#
# `--force` cannot do 4: a donor stage refuses to write over an existing partition. Nothing is deleted: what is
# moved aside goes to results/v3/senior_superseded/<stamp>_before_context/.
#
#   reingest_with_context.sh            (no options; stops before moving anything if an image has no context table)
set -euo pipefail
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
WT=/home/smith6jt/IO60panc2nd/Islet-Explorer-Senior-v3
DATA=$WT/core/phenocycler/data
V3=/home/smith6jt/IO60panc2nd/results/v3
EXPORT=$V3/export
PY=/home/smith6jt/IO60panc2nd/Islet-Explorer-Senior/.venv/bin/python
say() { echo "$(date +%H:%M:%S) [reingest] $*"; }
die() { say "STOPPED: $*"; exit 1; }

# ---- every image of the cohort has a context table at least as new as its cells table --------------------------
"$PY" - "$DATA/qupath_manifest.json" "$EXPORT" <<'EOF' || die "pass J has not finished for every image"
import json, re, sys
from pathlib import Path
images = json.load(open(sys.argv[1]))["images"]
export = Path(sys.argv[2])
assert len(images) == 25 and not any("6457" in i["image_id"] for i in images)
bad = []
for i in images:
    tag = re.sub(r"[^A-Za-z0-9._-]+", "_", i["image_id"])
    cells, context = export / "cells" / f"{tag}.parquet", export / "context" / f"{tag}.parquet"
    if not context.is_file():
        bad.append(f"{i['donor_id']}: no context table")
    elif context.stat().st_mtime < cells.stat().st_mtime:
        bad.append(f"{i['donor_id']}: the context table is older than the cells table")
print("\n".join(bad) if bad else "25 images: context tables present and newer than their cells tables")
sys.exit(1 if bad else 0)
EOF

SUP=$V3/senior_superseded/$(date +%Y%m%dT%H%M%S)_before_context
mkdir -p "$SUP/audit/donors"
touch "$SUP/.started"
cp -p "$DATA/qupath_manifest.json" "$SUP/qupath_manifest.json"
cp -p "$DATA/cohort-spec-v3.json" "$SUP/cohort-spec-v3.json"
for f in ingest_rejects.parquet panel_availability.json; do [[ -e $DATA/$f ]] && cp -p "$DATA/$f" "$SUP/$f"; done

# ---- 1. specification and manifest with the context tables ----------------------------------------------------
"$PY" "$HERE/make_cohort_spec.py" --senior "$WT" --export "$EXPORT" --out "$DATA/cohort-spec-v3.json" --context \
    || die "the specification could not be written (an image's files are missing)"
say "manifest: creating (hashes every cells table, context table and outline file; about 25 minutes)"
rm -f -- "$DATA/qupath_manifest.json.new"
"$HERE/senior_v3.sh" --cd core/phenocycler python -m phenocycler.pipeline manifest create \
    --spec data/cohort-spec-v3.json --out data/qupath_manifest.json.new > "$SUP/manifest_create.log" 2>&1 \
    || die "the manifest could not be created ($SUP/manifest_create.log)"
python3 - "$DATA/qupath_manifest.json.new" <<'EOF' || die "the new manifest is not what was asked for"
import json, sys
m = json.load(open(sys.argv[1]))
assert len(m["images"]) == 25, f"{len(m['images'])} images"
assert not any("6457" in i["image_id"] or str(i["donor_id"]) == "6457" for i in m["images"]), "6457 is not processed"
without = [str(i["donor_id"]) for i in m["images"] if not i.get("measurement_context")]
assert not without, f"no context table recorded for {without}"
print(f"  manifest: 25 images, each with its context table; content_id {m['content_id'][:16]}")
EOF
mv -- "$DATA/qupath_manifest.json.new" "$DATA/qupath_manifest.json"

# ---- 2 + 4. move aside what is rebuilt ------------------------------------------------------------------------
[[ -d $DATA/cells && -d $DATA/expression ]] || die "no cells / expression root in $DATA"
mv -- "$DATA/cells" "$SUP/cells"
mv -- "$DATA/expression" "$SUP/expression"
[[ -d $DATA/audit/donors/expression_availability ]] && mv -- "$DATA/audit/donors/expression_availability" "$SUP/audit/donors/expression_availability"
say "moved aside to $SUP: the cells root, the expression root and the expression audit files (geometry QC and REDSEA stay)"

# ---- the declared workflow ------------------------------------------------------------------------------------
"$HERE/run_source_workflow.sh" --redsea-workers 2 --geometry-workers 2 || die "the source workflow failed (results/v3/logs/senior_source_workflow_*.log)"

# ---- checks ---------------------------------------------------------------------------------------------------
for stage in cells geometry_qc redsea expression; do
    n=$(find "$DATA/$stage" -mindepth 1 -maxdepth 1 -type d -name 'donor_id=*' | wc -l)
    [[ $n -eq 25 ]] || die "$stage holds $n donor partitions, not 25"
done
for stage in geometry_qc redsea; do
    [[ -z $(find "$DATA/$stage" -type f -newer "$SUP/.started" | head -1) ]] || die "$stage changed during this run; it must not"
done
"$HERE/senior_v3.sh" --cd core/phenocycler python - "$SUP" <<'EOF'
import sys
from pathlib import Path
import pandas as pd
import pyarrow.parquet as pq
sup = Path(sys.argv[1])
rows = []
for d in sorted(Path("data/expression").glob("donor_id=*")):
    donor = d.name.split("=")[1]
    e = pq.read_table(next(d.glob("*.parquet")), columns=["object_id", "cell_region"]).to_pandas()
    cells = sum(pq.ParquetFile(f).metadata.num_rows for f in (Path("data/cells") / d.name).glob("*.parquet"))
    before = sum(pq.ParquetFile(f).metadata.num_rows for f in (sup / "expression" / d.name).glob("*.parquet"))
    by = e["cell_region"].value_counts().to_dict()
    rows.append({"donor": donor, "cells": cells, "expression": len(e), "expression_before": before,
                 **{k: by.get(k, 0) for k in ("core", "cluster", "single", "pLN", "tissue", "artifact", "other")}})
t = pd.DataFrame(rows)
pd.set_option("display.width", 220)
print(t.to_string(index=False))
bad = t[(t.cells != t.expression) | (t.expression != t.expression_before)]
assert bad.empty, f"row counts differ:\n{bad}"
print(f"25 donors: {int(t.expression.sum()):,} cells; core {int(t.core.sum()):,}, cluster {int(t.cluster.sum()):,}, "
      f"single {int(t.single.sum()):,}, pLN {int(t.pLN.sum()):,}, tissue {int(t.tissue.sum()):,}, other {int(t.other.sum() + t.artifact.sum()):,}")
EOF
say "done"
