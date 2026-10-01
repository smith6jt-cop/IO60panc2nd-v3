#!/usr/bin/env python3
"""islet_cell_test.py -- is this cell itself positive for a marker that made the islet object?

WHY THIS EXISTS (maintainer, 2026-10-01, on the 3-image comparison: "Too many false positives handling
segmentations around endocrine-positive cells"). An islet pixel object overshoots the islet: a bright
cell lights the cells beside it, the detector's threshold sits far below that light, and the cells of
the next ring end up more than half covered. Area alone therefore admits them. Measured on 6539 (raw
cell means): non-member cells within 5 um of an object read INS 55 at the median against 9 far from
any object, and the cells area alone admits beyond the real ones read INS 191 -- three to four times
the level that 1 % of ordinary tissue reaches (45-66). So NO level drawn from ordinary tissue can
turn them away on a hormone: 93 % of them clear it even after REDSEA correction.

So the area rule needs the second half the structure rule already has ("a vessel crossing a cell does
not make the cell a vessel"): the cell's own markers must say it. This script asks that of every cell
and writes the answer as a table; resolve_islets_to_cells.py reads it with --cell-test and stays purely
geometric.

RULES (--rule). Which one is the maintainer's choice; there is no default.

  top-marker            the structure rule's third step as worded: the marker that made the object is
                        simply the highest in the cell. The cell's highest identity marker
                        (lineages.LINEAGE_MARKERS, plus CD99) must be INS, GCG, SST or CD99.
                        [measured] on the 70-channel panel E-cadherin outranks the hormones in raw
                        counts on real endocrine cells (6533: 68 -> 21 cells in one islet).
  cd99-level            the cell's own CD99 is at or above this image's CD99 reference level.
                        [measured] a hormone-bright cell with little CD99 fails (6523panc: single
                        SST and GCG cells).
  cd99-or-hormone-top   either of the two: CD99 at or above the level, OR a hormone (INS, GCG, SST) is
                        the cell's highest identity marker.

THE REFERENCE LEVEL is the typing lane's, not a new one: the REFERENCE_QUANTILE-th percentile (99) of
the marker over a population that cannot contain the structure -- here the cells with a nucleus, in
Tissue (not the lymph node region), whose centroid lies in no islet pixel object and, when --structures
is given, in no Nerve / Blood vessel / Lymphatic object. Fitted on this image alone, never pooled. CD99
is the only marker judged against a level: it is the one of the four whose light does not reach the
neighbours enough to cross it (6539, raw: adjacent non-members read 228 at the median against a level
of about 530 and 2,274 in the object's own cells). It is not a clean cut either: endothelium carries
CD99, and a neighbour shares a CD99-bright membrane with the endocrine cell beside it.

VALUES. --values raw reads QuPath's own 'Cell: <marker>: Mean'. --values redsea reads the source
workflow's REDSEA-corrected cell means ('Cell__<marker>' in data/redsea/donor_id=<d>/data_0.parquet),
which is what the marker registry declares for these markers; the level is then fitted on the same
values. A missing measurement is unavailable: a cell with no CD99 cannot pass on CD99, and a cell with
no identity marker at all has no highest marker.

    islet_cell_test.py --rule cd99-or-hormone-top --cells cells.parquet --objects pixel_objects.geojson \\
        [--structures <export>/structures --tag <image tag>] [--values redsea --redsea data_0.parquet] \\
        --out cell_test.parquet --summary cell_test.json
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import shapely
from shapely import STRtree
from shapely.geometry import shape

SENIOR_TYPING = Path("/home/smith6jt/IO60panc2nd/Islet-Explorer-Senior/core/phenocycler/scripts/typing")
HORMONES = ("INS", "GCG", "SST")
LEVEL_MARKER = "CD99"
NOT_OBJECTS = ("Tissue", "pLN")
TISSUE_PARENT = "Annotation (Tissue)"
STRUCTURES = ("Nerve", "Blood_vessel", "Lymphatic")
RULES = {
    "top-marker": "the cell's highest identity marker (+ CD99) is INS, GCG, SST or CD99",
    "cd99-level": "the cell's CD99 is at or above the image's CD99 reference level",
    "cd99-or-hormone-top": "the cell's CD99 is at or above the image's CD99 reference level, or INS / GCG / SST is its highest identity marker",
}


def sanitize(name: str) -> str:
    """Channel name -> marker name, as phenocycler.redsea.SANITIZE does."""
    return re.sub(r"[ /\-]", "_", name)


def senior(typing_dir: Path):
    sys.path.insert(0, str(typing_dir))
    import build_crossed_fraction as bcf  # noqa: PLC0415
    import lineages  # noqa: PLC0415
    import subpanel_io  # noqa: PLC0415
    identity = sorted({m for markers in lineages.LINEAGE_MARKERS.values() for m in markers})
    return bcf, identity, float(subpanel_io.REFERENCE_QUANTILE)


def read_objects(path: Path) -> list:
    raw = json.loads(path.read_text())
    feats = raw["features"] if isinstance(raw, dict) else raw
    out = []
    for f in feats:
        if not isinstance(f, dict) or not f.get("geometry"):
            continue
        if ((f.get("properties") or {}).get("classification") or {}).get("name") in NOT_OBJECTS:
            continue
        g = shape(f["geometry"])
        if not g.is_valid:
            g = shapely.make_valid(g)
        if not g.is_empty:
            out.append(g)
    return out


def centroid_in_any(x_px: np.ndarray, y_px: np.ndarray, geoms: list, subdivide) -> np.ndarray:
    """True where the point lies within one of `geoms` (large objects are cut first, so the test stays fast)."""
    inside = np.zeros(len(x_px), dtype=bool)
    if not geoms:
        return inside
    pieces = list(subdivide(geoms))
    hit, _ = STRtree(pieces).query(shapely.points(x_px, y_px), predicate="within")
    inside[np.unique(hit)] = True
    return inside


def top_marker(values: pd.DataFrame) -> np.ndarray:
    """The highest of the columns per row; "" where none was measured (never a tie-break winner)."""
    m = values.to_numpy(dtype=float)
    none = np.isnan(m).all(axis=1)
    top = np.asarray(values.columns, dtype=object)[np.argmax(np.where(np.isnan(m), -np.inf, m), axis=1)]
    return np.where(none, "", top)


def reference_level(values: np.ndarray, quantile: float) -> float:
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        raise SystemExit(f"{LEVEL_MARKER}: the reference population carries no measurement -- refusing to fit a level")
    level = float(np.percentile(finite, quantile))
    if not np.isfinite(level) or level <= 0:
        # a level at the floor admits every measured cell while looking fitted (subpanel_io.endocrine_level)
        raise SystemExit(f"{LEVEL_MARKER}: reference level {level!r} is at the floor -- refusing to classify")
    return level


def decide(rule: str, top: np.ndarray, level_value: np.ndarray, level: float | None) -> tuple[np.ndarray, np.ndarray]:
    """(passes, why) per cell. `why` names what passed the cell, or what stood in its way."""
    n = len(top)
    top = np.asarray(top, dtype=object)
    hormone_top = np.isin(top, HORMONES)
    measured = np.isfinite(level_value)
    clears = (measured & (level_value >= level)) if level is not None else np.zeros(n, dtype=bool)
    top_text = np.array(["no identity marker measured" if t == "" else f"highest marker {t}" for t in top], dtype=object)
    if rule == "top-marker":
        return hormone_top | (top == LEVEL_MARKER), top_text
    if rule not in ("cd99-level", "cd99-or-hormone-top"):
        raise SystemExit(f"unknown rule: {rule}")
    by_top = hormone_top if rule == "cd99-or-hormone-top" else np.zeros(n, dtype=bool)
    below = np.where(measured, f"{LEVEL_MARKER} below the level", f"{LEVEL_MARKER} not measured").astype(object)
    why = np.where(clears, f"{LEVEL_MARKER} at or above the level", np.where(by_top, top_text, below + "; " + top_text))
    return clears | by_top, why.astype(object)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--rule", choices=sorted(RULES), required=True, help="the maintainer's choice; there is no default")
    ap.add_argument("--cells", type=Path, required=True, help="the image's cells Parquet (QuPath measurements)")
    ap.add_argument("--objects", type=Path, required=True, help="the detector's islet pixel objects (GeoJSON)")
    ap.add_argument("--structures", type=Path, default=None, help="<export>/structures; with --tag, their cells leave the reference population")
    ap.add_argument("--tag", default=None)
    ap.add_argument("--values", choices=("raw", "redsea"), default="raw")
    ap.add_argument("--redsea", type=Path, default=None, help="data/redsea/donor_id=<d>/data_0.parquet (with --values redsea)")
    ap.add_argument("--pixel-um", type=float, default=0.5078, help="the calibration the export's um centroids were written in")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--summary", type=Path, required=True)
    ap.add_argument("--senior-typing", type=Path, default=SENIOR_TYPING)
    args = ap.parse_args()
    if (args.values == "redsea") != (args.redsea is not None):
        raise SystemExit("--values redsea and --redsea <partition> go together")
    if (args.structures is None) != (args.tag is None):
        raise SystemExit("--structures and --tag go together")
    bcf, identity, quantile = senior(args.senior_typing)
    wanted = [*identity, LEVEL_MARKER]

    schema = pq.ParquetFile(args.cells).schema_arrow.names
    raw_column = {sanitize(c[len("Cell: "):-len(": Mean")]): c for c in schema if c.startswith("Cell: ") and c.endswith(": Mean")}
    base = pq.read_table(args.cells, columns=["Object ID", "Parent", "Centroid X µm", "Centroid Y µm", "Nucleus: Area µm^2"]
                         + ([raw_column[m] for m in wanted if m in raw_column] if args.values == "raw" else [])).to_pandas()
    if base["Object ID"].duplicated().any():
        raise SystemExit(f"{args.cells}: duplicate Object ID")
    base = base.set_index("Object ID")
    if args.values == "raw":
        have = [m for m in wanted if m in raw_column]
        values = base[[raw_column[m] for m in have]].rename(columns={raw_column[m]: m for m in have})
    else:
        red = pq.read_table(args.redsea).to_pandas().set_index("object_id")
        have = [m for m in wanted if f"Cell__{m}" in red.columns and red[f"Cell__{m}"].notna().any()]
        values = red[[f"Cell__{m}" for m in have]].rename(columns=lambda c: c.split("__", 1)[1]).reindex(base.index)
    if LEVEL_MARKER not in have and args.rule != "top-marker":
        raise SystemExit(f"{LEVEL_MARKER} is not in this image's measurements: rule {args.rule} cannot be asked")
    not_in_panel = [m for m in wanted if m not in have]

    x = base["Centroid X µm"].to_numpy(dtype=float) / args.pixel_um
    y = base["Centroid Y µm"].to_numpy(dtype=float) / args.pixel_um
    nucleated = base["Nucleus: Area µm^2"].notna().to_numpy()
    in_tissue = (base["Parent"] == TISSUE_PARENT).to_numpy()
    in_islet_object = centroid_in_any(x, y, read_objects(args.objects), bcf.subdivide)
    in_structure = np.zeros(len(base), dtype=bool)
    structure_files = []
    if args.structures is not None:
        for name in STRUCTURES:
            f = args.structures / name / f"{name}__{args.tag}.geojson"
            if not f.exists():
                raise SystemExit(f"missing structure objects: {f}")
            structure_files.append(str(f))
            in_structure |= centroid_in_any(x, y, read_objects(f), bcf.subdivide)
    reference = nucleated & in_tissue & ~in_islet_object & ~in_structure

    level = None
    if LEVEL_MARKER in have:
        level = reference_level(values[LEVEL_MARKER].to_numpy(dtype=float)[reference], quantile)
    top = top_marker(values[have])
    level_value = values[LEVEL_MARKER].to_numpy(dtype=float) if LEVEL_MARKER in have else np.full(len(base), np.nan)
    passes, why = decide(args.rule, top, level_value, level)

    table = pa.table({"object_id": base.index.to_numpy(dtype=object), "passes": passes, "why": why.astype(object),
                      "top": top.astype(object), LEVEL_MARKER: level_value})
    meta = {"rule": args.rule, "rule_text": RULES[args.rule], "values": args.values,
            "level_marker": LEVEL_MARKER, "level": "" if level is None else repr(level)}
    table = table.replace_schema_metadata({k: str(v) for k, v in meta.items()})
    args.out.parent.mkdir(parents=True, exist_ok=True)
    pq.write_table(table, args.out)
    summary = {
        **meta, "level": level, "reference_quantile": quantile, "cells_file": str(args.cells), "objects_file": str(args.objects),
        "redsea_file": str(args.redsea) if args.redsea else None, "structure_files": structure_files,
        "cells": int(len(base)), "cells_with_nucleus": int(nucleated.sum()), "cells_in_tissue": int(in_tissue.sum()),
        "reference_population": int(reference.sum()),
        "identity_markers_used": len(have) - (LEVEL_MARKER in have), "identity_markers_declared": len(identity),
        "markers_not_in_this_panel": not_in_panel,
        "cells_without_values": int((top == "").sum()),
        "passing_cells_with_nucleus": int((passes & nucleated).sum()),
        "highest_marker_share_percent_all_nucleated": {k: round(100 * v, 1) for k, v in
                                                       pd.Series(top[nucleated]).value_counts(normalize=True).head(8).items()},
    }
    args.summary.write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
