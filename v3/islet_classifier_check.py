#!/usr/bin/env python3
"""islet_classifier_check.py -- what QuPath made of the image's classifier, counted the way the fit counts.

Two things, both from the objects QuPath itself exported (1_panc_objects_perimage.groovy on a private copy):

  1. the build: islet_classifier_fit.py reports from an emulation of script 1 on the classifier's own 1 um grid. Are
     the same cells inside QuPath's objects as inside the emulated ones (by centroid)? This checks the classifier
     file, the op chain and the emulation -- not whether the classifier is right.
  2. the same counts for QuPath's objects, and for the objects script 1 makes with Islets2 when --baseline is given,
     so the two can be read side by side: cells inside, neighbours inside (`other`, at the edge), endocrine cells
     left outside.

The counts rest on the image's own cells as islet_classifier_fit.py groups them (endocrine / other / left out). They
describe the classifier against the cells that set it; the check that counts is looking at the image.

    islet_classifier_check.py --tag TAG --fit fit/TAG.json --fit-objects fit/TAG.objects.npy \\
        --qupath-objects trial/out/reviewed_islets/TAG.geojson [--baseline trial/out_islets2/reviewed_islets/TAG.geojson] \\
        --cells cells/TAG.parquet --redsea redsea/donor_id=D/data_0.parquet --structures <export>/structures --out fit/TAG.check.json
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import shapely
from shapely import STRtree

sys.path.insert(0, str(Path(__file__).resolve().parent))
import islet_cell_test as ict  # noqa: E402
import islet_classifier_fit as fit  # noqa: E402

PX0 = fit.CLASSIFIER_UM / 2.0


def inside_objects(geoms: list, x_px: np.ndarray, y_px: np.ndarray) -> np.ndarray:
    """1-based index of the object each point lies in (0 = none). An object is one exported annotation."""
    which = np.zeros(len(x_px), dtype=np.int64)
    if geoms:
        hit, obj = STRtree(geoms).query(shapely.points(x_px, y_px), predicate="within")
        which[hit] = obj + 1
    return which


def score_geoms(geoms: list, x_px: np.ndarray, y_px: np.ndarray, lab: dict, eligible: np.ndarray) -> dict:
    which = inside_objects(geoms, x_px, y_px)
    inside = eligible & (which > 0)
    per = np.bincount(which[inside], minlength=len(geoms) + 1)[1:]
    pct = lambda a, b: round(100.0 * a / b, 2) if b else None  # noqa: E731
    n_in = int(inside.sum())
    return {
        "objects": len(geoms), "objects_with_cells": int((per > 0).sum()), "objects_by_cells": fit.strata(per[per > 0]),
        "area_mm2": round(float(sum(g.area for g in geoms) * PX0 * PX0 / 1e6), 4),
        "cells_inside": n_in,
        "inside_endocrine": int((inside & lab["endocrine"]).sum()),
        "inside_other_at_the_edge": int((inside & lab["other"]).sum()),
        "inside_other_at_the_edge_percent": pct(int((inside & lab["other"]).sum()), n_in),
        "inside_other_enclosed": int((inside & lab["enclosed"]).sum()),
        "inside_left_out": int((inside & lab["left_out"]).sum()),
        "endocrine_outside": int((~inside & lab["endocrine"]).sum()),
        "endocrine_outside_percent": pct(int((~inside & lab["endocrine"]).sum()), int(lab["endocrine"].sum())),
        "confirmed_outside": int((~inside & lab["confirmed"]).sum()),
        "confirmed_outside_percent": pct(int((~inside & lab["confirmed"]).sum()), int(lab["confirmed"].sum())),
        "single_cell_objects": fit.single_cell_objects(which, inside, per, lab),
    }, inside


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tag", required=True)
    ap.add_argument("--fit", type=Path, required=True)
    ap.add_argument("--fit-objects", type=Path, required=True)
    ap.add_argument("--qupath-objects", type=Path, required=True)
    ap.add_argument("--baseline", type=Path, default=None, help="the objects script 1 makes with Islets2 on the same image")
    ap.add_argument("--cells", type=Path, required=True)
    ap.add_argument("--redsea", type=Path, required=True)
    ap.add_argument("--structures", type=Path, default=None)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    report = json.loads(args.fit.read_text())
    base, values, x_um, y_um, eligible, in_structure, other_markers, quantile, _ = fit.load_cells(
        args.cells, args.redsea, args.structures, args.tag, ict.SENIOR_TYPING)
    lab = fit.label_cells(values, x_um, y_um, eligible, in_structure, other_markers, quantile)
    x_px, y_px = x_um / PX0, y_um / PX0

    new, in_qupath = score_geoms(ict.read_objects(args.qupath_objects), x_px, y_px, lab, eligible)
    emulated = np.load(args.fit_objects)
    h, w = emulated.shape
    cx = np.clip((x_um / fit.CLASSIFIER_UM).astype(int), 0, w - 1); cy = np.clip((y_um / fit.CLASSIFIER_UM).astype(int), 0, h - 1)
    in_emulated = eligible & (emulated[cy, cx] > 0)
    both = int((in_qupath & in_emulated).sum())
    out = {
        "tag": args.tag, "thresholds": report["thresholds"],
        "cells": {"with_nucleus_in_tissue": int(eligible.sum()), "endocrine": int(lab["endocrine"].sum()), "confirmed": int(lab["confirmed"].sum()),
                  "other": int(lab["other"].sum()), "other_enclosed": int(lab["enclosed"].sum()), "left_out": int(lab["left_out"].sum())},
        "qupath": new,
        "emulation": {"cells_inside": int(in_emulated.sum()), "cells_in_both": both, "only_qupath": int((in_qupath & ~in_emulated).sum()),
                      "only_emulated": int((in_emulated & ~in_qupath).sum()),
                      "agreement_percent": round(100.0 * both / max(int((in_qupath | in_emulated).sum()), 1), 2),
                      "area_mm2": report["after_script_1"]["area_mm2"]},
    }
    if args.baseline is not None:
        out["islets2"], _ = score_geoms(ict.read_objects(args.baseline), x_px, y_px, lab, eligible)
    args.out.write_text(json.dumps(out, indent=2))
    print(json.dumps(out, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
