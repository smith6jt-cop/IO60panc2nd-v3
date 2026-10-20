#!/usr/bin/env python3
"""islet_classifier_check.py -- what QuPath made of the image's classifier, counted the way the fit counts.

From the objects QuPath itself exported on a private copy (islet_classifiers.sh trial):

  1. the build of script 1: islet_classifier_fit.py reports from an emulation of script 1 on the classifier's own
     1 um grid. Are the same cells inside QuPath's script-1 outlines as inside the emulated ones (by centroid)?
  2. the build of step 1b (outlines erode to the cell outlines; Hormone Fragments -- maintainer, 2026-10-02):
     islet_outlines_to_cells.py applies the same rule in Python to script 1's outlines and the exported cell
     outlines. Same islets, same whole cells in each?
  3. the counts, side by side: the final islets, script 1's outlines before 1b, and the objects script 1 makes
     with Islets2 when --baseline is given: cells inside, neighbours inside (`other`, at the edge), endocrine cells
     left outside. And for the Hormone Fragments: how many, and what the cell under each one is.

1 and 2 check the classifier file, the op chain, the Groovy step and the emulations -- not whether any of it is
right. The counts rest on the image's own cells as islet_classifier_fit.py groups them (endocrine / other / left
out); the check that counts is looking at the image.

    islet_classifier_check.py --tag TAG --fit fit/TAG.json --fit-objects fit/TAG.objects.npy \\
        --qupath-objects trial/out/reviewed_islets/TAG.geojson [--script1-objects trial/out/script1_outlines/TAG.geojson] \\
        [--cell-outlines <export>/cell_geojson/TAG.geojson] [--baseline trial/out_islets2/reviewed_islets/TAG.geojson] \\
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
import islet_outlines_to_cells as oc  # noqa: E402

PX0 = fit.CLASSIFIER_UM / 2.0


def read_by_class(path: Path) -> dict:
    """An islets export split by class: {class name: [geometry, ...]}; and its endocrine objects (Islet,
    EndocrineCluster, EndocrineSingle) in file order with their class, 'Whole cells' and recorded diameter."""
    from shapely.geometry import shape  # noqa: PLC0415
    raw = json.loads(path.read_text())
    feats = raw["features"] if isinstance(raw, dict) else raw
    out, whole, endocrine = {}, [], []
    for f in feats:
        if not isinstance(f, dict) or not f.get("geometry"):
            continue
        props = f.get("properties") or {}
        name = (props.get("classification") or {}).get("name")
        if name in ict.NOT_OBJECTS:
            continue
        g = shape(f["geometry"])
        if not g.is_valid:
            g = shapely.make_valid(g)
        if g.is_empty:
            continue
        out.setdefault(name, []).append(g)
        if name in oc.ENDOCRINE_CLASSES:
            meas = props.get("measurements") or {}
            whole.append(meas.get(oc.KEY_CELLS))
            endocrine.append({"class": name, "geometry": g, "whole": meas.get(oc.KEY_CELLS), "diameter": meas.get(oc.KEY_DIAMETER)})
    return {"by_class": out, "whole_cells": whole, "endocrine": endocrine}


def designation_report(endocrine: list) -> dict:
    """The classes QuPath gave against the rule applied here to the same outline and the same count of whole cells."""
    counts = {c: 0 for c in oc.ENDOCRINE_CLASSES}
    agree = unknown = 0
    worst = 0.0
    near = {"30-35 um": 0, "35-40 um": 0}
    for e in endocrine:
        counts[e["class"]] += 1
        if e["whole"] is None:
            unknown += 1
            continue
        across = oc.mean_caliper(e["geometry"]) * PX0
        agree += oc.designation(float(e["whole"]), across) == e["class"]
        if e["diameter"] is not None:
            worst = max(worst, abs(float(e["diameter"]) - across))
        if float(e["whole"]) >= 2:
            near["30-35 um"] += 30.0 <= across < 35.0
            near["35-40 um"] += 35.0 <= across < 40.0
    return {"by_class": counts, "objects": len(endocrine), "same_class_by_the_rule_here": int(agree), "without_whole_cells": int(unknown),
            "largest_difference_in_diameter_um": round(worst, 6), "two_or_more_cells_near_the_line": {k: int(v) for k, v in near.items()},
            "islet_min_diameter_um": oc.ISLET_MIN_DIAMETER_UM}


def fragments_report(fragments: list, script1: list | None, x_px, y_px, lab: dict, eligible) -> dict:
    """What lies under the Hormone Fragments: the cell whose centre is in the outline script 1 drew there."""
    areas = np.array([g.area for g in fragments]) * PX0 * PX0
    out = {"hormone_fragments": len(fragments), "area_um2_p50": round(float(np.median(areas)), 1) if len(areas) else None,
           "area_um2_p90": round(float(np.percentile(areas, 90)), 1) if len(areas) else None}
    if script1 and fragments:
        s1 = np.asarray(script1, dtype=object)
        frag = np.asarray(fragments, dtype=object)
        hit_f, hit_s = STRtree(s1).query(frag, predicate="intersects")
        # the script-1 outline each fragment came from: the one it lies in (a neighbour may share an edge with it)
        overlap = shapely.area(shapely.intersection(frag[hit_f], s1[hit_s]))
        own = np.full(len(frag), -1); order = np.argsort(overlap)
        own[hit_f[order]] = hit_s[order]
        source = np.unique(own[own >= 0])
        ci, _ = STRtree(s1[source]).query(shapely.points(x_px, y_px), predicate="within")
        under = np.zeros(len(x_px), dtype=bool); under[np.unique(ci)] = True
        under &= eligible
        out["cells_with_centre_in_a_fragments_old_outline"] = {
            "all": int(under.sum()), "confirmed_endocrine": int((under & lab["confirmed"]).sum()),
            "endocrine_by_cd99": int((under & lab["endocrine"] & ~lab["confirmed"]).sum()),
            "other": int((under & (lab["other"] | lab["enclosed"])).sum()), "left_out": int((under & lab["left_out"] & ~lab["enclosed"]).sum())}
    return out


def reference_agreement(script1: list, cell_outlines: Path, qupath_islets: list, qupath_whole: list, n_fragments: int) -> dict:
    """islet_outlines_to_cells.py on script 1's outlines against what QuPath's step 1b made of them."""
    s1 = np.asarray(script1, dtype=object)
    _, cells, nucleus = oc.cells_near(cell_outlines, s1)
    ref = oc.erode_to_cells(s1, cells, nucleus)
    ref_geoms = np.asarray([i["geometry"] for i in ref["islets"]], dtype=object)
    ref_cells = np.array([len(i["cells"]) for i in ref["islets"]])
    q = np.asarray(qupath_islets, dtype=object)
    out = {"islets_qupath": len(q), "islets_reference": len(ref_geoms), "whole_cells_qupath": int(np.nansum(np.array(qupath_whole, dtype=float))),
           "whole_cells_reference": int(ref_cells.sum()), "fragments_reference_before_marker_outline": len(ref["fragments"]),
           "fragments_qupath": n_fragments}
    if len(q) and len(ref_geoms):
        a, b = STRtree(ref_geoms).query(q, predicate="intersects")
        inter = shapely.area(shapely.intersection(q[a], ref_geoms[b]))
        best_area = np.zeros(len(q)); best = np.full(len(q), -1)
        order = np.argsort(inter)
        best_area[a[order]] = inter[order]; best[a[order]] = b[order]
        union = np.array([shapely.area(shapely.union(q[i], ref_geoms[best[i]])) if best[i] >= 0 else np.inf for i in range(len(q))])
        iou = best_area / union
        out.update({"iou_min": round(float(iou.min()), 4), "iou_p1": round(float(np.percentile(iou, 1)), 4), "iou_median": round(float(np.median(iou)), 6),
                    "islets_with_iou_over_0.99": int((iou > 0.99).sum())})
    return out


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
    ap.add_argument("--qupath-objects", type=Path, required=True, help="the final export: Islets on cell outlines + Hormone Fragments")
    ap.add_argument("--script1-objects", type=Path, default=None, help="script 1's outlines, exported before step 1b")
    ap.add_argument("--cell-outlines", type=Path, default=None, help="<export>/cell_geojson/TAG.geojson: also run the Python reference of 1b")
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

    final = read_by_class(args.qupath_objects)
    islets = [e["geometry"] for e in final["endocrine"]]          # Islet, EndocrineCluster and EndocrineSingle alike
    fragments = final["by_class"].get(oc.CLASS_FRAGMENT, [])
    new, _ = score_geoms(islets, x_px, y_px, lab, eligible)
    # script 1 makes Islets; anything else in that export (Hormone Fragments of an earlier run of 1b, which script 1
    # does not remove) is not one of its outlines
    script1 = read_by_class(args.script1_objects)["by_class"].get(oc.CLASS_ISLET, []) if args.script1_objects is not None else None
    # the emulation of script 1 is compared with script 1's own outlines (before 1b), when they were kept
    before, in_script1 = score_geoms(script1 if script1 is not None else islets + fragments, x_px, y_px, lab, eligible)
    emulated = np.load(args.fit_objects)
    h, w = emulated.shape
    cx = np.clip((x_um / fit.CLASSIFIER_UM).astype(int), 0, w - 1); cy = np.clip((y_um / fit.CLASSIFIER_UM).astype(int), 0, h - 1)
    in_emulated = eligible & (emulated[cy, cx] > 0)
    both = int((in_script1 & in_emulated).sum())
    out = {
        "tag": args.tag, "thresholds": report["thresholds"],
        "cells": {"with_nucleus_in_tissue": int(eligible.sum()), "endocrine": int(lab["endocrine"].sum()), "confirmed": int(lab["confirmed"].sum()),
                  "other": int(lab["other"].sum()), "other_enclosed": int(lab["enclosed"].sum()), "left_out": int(lab["left_out"].sum())},
        "qupath": new,
        "designation": designation_report(final["endocrine"]),
        "fragments": fragments_report(fragments, script1, x_px, y_px, lab, eligible),
        "emulation": {"cells_inside": int(in_emulated.sum()), "cells_in_both": both, "only_qupath": int((in_script1 & ~in_emulated).sum()),
                      "only_emulated": int((in_emulated & ~in_script1).sum()),
                      "agreement_percent": round(100.0 * both / max(int((in_script1 | in_emulated).sum()), 1), 2),
                      "area_mm2": report["after_script_1"]["area_mm2"]},
    }
    if script1 is not None:
        out["script_1"] = before
        if args.cell_outlines is not None:
            out["step_1b_reference"] = reference_agreement(script1, args.cell_outlines, islets, final["whole_cells"], len(fragments))
    if args.baseline is not None:
        out["islets2"], _ = score_geoms(ict.read_objects(args.baseline), x_px, y_px, lab, eligible)
    args.out.write_text(json.dumps(out, indent=2))
    print(json.dumps(out, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
