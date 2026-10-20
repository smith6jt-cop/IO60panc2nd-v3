#!/usr/bin/env python3
"""check_islet_cells_1c.py -- independent re-check of one 1c run (v3/qupath/1c_islet_objects_to_cells_v3.groovy).

Everything is recomputed with shapely/GEOS (not JTS) from what QuPath exported before and after the run
(tests/islet_cells_1c_export.groovy, OUT=<before> and OUT=<after>) and from the run's own files
(WRITE_GEOJSON=<run>.geojson: <run>.geojson, .members.tsv, .pixel_objects.tsv, .summary.json):

  hash       the Tissue/pLN/Artifact annotation hash (islet_anchors.annotation_hash) is unchanged
  cells      the cell IDs are unchanged
  saved      the saved islet-class objects are exactly the run's new objects (APPLY) or the old ones (--not-applied);
             none sits under anything but the root; no detection is left under the root
  cover      every cell with coverage >= COVER_MIN by some pixel object is a member, and only those; recorded
             coverages match
  members    every object has >= 1 member cell (n_cells = members) and every member lies inside its object
  whole      every other cell that touches an object is either <= 2% inside it or >= 98% (a cell in a filled hole),
             unless its own outline overlaps a member's (an InstanSeg segmentation overlap)
  extent     nothing of an object lies farther than TOUCH_GAP_PX + 0.1 px from its members, holes aside
  touch      members of different objects are > TOUCH_GAP_PX apart; an object's members form one touching group,
             or the extra groups lie inside the main group's filled outline
  shape      objects are valid and do not overlap; a multi-part object only with parts > TOUCH_GAP_PX apart,
             each holding member material
  caliper    Max diameter (2 x minimum bounding circle radius) and Min diameter (minimum caliper width) match
  class      SINGLE_RULE: 'cells' EndocrineSingle <=> n_cells == 1, else Islet >= ISLET_MIN_DIAM_UM mean caliper;
             'caliper' script 1's size rule; singles hold exactly one whole cell
  route      Endocrine route = OR of the feeding pixel objects' routes; Bridge flag = any of theirs
  summary    the run's summary agrees (class counts, dropped, merges, splits)
  expect     (--expect) the synthetic scenarios' known answers

Exit 0 when every check passes. Usage:
  check_islet_cells_1c.py --before B --after A --run R [--cover-min 0.5 --gap 1.0 --single-rule cells]
                          [--not-applied] [--expect expected.json] [--json-out report.json]
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import shapely
from shapely.geometry import MultiPolygon, Polygon, shape
from shapely.strtree import STRtree

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import islet_anchors as ia  # noqa: E402

ENDO = ("Islet", "EndocrineCluster", "EndocrineSingle")
REGIONS = ("Tissue", "pLN", "Artifact")
UM = "µm"


def load_fc(path: Path, repair: bool = False) -> list[dict]:
    """Features as dicts: id, cls, geom, m (measurements). repair=True (cells only) makes invalid outlines valid;
    objects are loaded as written, so the shape check sees them as they are."""
    if not path.exists():
        return []
    gj = json.load(open(path, encoding="utf-8"))
    out = []
    for f in gj.get("features", []):
        p = f.get("properties") or {}
        g = shape(f["geometry"]) if f.get("geometry") else None
        if repair and g is not None and not g.is_valid:
            g = shapely.make_valid(g)
        out.append({"id": f.get("id"), "cls": (p.get("classification") or {}).get("name"), "geom": g,
                    "m": p.get("measurements") or {}})
    return out


def region_hash(path: Path) -> str:
    polys = ia.read_annotation_polygons(path, REGIONS)
    return ia.annotation_hash(polys, REGIONS)


def fill_holes(g):
    if g.is_empty:
        return g
    polys = [g] if isinstance(g, Polygon) else [p for p in getattr(g, "geoms", []) if isinstance(p, Polygon)]
    return shapely.union_all([Polygon(p.exterior) for p in polys])


def min_width(g) -> float:
    """Minimum caliper width (JTS MinimumDiameter.getLength): rotating calipers over the convex hull edges."""
    hull = g.convex_hull
    if not isinstance(hull, Polygon):
        return 0.0
    pts = np.asarray(hull.exterior.coords)[:-1]
    best = math.inf
    for i in range(len(pts)):
        a, b = pts[i], pts[(i + 1) % len(pts)]
        d = b - a
        n = np.hypot(*d)
        if n == 0:
            continue
        normal = np.array([-d[1], d[0]]) / n
        best = min(best, float(np.max(np.abs((pts - a) @ normal))))
    return best


class Report:
    def __init__(self):
        self.rows = []

    def check(self, name: str, ok: bool, detail=""):
        self.rows.append({"check": name, "pass": bool(ok), "detail": detail})
        print(f"  {'PASS' if ok else 'FAIL'}  {name}: {detail}")

    @property
    def ok(self):
        return all(r["pass"] for r in self.rows)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--before", required=True, help="export prefix before 1c")
    ap.add_argument("--after", required=True, help="export prefix after 1c (the saved image)")
    ap.add_argument("--run", required=True, help="1c WRITE_GEOJSON base (without .geojson)")
    ap.add_argument("--cover-min", type=float, default=0.5)
    ap.add_argument("--gap", type=float, default=1.0, help="TOUCH_GAP_PX")
    ap.add_argument("--single-rule", choices=["cells", "caliper"], default="cells")
    ap.add_argument("--islet-min-um", type=float, default=35.0)
    ap.add_argument("--single-max-um", type=float, default=10.0)
    ap.add_argument("--pixel-um", type=float, default=0.5078)
    ap.add_argument("--not-applied", action="store_true", help="the run had APPLY=false")
    ap.add_argument("--expect", help="synthetic scenarios (islet_cells_1c_synthetic.groovy EXPECT=)")
    ap.add_argument("--json-out")
    a = ap.parse_args(argv)
    R = Report()
    px, px2 = a.pixel_um, a.pixel_um ** 2
    before, after, run = Path(a.before), Path(a.after), Path(a.run)
    summary = json.load(open(f"{run}.summary.json"))

    # ---- hash, cells, saved state
    hb, ha = region_hash(Path(f"{before}.regions.geojson")), region_hash(Path(f"{after}.regions.geojson"))
    R.check("hash", hb == ha and hb == summary.get("annotation_hash_before") == summary.get("annotation_hash_after"),
            f"before {hb}, after {ha}, run {summary.get('annotation_hash_before')} -> {summary.get('annotation_hash_after')}")
    ids_b = Path(f"{before}.cell_ids.txt").read_text().split()
    ids_a = Path(f"{after}.cell_ids.txt").read_text().split()
    R.check("cells", ids_b == ids_a, f"{len(ids_b)} cells before, {len(ids_a)} after")
    pix = [f for f in load_fc(Path(f"{before}.islets.geojson")) if f["cls"] in ENDO]
    saved = [f for f in load_fc(Path(f"{after}.islets.geojson")) if f["cls"] in ENDO]
    new = [f for f in load_fc(Path(f"{run}.geojson")) if f["cls"] in ENDO]
    info_a = json.load(open(f"{after}.info.json"))
    if a.not_applied:
        R.check("saved", sorted(f["id"] for f in saved) == sorted(f["id"] for f in pix),
                f"APPLY=false: {len(saved)} saved objects = the {len(pix)} pixel objects")
    else:
        same = sorted(f["id"] for f in saved) == sorted(f["id"] for f in new)
        keys = ("n_cells", "n_pixel_objects", "cover_min", "min_cell_cover", "Endocrine route")
        meas_same = same and all(all(s["m"].get(k) == n["m"].get(k) for k in keys) and s["cls"] == n["cls"]
                                 for s, n in zip(sorted(saved, key=lambda f: f["id"]), sorted(new, key=lambda f: f["id"])))
        R.check("saved", same and meas_same and not ({f["id"] for f in pix} & {f["id"] for f in saved})
                and set(info_a["islet_parents"]) <= {"ROOT"} and info_a["detections_under_root"] == 0,
                f"{len(saved)} saved = {len(new)} built, none of the {len(pix)} pixel objects left; parents "
                f"{info_a['islet_parents']}; detections under root {info_a['detections_under_root']}")

    # ---- cells (union of both exports) and membership
    cells = {}
    for p in (Path(f"{before}.cells.geojson"), Path(f"{after}.cells.geojson")):
        for f in load_fc(p, repair=True):
            cells[f["id"]] = f["geom"]
    cid = list(cells)
    cgeom = np.array([cells[i] for i in cid], dtype=object)
    carea = shapely.area(cgeom)
    cidx = {c: i for i, c in enumerate(cid)}
    ctree = STRtree(cgeom)
    members = defaultdict(list)       # object id -> cell indices
    rec_cov, kept_by = {}, {}
    with open(f"{run}.members.tsv", encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            members[r["object_id"]].append(cidx[r["cell_id"]])
            rec_cov[r["cell_id"]] = float(r["coverage"])
            kept_by[r["cell_id"]] = r["kept_by"].split(";")
    member_of = {i: oid for oid, lst in members.items() for i in lst}

    # ---- cover: recompute every (cell, pixel object) coverage
    pgeom = np.array([f["geom"] for f in pix], dtype=object)
    best = np.zeros(len(cid))
    # QuPath keeps ROI vertices as floats at 0.01 px and the GeoJSON export writes 2 decimals, so coverages
    # recomputed here differ from 1c's in-memory ones by ~1e-5: cells within COV_TOL of the threshold are ambiguous
    COV_TOL = 1e-4
    kept_expected, ambiguous = set(), set()
    p_kept = np.zeros(len(pix), dtype=int)
    for j, pg in enumerate(pgeom):
        idx = ctree.query(pg, predicate="intersects")
        if len(idx) == 0:
            continue
        ok = carea[idx] > 0
        idx = idx[ok]
        cov = shapely.area(shapely.intersection(cgeom[idx], pg)) / carea[idx]
        best[idx] = np.maximum(best[idx], cov)
        for i, c in zip(idx, cov):
            if abs(c - a.cover_min) < COV_TOL:
                ambiguous.add(int(i))
            elif c >= a.cover_min:
                kept_expected.add(int(i))
                p_kept[j] += 1
    got = set(member_of)
    diff = (kept_expected ^ got) - ambiguous
    cov_err = max((abs(rec_cov[cid[i]] - best[i]) for i in got), default=0.0)
    R.check("cover", not diff and cov_err < COV_TOL,
            f"{len(kept_expected)} cells with coverage >= {a.cover_min} (recomputed), {len(got)} members, "
            f"{len(diff)} differ, {len(ambiguous)} at the threshold; max |recorded - recomputed| {cov_err:.2e}")

    # ---- members, whole, extent, class, caliper, route
    newg = {f["id"]: f["geom"] for f in new}
    bad_members, bad_ncells, partial, holes, extent_bad, seg_overlap = [], [], [], 0, [], 0
    single_bad, class_bad, cal_bad, route_bad, whole_in = [], [], [], [], {}
    pix_by_id = {f["id"]: f for f in pix}
    for f in new:
        o, oid = f["geom"], f["id"]
        mem = members.get(oid, [])
        if len(mem) < 1 or f["m"].get("n_cells") != len(mem):
            bad_ncells.append((oid, f["m"].get("n_cells"), len(mem)))
        if mem:
            outside = shapely.area(shapely.difference(cgeom[mem], o))
            if np.any(outside > 1e-6 * carea[mem] + 1e-6):
                bad_members.append(oid)
        idx = ctree.query(o, predicate="intersects")
        frac = shapely.area(shapely.intersection(cgeom[idx], o)) / np.where(carea[idx] > 0, carea[idx], 1)
        n_whole = 0
        for i, fr in zip(idx, frac):
            if member_of.get(int(i)) == oid:
                n_whole += 1
                continue
            if fr >= 0.98:
                holes += 1
                n_whole += 1
            elif fr > 0.02:
                # allowed only where this cell's outline overlaps a member's (an InstanSeg segmentation overlap)
                if any(shapely.area(shapely.intersection(cgeom[k], cgeom[i])) > 1e-6 for k in mem if k != i):
                    seg_overlap += 1
                else:
                    partial.append((oid, cid[i], round(float(fr), 4)))
        whole_in[oid] = n_whole
        if mem:
            reach = fill_holes(shapely.union_all(cgeom[mem]).buffer(a.gap + 0.1))
            if shapely.area(shapely.difference(o, reach)) > 1e-6 * o.area + 1e-9:
                extent_bad.append(oid)
        # caliper and class
        d_max = 2 * shapely.minimum_bounding_radius(o) * px
        d_min = min_width(o) * px
        m_max, m_min = f["m"].get(f"Max diameter {UM}"), f["m"].get(f"Min diameter {UM}")
        mean_cal = f["m"].get(f"Mean caliper diameter {UM}")
        if (m_max is None or m_min is None or abs(d_max - m_max) > 1e-3 * max(1, m_max) or abs(d_min - m_min) > 1e-3 * max(1, m_min)
                or abs(mean_cal - (m_max + m_min) / 2) > 1e-9):
            cal_bad.append((oid, round(d_max, 3), m_max, round(d_min, 3), m_min))
        n = len(mem)
        if a.single_rule == "cells":
            want = "EndocrineSingle" if n == 1 else ("Islet" if mean_cal >= a.islet_min_um else "EndocrineCluster")
        else:
            want = "Islet" if mean_cal >= a.islet_min_um else ("EndocrineCluster" if mean_cal > a.single_max_um else "EndocrineSingle")
        if f["cls"] != want:
            class_bad.append((oid, f["cls"], want, n, mean_cal))
        if a.single_rule == "cells" and f["cls"] == "EndocrineSingle" and n_whole != 1:
            single_bad.append((oid, n_whole))
        feeders = sorted({p for i in mem for p in kept_by[cid[i]]})
        routes = [pix_by_id[p]["m"].get("Endocrine route") for p in feeders]
        routes = [int(r) for r in routes if r is not None]
        want_route = None if not routes else int(np.bitwise_or.reduce(routes))
        want_bridge = int(any(pix_by_id[p]["m"].get("Bridge flag") == 1 for p in feeders))
        got_route = f["m"].get("Endocrine route")
        if (want_route != (None if got_route is None else int(got_route)) or want_bridge != f["m"].get("Bridge flag")
                or f["m"].get("n_pixel_objects") != len(feeders) or f["m"].get("cover_min") != a.cover_min
                or f["m"].get("min_cell_cover", -1) < a.cover_min):
            route_bad.append((oid, got_route, want_route, f["m"].get("Bridge flag"), want_bridge, f["m"].get("n_pixel_objects"), len(feeders)))
    R.check("members", not bad_ncells and not bad_members,
            f"{len(new)} objects; n_cells != members or < 1: {bad_ncells[:3]}; members outside: {bad_members[:3]}")
    R.check("whole", not partial, f"{holes} non-member cells >= 98% inside (filled holes); {seg_overlap} partly inside through a "
                                  f"segmentation overlap with a member; partly inside otherwise: {partial[:5]}")
    R.check("extent", not extent_bad, f"objects reaching beyond gap + 0.1 px of their members: {extent_bad[:3]}")
    R.check("caliper", not cal_bad, f"mismatches: {cal_bad[:3]}")
    R.check("class", not class_bad and not single_bad,
            f"rule {a.single_rule}: wrong class {class_bad[:3]}; singles without exactly one whole cell {single_bad[:3]}")
    R.check("route", not route_bad, f"route/bridge/n_pixel_objects/cover_min mismatches: {route_bad[:3]}")

    # ---- touch: members of different objects apart; members of one object touching (or inside its outline)
    all_mem = np.array(sorted(member_of), dtype=int)
    cross, split_groups = 0, []
    if len(all_mem):
        mtree = STRtree(cgeom[all_mem])
        pairs = mtree.query(cgeom[all_mem], predicate="dwithin", distance=a.gap)
        parent = list(range(len(all_mem)))

        def find(i):
            while parent[i] != i:
                parent[i] = parent[parent[i]]
                i = parent[i]
            return i
        for i, j in zip(*pairs):
            if i >= j:
                continue
            if member_of[int(all_mem[i])] != member_of[int(all_mem[j])]:
                cross += 1
            ri, rj = find(i), find(j)
            if ri != rj:
                parent[ri] = rj
        groups = defaultdict(lambda: defaultdict(list))   # object -> root -> cell indices
        for k, i in enumerate(all_mem):
            groups[member_of[int(i)]][find(k)].append(int(i))
        for oid, gs in groups.items():
            if len(gs) == 1:
                continue
            gl = sorted(gs.values(), key=len, reverse=True)
            main_shell = fill_holes(shapely.union_all(cgeom[gl[0]]).buffer(a.gap / 2 + 0.05))
            if not all(main_shell.covers(shapely.union_all(cgeom[g])) for g in gl[1:]):
                split_groups.append((oid, [len(g) for g in gl]))
    R.check("touch", cross == 0 and not split_groups,
            f"member pairs within {a.gap} px in different objects: {cross}; objects with separate member groups "
            f"outside the main outline: {split_groups[:3]}")

    # ---- shape
    geoms = [f["geom"] for f in new]
    invalid = sum(1 for g in geoms if not g.is_valid)
    # a multi-part object is allowed only when its parts are > TOUCH_GAP_PX apart and each holds member material
    multi, multi_bad = 0, []
    for f in new:
        g = f["geom"]
        if isinstance(g, Polygon):
            continue
        multi += 1
        parts = list(g.geoms)
        memu = shapely.union_all(cgeom[members.get(f["id"], [])])
        near = any(parts[i].distance(parts[j]) <= a.gap for i in range(len(parts)) for j in range(i + 1, len(parts)))
        if near or any(p.intersection(memu).area <= 0 for p in parts):
            multi_bad.append(f["id"])
    overl = 0
    if geoms:
        gt = STRtree(geoms)
        for i, j in zip(*gt.query(geoms, predicate="intersects")):
            if i < j and geoms[i].intersection(geoms[j]).area > 1e-6:
                overl += 1
    R.check("shape", invalid == 0 and not multi_bad and overl == 0,
            f"invalid {invalid}, multi-part {multi} (bad: parts within the gap or without members {multi_bad[:3]}), overlapping pairs {overl}")

    # ---- summary consistency
    counts = {c: sum(1 for f in new if f["cls"] == c) for c in ENDO}
    dropped = int(np.sum(p_kept == 0))
    feeders_of = {f["id"]: sorted({p for i in members.get(f["id"], []) for p in kept_by[cid[i]]}) for f in new}
    merges = sum(1 for v in feeders_of.values() if len(v) >= 2)
    feeds = defaultdict(set)
    for oid, ps in feeders_of.items():
        for p in ps:
            feeds[p].add(oid)
    splits = sum(1 for v in feeds.values() if len(v) >= 2)
    s_ok = (summary["counts_after"] == counts and summary["n_dropped_halo_only"] == dropped
            and summary["merges"]["n_new_objects"] == merges and summary["splits"]["n_pixel_objects"] == splits
            and summary["n_pixel_objects"] == len(pix) and summary["n_hole_cells"] == holes)
    R.check("summary", s_ok, f"after {counts}, dropped {dropped}, merges {merges}, splits {splits}, hole cells {holes} "
                             f"(run: {summary['counts_after']}, {summary['n_dropped_halo_only']}, {summary['merges']['n_new_objects']}, "
                             f"{summary['splits']['n_pixel_objects']}, {summary['n_hole_cells']})")

    # ---- synthetic expectations
    if a.expect:
        exp = json.load(open(a.expect))
        pmap = {r["pixel_object_id"]: r for r in csv.DictReader(open(f"{run}.pixel_objects.tsv", encoding="utf-8"), delimiter="\t")}
        for sc in exp["scenarios"]:
            x0, y0, x1, y1 = sc["box"]
            inbox = lambda g: x0 <= g.centroid.x <= x1 and y0 <= g.centroid.y <= y1  # noqa: E731
            objs = sorted([f for f in new if inbox(f["geom"])], key=lambda f: (f["geom"].bounds[0], f["geom"].bounds[1]))
            got = {"n_objects": len(objs), "classes": sorted(f["cls"] for f in objs),
                   "n_cells": sorted(int(f["m"]["n_cells"]) for f in objs),
                   "n_pixel_objects": sorted(int(f["m"]["n_pixel_objects"]) for f in objs),
                   "route": sorted(int(f["m"]["Endocrine route"]) for f in objs if f["m"].get("Endocrine route") is not None),
                   "bridge": sorted(int(f["m"]["Bridge flag"]) for f in objs),
                   "border": sorted(int(f["m"]["Tissue border flag"]) for f in objs),
                   "parts": sorted(len(getattr(f["geom"], "geoms", [f["geom"]])) for f in objs)}
            pin = [f for f in pix if inbox(f["geom"])]
            got["dropped"] = sum(1 for f in pin if pmap[f["id"]]["outcome"] == "dropped")
            got["split"] = sum(1 for f in pin if pmap[f["id"]]["outcome"] == "split")
            hc = 0
            for f in objs:
                idx = ctree.query(f["geom"], predicate="intersects")
                fr = shapely.area(shapely.intersection(cgeom[idx], f["geom"])) / np.where(carea[idx] > 0, carea[idx], 1)
                hc += sum(1 for i, x in zip(idx, fr) if x >= 0.98 and member_of.get(int(i)) != f["id"])
            got["hole_cells"] = hc
            if "min_cell_cover" in sc:
                got["min_cell_cover"] = round(min(f["m"]["min_cell_cover"] for f in objs), 6) if objs else None
            # floats: the synthetic outlines are stored at QuPath's 0.01 px, so 0.6 is 0.5999 in the image
            bad = {k: (got.get(k), v) for k, v in sc.items() if k not in ("name", "box")
                   and ((got.get(k) is None or abs(got.get(k) - v) > 1e-3) if isinstance(v, float) else got.get(k) != v)}
            R.check(f"expect {sc['name']}", not bad, "ok" if not bad else f"got/expected {bad}")

    if a.json_out:
        Path(a.json_out).write_text(json.dumps({"ok": R.ok, "checks": R.rows}, indent=1))
    print(f"[check1c] {'ALL PASS' if R.ok else 'FAILED'} ({sum(r['pass'] for r in R.rows)}/{len(R.rows)} checks)")
    return 0 if R.ok else 1


if __name__ == "__main__":
    sys.exit(main())
