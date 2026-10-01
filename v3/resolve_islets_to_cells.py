#!/usr/bin/env python3
"""resolve_islets_to_cells.py -- turn pixel-classifier islet objects into objects made of whole cells.

The maintainer's complaint (2026-09-29): an islet object is a thresholded pixel region, not resolved
to the cell segmentation, so one bright cell's halo over parts of its neighbours is read as a
multi-cell object. The rule applied here is the one the Senior structure call already uses and
tests, not a new one:

    a cell belongs to an object when the MAJORITY of the cell's area lies inside it
    (area fraction > subpanel.MAJORITY_AREA, strict), and only a cell with a nucleus is a cell.

Nothing here detects anything. Each pixel object keeps its identity:
    members(object)  = cells with a nucleus and > half their area inside that object
    resolved object  = the union of its members' own outlines, as one region: the one-pixel line the
                       segmentation leaves between two touching cells is closed (--seam-px), and holes
                       enclosed by members are closed, so the outline is the outer boundary of the
                       cell set
    an object with no member is dropped (it was halo, debris, or cells without a nucleus)
Nothing wider than the seam is bridged and no two pixel objects are merged. The area arithmetic is
Senior's `build_crossed_fraction.covered_fraction` / `subdivide`, imported, not copied.

The seam [measured 2026-10-01, 6523panc]: InstanSeg outlines are rectilinear on the model's 0.5 um
grid (0.985 image px), and neighbouring cells are separated by one grid step or meet at corners --
93-95 % of the gaps between the parts of a plain union are <= 1 px. Without closing it, every islet
is a lace of separate cells and no enclosed hole can ever close. `--seam-px 0` gives the plain union.

Size classes follow Drotar 2025 as published, measured on the resolved outline:
    mean diameter = (minimum caliper width + maximum diameter) / 2, as QuPath defines them
    <= 10 um EndocrineSingle | < 35 um EndocrineCluster | >= 35 um Islet
`n_cells` is written beside the class so "single = exactly one cell" can be read off as well.

    resolve_islets_to_cells.py --objects pixel_objects.geojson --cells cell_outlines.geojson \
        --out resolved.geojson --members members.tsv --summary summary.json
"""
from __future__ import annotations

import argparse
import json
import math
import sys
import time
from collections import Counter
from pathlib import Path

import ijson
import numpy as np
import shapely
from shapely import STRtree
from shapely.geometry import MultiPolygon, Polygon, mapping, shape
from shapely.ops import unary_union

SENIOR_TYPING = Path("/home/smith6jt/IO60panc2nd/Islet-Explorer-Senior/core/phenocycler/scripts/typing")
ENDOCRINE = ("Islet", "EndocrineCluster", "EndocrineSingle")
NOT_OBJECTS = ("Tissue", "pLN")
SINGLE_MAX_UM, ISLET_MIN_UM = 10.0, 35.0          # Drotar 2025, p7
COLORS = {"Islet": [200, 0, 0], "EndocrineCluster": [255, 140, 0], "EndocrineSingle": [0, 200, 200]}


def senior(typing_dir: Path):
    """Senior's own area arithmetic and majority constant, imported so the rule cannot drift."""
    sys.path.insert(0, str(typing_dir))
    import build_crossed_fraction as bcf  # noqa: PLC0415
    import subpanel  # noqa: PLC0415
    return bcf, float(subpanel.MAJORITY_AREA)


def class_of(feature: dict) -> str:
    return ((feature.get("properties") or {}).get("classification") or {}).get("name") or "Unclassified"


def min_caliper_width(geom) -> float:
    """Minimum caliper width of the convex hull: the smallest distance between two parallel support
    lines, attained with one line flush to a hull edge (QuPath's 'Min diameter')."""
    hull = geom.convex_hull
    if hull.geom_type != "Polygon":
        return 0.0
    pts = np.asarray(hull.exterior.coords)[:-1]
    if len(pts) < 3:
        return 0.0
    edges = np.roll(pts, -1, axis=0) - pts
    lengths = np.hypot(edges[:, 0], edges[:, 1])
    keep = lengths > 0
    normals = np.stack([-edges[keep, 1], edges[keep, 0]], axis=1) / lengths[keep, None]
    # |.|: the hull ring may run either way (GEOS returns it clockwise), so the normal may point out
    widths = np.abs(normals @ pts.T - np.sum(normals * pts[keep], axis=1)[:, None]).max(axis=1)
    return float(widths.min())


def join_cells(cells: list, seam_px: float):
    """The union of `cells` as a region. With seam_px > 0, gaps no wider than seam_px are closed
    (a square closing: grow each outline by half the seam, unite, shrink by the same amount), which
    joins cells separated only by the segmentation's one-pixel boundary line and adds nothing at the
    outer boundary. seam_px = 0 is the plain union."""
    if seam_px <= 0 or len(cells) < 2:
        return unary_union(cells)
    r = 0.51 * seam_px
    grown = unary_union([c.buffer(r, join_style="mitre") for c in cells])
    closed = grown.buffer(-r, join_style="mitre")
    # the shrink can only lose what the grow added; never return less than the cells themselves
    return unary_union([closed, *cells])


def outer_boundary(geom):
    """`geom` with every hole enclosed by it closed: the outer boundary of each part."""
    parts = list(geom.geoms) if geom.geom_type == "MultiPolygon" else [geom]
    filled = [Polygon(p.exterior) for p in parts if p.geom_type == "Polygon" and not p.is_empty]
    if not filled:
        return geom
    merged = unary_union(filled)
    return merged


def size_class(mean_diameter_um: float) -> str:
    if mean_diameter_um >= ISLET_MIN_UM:
        return "Islet"
    return "EndocrineCluster" if mean_diameter_um > SINGLE_MAX_UM else "EndocrineSingle"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--objects", type=Path, required=True, help="annotations GeoJSON holding the pixel objects")
    ap.add_argument("--cells", type=Path, required=True, help="cell-outline GeoJSON (feature id = QuPath UUID)")
    ap.add_argument("--out", type=Path, required=True, help="resolved objects, QuPath-importable GeoJSON")
    ap.add_argument("--members", type=Path, required=True, help="TSV: cell id, object id, area fraction")
    ap.add_argument("--summary", type=Path, required=True)
    ap.add_argument("--pixel-um", type=float, default=0.5078)
    ap.add_argument("--seam-px", type=float, default=1.0,
                    help="close gaps between member cells no wider than this (image px); 0 = plain union")
    ap.add_argument("--senior-typing", type=Path, default=SENIOR_TYPING)
    args = ap.parse_args()
    started = time.time()
    bcf, majority = senior(args.senior_typing)
    px = args.pixel_um

    # ---- pixel objects: every annotation that is not a region ----
    raw = json.loads(args.objects.read_text())
    feats = raw["features"] if isinstance(raw, dict) else raw
    objects = []
    for f in feats:
        if not isinstance(f, dict) or class_of(f) in NOT_OBJECTS:
            continue
        g = shape(f["geometry"])
        if not g.is_valid:
            g = shapely.make_valid(g)
        if g.is_empty:
            continue
        objects.append({"id": f.get("id"), "class": class_of(f), "geom": g,
                        "name": (f.get("properties") or {}).get("name"),
                        "measurements": (f.get("properties") or {}).get("measurements") or {}})
    if not objects:
        raise SystemExit(f"no pixel objects in {args.objects}")
    if len({o["id"] for o in objects}) != len(objects) or any(o["id"] is None for o in objects):
        raise SystemExit("pixel objects must carry unique ids")

    # pieces of an oversized object share its index, so its areas add exactly (bcf.subdivide)
    pieces, owner = [], []
    for k, o in enumerate(objects):
        for piece in bcf.subdivide([o["geom"]]):
            pieces.append(piece)
            owner.append(k)
    tree = STRtree(pieces)
    owner = np.asarray(owner)

    # ---- every cell against the objects it touches ----
    members: dict[int, list] = {}
    rows = []
    n_cells = n_touching = n_majority_no_nucleus = 0
    with open(args.cells, "rb") as handle:
        for f in ijson.items(handle, "features.item", use_float=True):
            n_cells += 1
            outline = shape(f["geometry"])
            hits = tree.query(outline, predicate="intersects")
            if len(hits) == 0:
                continue
            n_touching += 1
            cell_area = float(outline.area)
            if cell_area <= 0:
                continue
            if not outline.is_valid:
                outline = shapely.make_valid(outline)
            inside: dict[int, float] = {}
            for i in hits:
                a = float(pieces[i].intersection(outline).area)
                if a > 0:
                    inside[int(owner[i])] = inside.get(int(owner[i]), 0.0) + a
            if not inside:
                continue
            k, a = max(inside.items(), key=lambda kv: kv[1])
            frac = a / cell_area
            if not frac > majority:                    # strict, as subpanel.py:871
                continue
            if "nucleusGeometry" not in f:             # no nucleus: not a cell for any calculation
                n_majority_no_nucleus += 1
                continue
            total = bcf.covered_fraction(tree, pieces, outline, cell_area)
            members.setdefault(k, []).append(outline)
            rows.append((f["id"], objects[k]["id"], frac, total))

    # ---- resolved objects ----
    out_feats, dropped = [], Counter()
    by_class_pixel = Counter(o["class"] for o in objects)
    by_class_resolved, by_cells = Counter(), Counter()
    multipart = 0
    for k, o in enumerate(objects):
        cells = members.get(k)
        if not cells:
            dropped[o["class"]] += 1
            continue
        geom = outer_boundary(join_cells(cells, args.seam_px))
        n_parts = len(geom.geoms) if geom.geom_type == "MultiPolygon" else 1
        multipart += n_parts > 1
        max_d = 2.0 * shapely.minimum_bounding_radius(geom) * px
        min_d = min_caliper_width(geom) * px
        mean_d = (max_d + min_d) / 2.0
        cls = size_class(mean_d)
        by_class_resolved[cls] += 1
        by_cells["1 cell" if len(cells) == 1 else "2-3" if len(cells) <= 3 else "4-10" if len(cells) <= 10 else ">10"] += 1
        meas = dict(o["measurements"])
        meas.update({"n_cells": len(cells), "n_parts": n_parts,
                     "Resolved area um^2": geom.area * px * px,
                     "Pixel object area um^2": o["geom"].area * px * px,
                     "Max diameter um": max_d, "Min diameter um": min_d, "Mean diameter um": mean_d})
        out_feats.append({
            "type": "Feature", "id": o["id"], "geometry": mapping(geom),
            "properties": {"objectType": "annotation", "name": o["name"],
                           "classification": {"name": cls, "color": COLORS[cls]},
                           "measurements": meas, "pixel_object_class": o["class"]}})
    for f in out_feats:
        if f["properties"]["name"] is None:
            del f["properties"]["name"]

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps({"type": "FeatureCollection", "features": out_feats}))
    with open(args.members, "w") as h:
        h.write("cell_id\tobject_id\tarea_fraction_in_object\tarea_fraction_in_any_object\n")
        for cid, oid, frac, total in rows:
            h.write(f"{cid}\t{oid}\t{frac:.6f}\t{total:.6f}\n")
    summary = {
        "objects_file": str(args.objects), "cells_file": str(args.cells), "pixel_um": px, "seam_px": args.seam_px,
        "majority_area": majority, "rule": "cell with a nucleus and area fraction > majority_area inside one pixel object",
        "cells_in_file": n_cells, "cells_touching_an_object": n_touching, "member_cells": len(rows),
        "majority_cells_without_nucleus": n_majority_no_nucleus,
        "pixel_objects": len(objects), "pixel_objects_by_class": dict(by_class_pixel),
        "resolved_objects": len(out_feats), "resolved_by_class_drotar": dict(by_class_resolved),
        "resolved_by_cell_count": dict(by_cells), "dropped_no_member_by_pixel_class": dict(dropped),
        "resolved_multipart": multipart, "seconds": round(time.time() - started, 1),
    }
    args.summary.write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
