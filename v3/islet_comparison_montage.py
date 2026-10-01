#!/usr/bin/env python3
"""islet_comparison_montage.py -- a quick look at the islet comparison before opening QuPath.

For one image: twelve windows, each drawn twice on the same crop -- detector A's objects on the left,
detector B's on the right -- so the two can be compared on the same tissue.

  windows   2 per size stratum of A's cell-based objects (1, 2-3, 4-10, >10 cells), 2 where only A has
            an object and 2 where only B has one (chosen at random with a fixed seed, not by eye)
  picture   red INS, green CD99, blue GCG, each scaled to its own 99.5th percentile IN THE CROP: the
            colours say where a marker is brightest in that window, not how bright it is
  yellow dashed  the detector's pixel object      cyan  the cell-based object made from it
  white  cell outlines (grey dotted: no nucleus, so not a cell for any calculation)

This is a preview for orientation. The review itself is done in QuPath, where every object is there.

    islet_comparison_montage.py <image tag> <image file> <out.png>
"""
from __future__ import annotations

import json
import re
import sys

import ijson
import matplotlib
import numpy as np
import tifffile
from shapely import STRtree
from shapely.geometry import box, shape

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

REV = "/home/smith6jt/IO60panc2nd/results/v3/review/islet_comparison"
EXPORT = "/home/smith6jt/IO60panc2nd/results/v3/export"
PX = 0.5078
PAD = 40
LABEL = {"A": "A threshold detector", "B": "B Islets2"}


def parts(g):
    return list(g.geoms) if g.geom_type in ("MultiPolygon", "GeometryCollection") else [g]


def stratum(n: int) -> str:
    return "1" if n == 1 else "2-3" if n <= 3 else "4-10" if n <= 10 else ">10"


def load(tag: str, v: str):
    res = json.load(open(f"{REV}/{v}/out/resolved/{tag}.geojson"))["features"]
    raw = json.load(open(f"{REV}/{v}/out/pixel_objects/{tag}.geojson"))
    raw = raw["features"] if isinstance(raw, dict) else raw
    pix = [f for f in raw if ((f.get("properties") or {}).get("classification") or {}).get("name") not in ("Tissue", "pLN")]
    return ({"feat": res, "geom": [shape(f["geometry"]) for f in res]},
            {"feat": pix, "geom": [shape(f["geometry"]) for f in pix]})


def main() -> int:
    tag, image, out = sys.argv[1:4]
    rng = np.random.default_rng(3)
    data = {v: load(tag, v) for v in "AB"}
    trees = {v: (STRtree(data[v][0]["geom"]), STRtree(data[v][1]["geom"])) for v in "AB"}

    resA = data["A"][0]
    windows = []                                    # (why, geometry the window is centred on)
    for s in ("1", "2-3", "4-10", ">10"):
        pool = [i for i, f in enumerate(resA["feat"]) if stratum(f["properties"]["measurements"]["n_cells"]) == s]
        for i in rng.choice(pool, size=min(2, len(pool)), replace=False):
            windows.append((f"A object, {s} cell(s)", resA["geom"][i]))
    for v, other in (("A", "B"), ("B", "A")):
        res = data[v][0]
        alone = [i for i, g in enumerate(res["geom"])
                 if res["feat"][i]["properties"]["measurements"]["n_cells"] >= 2
                 and len(trees[other][1].query(g.buffer(10 / PX), predicate="intersects")) == 0]
        for i in rng.choice(alone, size=min(2, len(alone)), replace=False) if alone else []:
            windows.append((f"only {v} has an object here (2+ cells)", res["geom"][i]))

    boxes = []
    for _, g in windows:
        x0, y0, x1, y1 = g.bounds
        side = max(x1 - x0, y1 - y0) + 2 * PAD
        cx, cy = (x0 + x1) / 2, (y0 + y1) / 2
        boxes.append((int(max(0, cx - side / 2)), int(max(0, cy - side / 2)), int(cx + side / 2), int(cy + side / 2)))
    wtree = STRtree([box(*b) for b in boxes])
    cells = [[] for _ in boxes]
    with open(f"{EXPORT}/cell_geojson/{tag}.geojson", "rb") as h:
        for f in ijson.items(h, "features.item", use_float=True):
            g = shape(f["geometry"])
            for k in wtree.query(g, predicate="intersects"):
                cells[k].append((g, "nucleusGeometry" in f))

    with tifffile.TiffFile(image) as tf:
        pages = tf.series[0].levels[0].pages
        names = []
        for p in pages:
            p = p.aspage() if hasattr(p, "aspage") else p
            m = re.search(r"<Biomarker>(.*?)</Biomarker>", p.description or "")
            names.append(m.group(1) if m else "?")
        planes = {k: pages[names.index(k)].asarray() for k in ("INS", "CD99", "GCG")}

    n = len(boxes)
    rows = (n + 1) // 2
    fig, axes = plt.subplots(rows, 4, figsize=(24, 6 * rows), squeeze=False)
    for ax in axes.ravel():
        ax.axis("off")
    for k, ((why, _), (x0, y0, x1, y1)) in enumerate(zip(windows, boxes)):
        rgb = np.zeros((y1 - y0, x1 - x0, 3), dtype=np.float32)
        for c, key in enumerate(("INS", "CD99", "GCG")):
            a = planes[key][y0:y1, x0:x1].astype(np.float32)
            rgb[: a.shape[0], : a.shape[1], c] = np.clip(a / max(np.percentile(a, 99.5), 1.0), 0, 1)
        win = box(x0, y0, x1, y1)
        for j, v in enumerate("AB"):
            ax = axes[k // 2][(k % 2) * 2 + j]
            ax.imshow(rgb, extent=(x0, x1, y1, y0), interpolation="nearest")
            for g, nuc in cells[k]:
                for part in parts(g):
                    if part.geom_type == "Polygon":
                        xs, ys = part.exterior.xy
                        ax.plot(xs, ys, color="white" if nuc else "grey", lw=0.35, ls="-" if nuc else ":")
            res, pix = data[v]
            n_pix = n_res = n_cells = 0
            for i in trees[v][1].query(win, predicate="intersects"):
                n_pix += 1
                for part in parts(pix["geom"][i]):
                    if part.geom_type == "Polygon":
                        xs, ys = part.exterior.xy
                        ax.plot(xs, ys, color="yellow", lw=1.1, ls="--")
            for i in trees[v][0].query(win, predicate="intersects"):
                n_res += 1
                n_cells += res["feat"][i]["properties"]["measurements"]["n_cells"]
                for part in parts(res["geom"][i]):
                    if part.geom_type == "Polygon":
                        xs, ys = part.exterior.xy
                        ax.plot(xs, ys, color="cyan", lw=1.5)
            ax.set_title(f"{LABEL[v]}  |  window {k + 1}: {why}\n{n_pix} pixel object(s) -> {n_res} cell-based object(s), "
                         f"{n_cells} member cells   ({(x1 - x0) * PX:.0f} um wide)", fontsize=9)
            ax.set_xlim(x0, x1)
            ax.set_ylim(y1, y0)
    fig.suptitle(f"{tag.split('_Scan')[0]}: red INS, green CD99, blue GCG (each scaled to its own 99.5th percentile in the crop). "
                 "Yellow dashed = pixel object, cyan = cell-based object, white = cells (grey dotted: no nucleus)", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.985))
    fig.savefig(out, dpi=80)
    print("->", out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
