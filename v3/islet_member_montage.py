#!/usr/bin/env python3
"""islet_member_montage.py -- which cells an islet object takes in under each candidate member rule.

For one image and one detector: windows chosen at random (fixed seed, not by eye) from the objects of the
FIRST rule, the same number per size stratum of that rule's objects (1, 2-3, 4-10, >10 cells).  Each window
is one row, drawn once per rule on the same crop, and once more for CD99:

  1..n  hormones (red INS, green GCG, blue SST) with that rule's object; from the second rule on, a cell the
        first rule had and this rule drops is outlined in orange and labelled with its highest marker
  last  CD99 alone (grey) with the last rule's object

Brightness is NOT stretched per crop.  For every channel the black point is that image's background anchor
(islet_anchors.csv: the level the detector itself measured over Tissue) and the white point is one number
per image and channel (the 90th percentile of the cell means of the cells of --levels-from, either detector,
whose highest marker is that channel; CD99: all of them).  Linear, gamma 1.  So a cell at background is
black in every window, and a window with nothing above background is black.

  yellow dashed  the detector's pixel object       cyan  the cell-based object
  white  cell outlines (grey dotted: no nucleus, so not a cell for any calculation)

    islet_member_montage.py --tag <image tag> --image <qptiff> --detector A|B --top-marker <cell_test.parquet>
                            --sets "area only=resolved" "CD99 level or hormone highest=resolved_cd99-or-hormone-top" [...]
                            [--per-stratum 3] [--rows-per-figure 4] --out fig.png      (writes fig_1.png, fig_2.png, ...)
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import ijson
import matplotlib
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import tifffile
from shapely import STRtree
from shapely.geometry import box, shape

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

REV = Path("/home/smith6jt/IO60panc2nd/results/v3/review/islet_comparison")
EXPORT = Path("/home/smith6jt/IO60panc2nd/results/v3/export")
PX = 0.5078
PAD = 30
HORMONES = ("INS", "GCG", "SST")
CHANNELS = (*HORMONES, "CD99")
STRATA = ("1", "2-3", "4-10", ">10")
SHORT = {"E_cadherin": "ECad", "Ker8_18": "K8/18", "Pan_Cytokeratin": "PanCK", "Vimentin": "Vim", "Caveolin": "Cav",
         "Collagen_IV": "ColIV", "Podoplanin": "PDPN"}


def parts(g):
    return list(g.geoms) if g.geom_type in ("MultiPolygon", "GeometryCollection") else [g]


def stratum(n: int) -> str:
    return "1" if n == 1 else "2-3" if n <= 3 else "4-10" if n <= 10 else ">10"


def objects(path: Path):
    feats = json.loads(path.read_text())
    feats = feats["features"] if isinstance(feats, dict) else feats
    feats = [f for f in feats if ((f.get("properties") or {}).get("classification") or {}).get("name") not in ("Tissue", "pLN")]
    return feats, [shape(f["geometry"]) for f in feats]


def display_levels(tag: str, after: str, top: pd.Series) -> tuple[dict, dict, dict]:
    """Black and white point per channel, one pair per image; see the module text."""
    anchors = pd.read_csv(EXPORT / "islet_anchors.csv")
    black = anchors[anchors["tag"] == tag].set_index("channel")["anchor"].to_dict()
    missing = [c for c in CHANNELS if c not in black]
    if missing:
        raise SystemExit(f"no background anchor for {missing} on {tag}: the picture cannot be windowed on a measured background")
    kept = set()
    for v in "AB":
        f = REV / v / "out" / after / f"{tag}.members.tsv"
        if f.exists():
            kept |= set(pd.read_csv(f, sep="\t", usecols=["cell_id"])["cell_id"])
    means = pq.read_table(EXPORT / "cells" / f"{tag}.parquet",
                          columns=["Object ID"] + [f"Cell: {c}: Mean" for c in CHANNELS]).to_pandas().set_index("Object ID")
    means = means.loc[means.index.intersection(list(kept))]
    t = top.reindex(means.index)
    white, n_used = {}, {}
    for c in CHANNELS:
        pool = means[f"Cell: {c}: Mean"] if c == "CD99" else means.loc[t == c, f"Cell: {c}: Mean"]
        pool = pool.dropna()
        if len(pool) < 50:
            raise SystemExit(f"{tag}: only {len(pool)} kept cells to set the white point of {c}")
        white[c], n_used[c] = float(np.percentile(pool, 90)), int(len(pool))
        if white[c] <= black[c]:
            raise SystemExit(f"{tag}: white point of {c} ({white[c]:.0f}) is not above its background anchor ({black[c]:.0f})")
    return black, white, n_used


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tag", required=True)
    ap.add_argument("--image", required=True)
    ap.add_argument("--detector", choices=("A", "B"), required=True)
    ap.add_argument("--top-marker", type=Path, required=True)
    ap.add_argument("--sets", nargs="+", required=True,
                    help="label=subdirectory of <detector>/out, first the rule the windows are drawn from")
    ap.add_argument("--levels-from", default="resolved_cd99-or-hormone-top",
                    help="the member set that fixes the white points (same for every figure of an image)")
    ap.add_argument("--per-stratum", type=int, default=3)
    ap.add_argument("--rows-per-figure", type=int, default=4)
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()

    v, tag = args.detector, args.tag
    sets = [tuple(x.split("=", 1)) for x in args.sets]
    top = pd.read_parquet(args.top_marker)                  # islet_cell_test.py's table, or any table with a `top` column
    if "object_id" in top.columns:
        top = top.set_index("object_id")
    top = top["top"]
    black, white, n_used = display_levels(tag, args.levels_from, top)

    loaded = []
    for label, sub in sets:
        feats, geoms = objects(REV / v / "out" / sub / f"{tag}.geojson")
        ids = set(pd.read_csv(REV / v / "out" / sub / f"{tag}.members.tsv", sep="\t", usecols=["cell_id"])["cell_id"])
        loaded.append({"label": label, "feats": feats, "geoms": geoms, "tree": STRtree(geoms), "ids": ids})
    base = loaded[0]
    pf, pg = objects(REV / v / "out" / "pixel_objects" / f"{tag}.geojson")
    ptree = STRtree(pg)

    rng = np.random.default_rng(args.seed)
    windows = []
    for s_ in STRATA:
        pool = [i for i, f in enumerate(base["feats"]) if stratum(f["properties"]["measurements"]["n_cells"]) == s_]
        for i in sorted(rng.choice(pool, size=min(args.per_stratum, len(pool)), replace=False)):
            windows.append((s_, i))
    boxes = []
    for _, i in windows:
        x0, y0, x1, y1 = base["geoms"][i].bounds
        side = max(x1 - x0, y1 - y0) + 2 * PAD
        cx, cy = (x0 + x1) / 2, (y0 + y1) / 2
        boxes.append((int(max(0, cx - side / 2)), int(max(0, cy - side / 2)), int(cx + side / 2), int(cy + side / 2)))

    wtree = STRtree([box(*b) for b in boxes])
    cells = [[] for _ in boxes]
    with open(EXPORT / "cell_geojson" / f"{tag}.geojson", "rb") as h:
        for f in ijson.items(h, "features.item", use_float=True):
            g = shape(f["geometry"])
            for k in wtree.query(g, predicate="intersects"):
                cells[k].append((g, "nucleusGeometry" in f, f.get("id")))

    with tifffile.TiffFile(args.image) as tf:
        pages = tf.series[0].levels[0].pages
        names = []
        for p in pages:
            p = p.aspage() if hasattr(p, "aspage") else p
            m = re.search(r"<Biomarker>(.*?)</Biomarker>", p.description or "")
            names.append(m.group(1) if m else "?")
        planes = {c: pages[names.index(c)].asarray() for c in CHANNELS}

    def scaled(c: str, x0: int, y0: int, x1: int, y1: int) -> np.ndarray:
        a = planes[c][y0:y1, x0:x1].astype(np.float32)
        return np.clip((a - black[c]) / (white[c] - black[c]), 0, 1)

    levels = "; ".join(f"{c} black {black[c]:.0f} white {white[c]:.0f}" for c in CHANNELS)
    ncol = len(loaded) + 1
    figures = [list(range(i, min(i + args.rows_per_figure, len(boxes)))) for i in range(0, len(boxes), args.rows_per_figure)]
    for n_fig, rows in enumerate(figures, 1):
        fig, axes = plt.subplots(len(rows), ncol, figsize=(5.6 * ncol, 5.85 * len(rows)), squeeze=False)
        for r, k in enumerate(rows):
            (s_, i), (x0, y0, x1, y1) = windows[k], boxes[k]
            hormones = np.dstack([scaled(c, x0, y0, x1, y1) for c in HORMONES])
            cd99 = scaled("CD99", x0, y0, x1, y1)
            win = box(x0, y0, x1, y1)
            in_view = {cid for _, _, cid in cells[k]}
            for col in range(ncol):
                ax = axes[r][col]
                ax.axis("off")
                rule = loaded[min(col, len(loaded) - 1)]
                ax.imshow(hormones if col < len(loaded) else np.dstack([cd99] * 3), extent=(x0, x1, y1, y0), interpolation="nearest")
                for g, nuc, cid in cells[k]:
                    dropped = col > 0 and cid in base["ids"] and cid not in rule["ids"]
                    for part in parts(g):
                        if part.geom_type == "Polygon":
                            xs, ys = part.exterior.xy
                            ax.plot(xs, ys, color="orange" if dropped else ("white" if nuc else "grey"),
                                    lw=1.1 if dropped else 0.35, ls="-" if nuc else ":")
                    if dropped and col < len(loaded):
                        name = top.get(cid, "?")
                        ax.text(g.centroid.x, g.centroid.y, SHORT.get(name, name), color="orange", fontsize=7, ha="center", va="center")
                for j in ptree.query(win, predicate="intersects"):
                    for part in parts(pg[j]):
                        if part.geom_type == "Polygon":
                            xs, ys = part.exterior.xy
                            ax.plot(xs, ys, color="yellow", lw=0.9, ls="--")
                for j in rule["tree"].query(win, predicate="intersects"):
                    for part in parts(rule["geoms"][j]):
                        if part.geom_type == "Polygon":
                            xs, ys = part.exterior.xy
                            ax.plot(xs, ys, color="cyan", lw=1.7)
                n_members = len(in_view & rule["ids"])
                title = f"{rule['label']}: {n_members} member cells in view" if col < len(loaded) else f"CD99, object of: {rule['label']}"
                if col == 0:
                    title = f"window {k + 1} ({s_} cell(s) under the first rule, {(x1 - x0) * PX:.0f} um wide)\n" + title
                ax.set_title(title, fontsize=9)
                ax.set_xlim(x0, x1)
                ax.set_ylim(y1, y0)
        fig.suptitle(f"{tag.split('_Scan')[0]}, detector {v}.  Red INS, green GCG, blue SST; last panel CD99.  Linear; black = the image's background "
                     f"anchor, white = one level per image and channel ({levels}).\nYellow dashed = pixel object, cyan = cell-based object, "
                     "orange = a cell the first rule had and this rule drops (label: its highest marker)", fontsize=10)
        fig.tight_layout(rect=(0, 0, 1, 0.955))
        out = args.out.with_name(f"{args.out.stem}_{n_fig}{args.out.suffix}")
        fig.savefig(out, dpi=80)
        plt.close(fig)
        print("->", out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
