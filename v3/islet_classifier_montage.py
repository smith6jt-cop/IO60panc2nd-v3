#!/usr/bin/env python3
"""islet_classifier_montage.py -- the islet objects on the image, for looking without QuPath.

Each window is shown twice: the three hormones (INS red, GCG green, SST blue) and CD99 (grey), both windowed on the
image's own levels (black = the image's background anchor, white = 1.5 x the median of its confirmed endocrine
cells), so a dim image is not made to look bright.

  the objects as they end, on cell segmentation outlines (step 1b; maintainer, 2026-10-02), in the project's own
  class colours:
  red line       Islet: two or more whole cells, 35 um across or more
  orange line    EndocrineCluster: two or more whole cells, under 35 um across
  cyan line      EndocrineSingle: one whole cell; the outline is that cell's own
  yellow line    Hormone Fragment: no whole cell fits; the outline is the hormone's own
  white dashes   the outline script 1 drew with the image's own classifier, before 1b (--script1), or else the
                 objects script 1 makes with Islets2 (--baseline)
  thin grey      cell segmentation outlines (--cell-outlines; dotted = a detection without a nucleus)
  white dot      a cell inside an Islet, EndocrineCluster or EndocrineSingle
  magenta +      a confirmed endocrine cell (a hormone is its highest marker AND its CD99 is at the image's
                 reference level or above) that those objects leave outside

Windows: --per-stratum Islets for each size (1 / 2-3 / 4-10 / >10 cells), chosen at random with a fixed seed, then
--fragments Hormone Fragments, then --misses windows on confirmed cells left outside -- the failures are shown, not
hidden.

    islet_classifier_montage.py --tag TAG --image IMG.qptiff --objects trial/out/reviewed_islets/TAG.geojson \\
        [--script1 trial/out/script1_outlines/TAG.geojson] [--cell-outlines <export>/cell_geojson/TAG.geojson] \\
        [--baseline islets2/TAG.geojson] --cells cells/TAG.parquet --redsea redsea/donor_id=D/data_0.parquet \\
        --structures <export>/structures --anchors <export>/islet_anchors.csv --out TAG.png
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
import shapely
from shapely import STRtree

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
import islet_cell_test as ict  # noqa: E402
import islet_classifier_check as check  # noqa: E402
import islet_classifier_fit as fit  # noqa: E402
import islet_outlines_to_cells as oc  # noqa: E402

PX0 = fit.CLASSIFIER_UM / 2.0


def outline(ax, geoms, box, **style):
    x0, y0, x1, y1 = box
    for g in geoms:
        b = g.bounds
        if b[2] < x0 or b[0] > x1 or b[3] < y0 or b[1] > y1:
            continue
        for part in (g.geoms if hasattr(g, "geoms") else [g]):
            if part.geom_type != "Polygon":
                continue
            xs, ys = part.exterior.xy
            ax.plot(xs, ys, **style)
            for ring in part.interiors:
                xs, ys = ring.xy
                ax.plot(xs, ys, **style)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tag", required=True)
    ap.add_argument("--image", type=Path, required=True)
    ap.add_argument("--objects", type=Path, required=True)
    ap.add_argument("--script1", type=Path, default=None, help="script 1's outlines, exported before step 1b")
    ap.add_argument("--cell-outlines", type=Path, default=None, help="<export>/cell_geojson/TAG.geojson")
    ap.add_argument("--fragments", type=int, default=6)
    ap.add_argument("--baseline", type=Path, default=None)
    ap.add_argument("--cells", type=Path, required=True)
    ap.add_argument("--redsea", type=Path, required=True)
    ap.add_argument("--structures", type=Path, default=None)
    ap.add_argument("--anchors", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--per-stratum", type=int, default=3)
    ap.add_argument("--misses", type=int, default=6)
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--label", default=None)
    args = ap.parse_args()

    base, values, x_um, y_um, eligible, in_structure, other_markers, quantile, _ = fit.load_cells(
        args.cells, args.redsea, args.structures, args.tag, ict.SENIOR_TYPING)
    lab = fit.label_cells(values, x_um, y_um, eligible, in_structure, other_markers, quantile)
    x0c, y0c = x_um / PX0, y_um / PX0
    read = check.read_by_class(args.objects)
    final = read["by_class"]
    new = [e["geometry"] for e in read["endocrine"]]
    new_class = np.array([e["class"] for e in read["endocrine"]], dtype=object)
    new_across = np.array([oc.mean_caliper(e["geometry"]) * PX0 for e in read["endocrine"]])
    fragments = final.get(oc.CLASS_FRAGMENT, [])
    old = (check.read_by_class(args.script1)["by_class"].get(oc.CLASS_ISLET, []) if args.script1
           else (ict.read_objects(args.baseline) if args.baseline else []))
    hit, which = STRtree(new).query(shapely.points(x0c, y0c), predicate="within")
    inside = np.zeros(len(base), dtype=bool); inside[hit] = True; inside &= eligible
    per = np.bincount(which[eligible[hit]], minlength=len(new))

    import tifffile  # noqa: PLC0415
    with tifffile.TiffFile(args.image) as tf:
        pages = tf.series[0].levels[0].pages
        names = []
        for p in pages:
            p = p.aspage() if hasattr(p, "aspage") else p
            m = re.search(r"<Biomarker>(.*?)</Biomarker>", p.description or "")
            names.append(m.group(1) if m else "?")
        planes = {c: pages[names.index(c)].asarray() for c in fit.CHANNELS}
    anchors = pd.read_csv(args.anchors)
    black = anchors[anchors["tag"] == args.tag].set_index("channel")["anchor"].to_dict()
    raw = pd.read_parquet(args.cells, columns=["Object ID"] + [f"Cell: {c}: Mean" for c in fit.CHANNELS]).set_index("Object ID").reindex(base.index)
    white = {}
    for c in fit.CHANNELS:
        on = lab["confirmed"] & ((lab["top"] == c) if c != "CD99" else True)
        white[c] = 1.5 * float(np.nanmedian(raw[f"Cell: {c}: Mean"].to_numpy(dtype=float)[on])) if on.any() else float(np.percentile(planes[c][::16, ::16], 99.9))

    rng = np.random.default_rng(args.seed)
    wins = []
    for name, lo, hi in (("1 cell", 1, 1), ("2-3 cells", 2, 3), ("4-10 cells", 4, 10), ("more than 10 cells", 11, 10 ** 9)):
        pool = np.flatnonzero((per >= lo) & (per <= hi))
        for i in sorted(rng.choice(pool, size=min(args.per_stratum, len(pool)), replace=False)):
            wins.append((f"{name} ({per[i]}): {new_class[i]}, {new_across[i]:.0f} um", new[i].bounds))
    for i in sorted(rng.choice(len(fragments), size=min(args.fragments, len(fragments)), replace=False)):
        wins.append(("Hormone Fragment", fragments[i].bounds))
    missed = np.flatnonzero(lab["confirmed"] & ~inside)
    for i in sorted(rng.choice(missed, size=min(args.misses, len(missed)), replace=False)):
        wins.append(("a confirmed endocrine cell left outside", (x0c[i] - 4, y0c[i] - 4, x0c[i] + 4, y0c[i] + 4)))
    cell_polys, cell_nucleus = [], np.zeros(0, dtype=bool)
    if args.cell_outlines is not None:
        boxes = [shapely.box(b[0] - 80, b[1] - 80, b[2] + 80, b[3] + 80) for _, b in wins]
        _, cell_polys, cell_nucleus = oc.cells_near(args.cell_outlines, boxes)

    cols = 3
    rows = int(np.ceil(len(wins) / cols))
    fig, axes = plt.subplots(rows, cols * 2, figsize=(cols * 10, rows * 5.1), squeeze=False)
    for ax in axes.ravel():
        ax.axis("off")
    hgt, wid = planes["CD99"].shape
    for k, (title, (bx0, by0, bx1, by1)) in enumerate(wins):
        side = max(bx1 - bx0, by1 - by0, 60) + 90
        cxw, cyw = (bx0 + bx1) / 2, (by0 + by1) / 2
        x0, y0 = int(max(0, cxw - side / 2)), int(max(0, cyw - side / 2))
        x1, y1 = int(min(wid, cxw + side / 2)), int(min(hgt, cyw + side / 2))
        scaled = lambda c: np.clip((planes[c][y0:y1, x0:x1].astype(np.float32) - black[c]) / max(white[c] - black[c], 1.0), 0, 1)  # noqa: E731
        for j, picture in enumerate((np.dstack([scaled("INS"), scaled("GCG"), scaled("SST")]), np.dstack([scaled("CD99")] * 3))):
            ax = axes[k // cols, (k % cols) * 2 + j]
            ax.imshow(picture, extent=(x0, x1, y1, y0), interpolation="nearest")
            if len(cell_polys):
                outline(ax, cell_polys[cell_nucleus], (x0, y0, x1, y1), color="#b0b0b0", lw=0.4)
                outline(ax, cell_polys[~cell_nucleus], (x0, y0, x1, y1), color="#b0b0b0", lw=0.4, ls=":")
            outline(ax, old, (x0, y0, x1, y1), color="white", lw=0.9, ls="--")
            outline(ax, fragments, (x0, y0, x1, y1), color="yellow", lw=1.5)
            for cls, colour in ((oc.CLASS_SINGLE, "cyan"), (oc.CLASS_CLUSTER, "orange"), (oc.CLASS_ISLET, "red")):
                outline(ax, [g for g, c in zip(new, new_class) if c == cls], (x0, y0, x1, y1), color=colour, lw=1.5)
            here = (x0c >= x0) & (x0c < x1) & (y0c >= y0) & (y0c < y1) & eligible
            ax.scatter(x0c[here & ~inside], y0c[here & ~inside], s=4, c="grey", marker=".")
            ax.scatter(x0c[here & inside], y0c[here & inside], s=14, c="white", marker=".")
            out_conf = here & ~inside & lab["confirmed"]
            ax.scatter(x0c[out_conf], y0c[out_conf], s=60, c="magenta", marker="+")
            ax.set_xlim(x0, x1); ax.set_ylim(y1, y0)
            ax.set_title(f"{k + 1}. {title} -- {(x1 - x0) * PX0:.0f} um across" if j == 0 else "CD99", fontsize=10)
    before = "the outline script 1 drew" if args.script1 else "Islets2"
    fig.suptitle(f"{args.label or args.tag}: red = Islet (35 um or more), orange = EndocrineCluster, cyan = EndocrineSingle (one cell), yellow = Hormone Fragment; "
                 f"white dashes = {before}; grey = cell outlines. Left INS red / GCG green / SST blue, right CD99. "
                 "White dot = cell inside one of them; magenta + = confirmed endocrine cell left outside.", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.975))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=70)
    print("->", args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
