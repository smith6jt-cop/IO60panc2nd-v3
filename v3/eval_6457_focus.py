#!/usr/bin/env python3
"""eval_6457_focus.py -- where is DAPI in focus, per cycle? READ-ONLY on the share; writes only --out.

Every cycle of a PhenoCycler run images DAPI again, so a region that is hazy in the cycle the final
image took its DAPI from may be sharp in another. Per raw cycle file and per tile over the tissue:

  focus   var(Laplacian of the tile) / var(the tile): the share of the picture's contrast that sits
          in fine detail. A blurred or hazy tile has little; the measure does not move with exposure
          or background level (it is unchanged by gain and offset), unlike the Laplacian variance
          alone or over mean^2, which mostly follows how bright and how covered the tile is.
  contrast  (p95 - p5) / (p95 + p5) of the tile: how far the nuclei stand above their surroundings.
          A veil of background light lowers it without blurring anything, which `focus` cannot see.
  mean, std

Uses the tissue tiles and the file selection of eval_6457_cycles.py (run that first: its
cycle0_dapi_tiles.parquet and files.json are read from --eval).

    eval_6457_focus.py --eval <eval output dir> --out <dir> [--workers 4]
"""
from __future__ import annotations

import argparse
import json
import math
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.ndimage import laplace

import eval_6457_cycles as ev


def focus(tile: np.ndarray) -> tuple[float, float, float, float]:
    """(focus, contrast, mean, std) of one tile."""
    t = tile.astype(np.float32, copy=False)
    var = float(t.var())
    lo, hi = np.percentile(t[::2, ::2], [5, 95])
    contrast = float((hi - lo) / (hi + lo)) if hi + lo > 0 else float("nan")
    if var <= 0:
        return float("nan"), contrast, float(t.mean()), 0.0
    return float(laplace(t).var()) / var, contrast, float(t.mean()), math.sqrt(var)


def one_file(cycle: int, path: str, tiles: list[tuple[int, int]], tile: int) -> pd.DataFrame:
    plane = ev.read_plane(Path(path), 0)
    rows = []
    for i, j in tiles:
        t = plane[i * tile:(i + 1) * tile, j * tile:(j + 1) * tile]
        f, k, m, s = focus(t)
        rows.append((cycle, i, j, f, k, m, s))
    return pd.DataFrame(rows, columns=["cycle", "row", "col", "focus", "contrast", "mean", "std"])


def maps(d: pd.DataFrame, shape, tile: int, out: Path, label: str, covered: set):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    n_rows, n_cols = math.ceil(shape[0] / tile), math.ceil(shape[1] / tile)
    cycles = sorted(d["cycle"].unique())
    best = d.groupby(["row", "col"])["focus"].max().rename("best")
    d = d.join(best, on=["row", "col"])
    d["share_of_best"] = d["focus"] / d["best"]
    inside = d[[(r, c) in covered for r, c in zip(d["row"], d["col"])]]
    lo, hi = inside["focus"].quantile(0.01), inside["focus"].quantile(0.99)
    for col, vmin, vmax, cmap, title in (
            ("focus", lo, hi, "viridis", "DAPI focus per raw cycle (share of contrast in fine detail)"),
            ("contrast", 0.0, 1.0, "cividis", "DAPI contrast per raw cycle, (p95 - p5) / (p95 + p5)"),
            ("share_of_best", 0.3, 1.0, "magma", "DAPI focus as a share of the sharpest cycle in the same tile")):
        ncol = 8
        nrow = math.ceil(len(cycles) / ncol)
        fig, axes = plt.subplots(nrow, ncol, figsize=(ncol * 1.9, nrow * 2.9 + 0.6), squeeze=False)
        for ax in axes.ravel():
            ax.axis("off")
        im = None
        for ax, c in zip(axes.ravel(), cycles):
            s = d[d["cycle"] == c]
            g = np.full((n_rows, n_cols), np.nan)
            g[s["row"].to_numpy(), s["col"].to_numpy()] = s[col].to_numpy()
            im = ax.imshow(g, vmin=vmin, vmax=vmax, cmap=cmap, interpolation="nearest")
            ax.set_title(f"cycle {c}", fontsize=8)
        fig.suptitle(f"{label}: {title}\n(blank = no tissue; one tile = {tile * ev.PX_UM:.0f} um)", fontsize=9)
        fig.colorbar(im, ax=axes.ravel().tolist(), shrink=0.6)
        fig.savefig(out / f"map_dapi_{col}.png", dpi=130)
        plt.close(fig)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--eval", type=Path, required=True, help="output dir of eval_6457_cycles.py for this acquisition")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--workers", type=int, default=4)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    started = time.time()
    meta = json.loads((args.eval / "files.json").read_text())
    tile = int(meta["tile_px"])
    c0 = pd.read_parquet(args.eval / "cycle0_dapi_tiles.parquet")
    tiles = list(map(tuple, c0[["row", "col"]].to_numpy()))
    files = {int(c): k["raw"] for c, k in meta["cycles"].items() if k.get("raw")}
    from joblib import Parallel, delayed
    parts = Parallel(n_jobs=args.workers, verbose=5)(
        delayed(one_file)(c, f, tiles, tile) for c, f in sorted(files.items()))
    d = pd.concat(parts, ignore_index=True)
    d.to_parquet(args.out / "dapi_focus_tiles.parquet", index=False)

    # fully covered tiles, as eval_6457_report.interior: at least half the DAPI of a covered tile, no missing neighbour
    full = c0["dapi_mean"].quantile(0.75)
    have = set(tiles)
    covered = {(r, c) for r, c in map(tuple, c0.loc[c0["dapi_mean"] >= 0.5 * full, ["row", "col"]].to_numpy())
               if all((r + a, c + b) in have for a in (-1, 0, 1) for b in (-1, 0, 1))}
    inside = d[[(r, c) in covered for r, c in zip(d["row"], d["col"])]].copy()
    best = inside.groupby(["row", "col"])["focus"].max().rename("best")
    inside = inside.join(best, on=["row", "col"])
    inside["share_of_best"] = inside["focus"] / inside["best"]
    rows = []
    for c, s in inside.groupby("cycle"):
        rows.append({"cycle": int(c), "tiles": len(s), "contrast_median": s["contrast"].median(),
                     "contrast_p5": s["contrast"].quantile(0.05), "focus_median": s["focus"].median(),
                     "focus_p5": s["focus"].quantile(0.05), "focus_p95": s["focus"].quantile(0.95),
                     "share_of_best_median": s["share_of_best"].median(),
                     "tiles_below_0.75_of_best": float((s["share_of_best"] < 0.75).mean()),
                     "tiles_below_0.5_of_best": float((s["share_of_best"] < 0.5).mean()),
                     "tiles_where_this_cycle_is_sharpest": float((s["share_of_best"] >= 0.999).mean())})
    summary = pd.DataFrame(rows)
    summary.to_csv(args.out / "dapi_focus_by_cycle.tsv", sep="\t", index=False)
    shape = tuple(meta["shape"])
    maps(d, shape, tile, args.out, Path(meta["reference_raw"]).name.split("_Cycle")[0], covered)
    pd.set_option("display.width", 220)
    print(summary.round(3).to_string(index=False))
    print(f"-> {args.out}  ({time.time() - started:.0f} s)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
