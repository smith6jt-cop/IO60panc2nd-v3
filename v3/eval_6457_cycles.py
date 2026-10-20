#!/usr/bin/env python3
"""eval_6457_cycles.py -- Track 6457, stage 1: measure, from the instrument's per-cycle files, how well
the cycles are registered and where DAPI is sharp or hazy. READ-ONLY on the share; writes only --out.

The PhenoCycler Fusion leaves, per cycle, a stitched 4-channel file (DAPI + three reporters, on the
final image's grid) as `.raw` (offset and shading corrected) and `.intermediate` (registered by the
instrument and autofluorescence-subtracted). Every cycle carries its own DAPI, so three comparisons
are made per cycle, tile by tile over the tissue:

  raw_vs_raw0   the cycle's raw DAPI against cycle 0's raw DAPI: how far the slide really moved
  int_vs_raw    the cycle's intermediate DAPI against ITS OWN raw DAPI: what the instrument's
                registration applied (same picture twice, so this one is measured very reliably)
  int_vs_int1   the cycle's intermediate DAPI against the first intermediate's DAPI: the
                misregistration left between the processed cycles, which the final image inherits

per tile:
  ncc           normalised cross-correlation (the computation of KINTSUGI
                notebooks/Kprocess._compute_spatial_ncc_grid, restated here because one pass over the
                tiles also yields the shift and the sharpness)
  dx_um, dy_um  the shift that aligns the tile to its reference (peak of the windowed
                cross-correlation, 0.1 px steps)
  sharp         Laplacian variance of the tile / mean^2 (the metric of kintsugi.qc.boundary_check
                ._sharpness, made relative so a dim tile is not read as a blurred one)

and per cycle and comparison an affine fit to the tile shifts (translation, rotation, scale, shear)
with what is left after it: a shift field that is affine can be undone by re-registering, one that is
not cannot. At pyramid level 1 the same sharpness is measured for the three reporter channels, to see
whether the markers are blurred where DAPI is hazy. Nothing is judged here: tables and maps.

6457's folder holds TWO acquisitions (Sept 2025 cycles 0-17 under plain names; April 2026 cycles
0-22, whose cycles 0-17 carry " (2)"). --month pins the acquisition by each file's TIFF DateTime, so
a file from the other session is refused rather than silently mixed in.

    eval_6457_cycles.py --temp <.../Scan1/.temp> --prefix 6457_Scan1 --month 2026:04 --out <dir>
"""
from __future__ import annotations

import argparse
import json
import math
import re
import time
from pathlib import Path

import numpy as np
import pandas as pd
import tifffile
from scipy.ndimage import laplace
from skimage.registration import phase_cross_correlation

PX_UM = 0.5078
COMPARISONS = ("raw_vs_raw0", "int_vs_raw", "int_vs_int1")


def cycle_file(temp: Path, prefix: str, cycle: int, kind: str, month: str) -> Path | None:
    """The file of `kind` for `cycle` acquired in `month` ("YYYY:MM"), or None if there is none."""
    for name in (f"{prefix}_Cycle{cycle}.{kind} (2).qptiff", f"{prefix}_Cycle{cycle}.{kind}.qptiff"):
        f = temp / name
        if not f.exists():
            continue
        with tifffile.TiffFile(f) as tf:
            stamp = tf.pages[0].tags.get("DateTime")
            if stamp is not None and str(stamp.value)[:7].replace("-", ":") == month:
                return f
    return None


def read_plane(f: Path, channel: int, level: int = 0) -> np.ndarray:
    with tifffile.TiffFile(f) as tf:
        return tf.series[0].levels[level].pages[channel].asarray()


def channel_names(f: Path) -> list[str]:
    names = []
    with tifffile.TiffFile(f) as tf:
        for page in tf.series[0].levels[0].pages:
            page = page.aspage() if hasattr(page, "aspage") else page
            m = re.search(r"<Name>(.*?)</Name>", page.description or "")
            names.append(m.group(1) if m else "?")
    return names


def rel_sharpness(tile: np.ndarray) -> float:
    t = tile.astype(np.float32, copy=False)
    mean = float(t.mean())
    if mean <= 0:
        return float("nan")
    return float(laplace(t).var()) / (mean * mean)


def tiles_of(shape, tile):
    h, w = shape
    for i in range(math.ceil(h / tile)):
        for j in range(math.ceil(w / tile)):
            yield i, j, slice(i * tile, min((i + 1) * tile, h)), slice(j * tile, min((j + 1) * tile, w))


def tissue_grid(ref: np.ndarray, tile: int) -> np.ndarray:
    """True for tiles that hold tissue: the tile's 95th-percentile DAPI clears the valley between the
    empty and the tissue tiles (Otsu on the log of the tile values)."""
    from skimage.filters import threshold_otsu
    h, w = ref.shape
    grid = np.zeros((math.ceil(h / tile), math.ceil(w / tile)), dtype=np.float32)
    for i, j, ys, xs in tiles_of(ref.shape, tile):
        grid[i, j] = np.percentile(ref[ys, xs][::4, ::4], 95)
    logv = np.log1p(grid)
    return logv > threshold_otsu(logv)


def measure_pair(ref: np.ndarray, mov: np.ndarray, tissue: np.ndarray, tile: int) -> pd.DataFrame:
    """Per tissue tile: NCC and the shift (um) that aligns `mov` to `ref`; sharpness of `mov`."""
    if mov.shape != ref.shape:
        raise ValueError(f"shape {mov.shape} != reference {ref.shape}")
    window = np.outer(np.hanning(tile), np.hanning(tile)).astype(np.float32)
    rows = []
    for i, j, ys, xs in tiles_of(ref.shape, tile):
        if not tissue[i, j]:
            continue
        r = ref[ys, xs].astype(np.float32)
        m = mov[ys, xs].astype(np.float32)
        rc, mc = r - r.mean(), m - m.mean()
        rn, mn = float(np.linalg.norm(rc)), float(np.linalg.norm(mc))
        ncc = float(np.dot(rc.ravel(), mc.ravel()) / (rn * mn)) if rn > 1e-10 and mn > 1e-10 else float("nan")
        if r.shape == (tile, tile) and rn > 1e-10 and mn > 1e-10:
            # windowed, NOT phase-normalised: on blurred content the phase-normalised estimate without
            # a window collapses to zero shift (the tile edges, which never move, win), so a hazy tile
            # would read as perfectly registered -- measured on synthetic tiles, tests/test_eval_6457_cycles.py
            shift, error, _ = phase_cross_correlation(rc * window, mc * window, upsample_factor=10, normalization=None)
            dy, dx = float(shift[0]) * PX_UM, float(shift[1]) * PX_UM
        else:
            dy = dx = error = float("nan")
        rows.append((i, j, ncc, dx, dy, float(np.hypot(dx, dy)), float(error), rel_sharpness(m), float(m.mean())))
    return pd.DataFrame(rows, columns=["row", "col", "ncc", "dx_um", "dy_um", "shift_um", "phase_error",
                                       "dapi_sharp", "dapi_mean"])


def marker_sharpness(f: Path, tissue: np.ndarray, tile: int, cycle: int, kind: str) -> pd.DataFrame:
    """Relative sharpness of every channel of one cycle file at pyramid level 1 (same tile footprint)."""
    names = channel_names(f)
    half = tile // 2
    rows = []
    for c, name in enumerate(names):
        plane = read_plane(f, c, level=1)
        for i, j, ys, xs in tiles_of(plane.shape, half):
            if i < tissue.shape[0] and j < tissue.shape[1] and tissue[i, j]:
                t = plane[ys, xs]
                rows.append((cycle, kind, name, i, j, rel_sharpness(t), float(t.mean())))
    return pd.DataFrame(rows, columns=["cycle", "kind", "channel", "row", "col", "sharp_l1", "mean_l1"])


def one_cycle(temp, prefix, month, cycle, raw0_file, int1_file, int1_cycle, tissue, tile, raw_only=False):
    """The three comparisons for one cycle, holding at most two full-resolution planes at a time.
    With raw_only, just the raw DAPI against the raw reference (which may be another cycle than 0)."""
    reg, mark = [], []
    if raw_only:
        raw_f = cycle_file(temp, prefix, cycle, "raw", month)
        if raw_f is None:
            return None, None
        d = measure_pair(read_plane(raw0_file, 0), read_plane(raw_f, 0), tissue, tile)
        d.insert(0, "comparison", "raw_vs_raw0")
        d.insert(0, "cycle", cycle)
        return d, None
    raw_f = cycle_file(temp, prefix, cycle, "raw", month)
    int_f = cycle_file(temp, prefix, cycle, "intermediate", month)

    def add(name, ref, mov):
        d = measure_pair(ref, mov, tissue, tile)
        d.insert(0, "comparison", name)
        d.insert(0, "cycle", cycle)
        reg.append(d)

    raw = read_plane(raw_f, 0) if raw_f else None
    if raw is not None:
        add("raw_vs_raw0", read_plane(raw0_file, 0), raw)
        mark.append(marker_sharpness(raw_f, tissue, tile, cycle, "raw"))
    inter = read_plane(int_f, 0) if int_f else None
    if inter is not None and raw is not None:
        add("int_vs_raw", raw, inter)
    del raw
    if inter is not None:
        if int1_file is not None and cycle != int1_cycle:
            add("int_vs_int1", read_plane(int1_file, 0), inter)
        mark.append(marker_sharpness(int_f, tissue, tile, cycle, "intermediate"))
    return (pd.concat(reg) if reg else None), (pd.concat(mark) if mark else None)


def affine_fit(d: pd.DataFrame, tile: int) -> dict:
    """Least-squares affine model of the tile shifts, refitted once without tiles further than three
    robust deviations from it. Lengths in um; the linear part is dimensionless, reported in ppm
    (1 ppm = 1 um per metre; 100 ppm over a 25 mm section = 2.5 um) and the rotation in millidegrees."""
    ok = d["dx_um"].notna() & d["dy_um"].notna()
    d = d[ok]
    out = {"tiles_fit": int(len(d))}
    if len(d) < 30:
        return out
    x = (d["col"].to_numpy() + 0.5) * tile * PX_UM
    y = (d["row"].to_numpy() + 0.5) * tile * PX_UM
    x0, y0 = x.mean(), y.mean()
    a = np.column_stack([np.ones_like(x), x - x0, y - y0])
    u = np.column_stack([d["dx_um"].to_numpy(), d["dy_um"].to_numpy()])
    keep = np.ones(len(d), dtype=bool)
    for _ in range(2):
        coef, *_ = np.linalg.lstsq(a[keep], u[keep], rcond=None)
        res = np.hypot(*(u - a @ coef).T)
        mad = np.median(np.abs(res[keep] - np.median(res[keep]))) * 1.4826
        keep = res <= np.median(res[keep]) + 3 * max(mad, 0.02)
    (tx, ty), (a1, b1), (a2, b2) = coef            # dx = tx + a1 (x-x0) + a2 (y-y0); dy = ty + b1 (x-x0) + b2 (y-y0)
    model = np.hypot(*(a @ coef).T)
    out.update({
        "tiles_fit": int(keep.sum()), "tiles_outlier_share": float(1 - keep.mean()),
        "translation_um": float(np.hypot(tx, ty)), "tx_um": float(tx), "ty_um": float(ty),
        "rotation_mdeg": float(np.degrees((b1 - a2) / 2) * 1000),
        "scale_ppm": float((a1 + b2) / 2 * 1e6), "stretch_ppm": float(np.hypot((a1 - b2) / 2, (a2 + b1) / 2) * 1e6),
        "affine_shift_um_p50": float(np.median(model)), "affine_shift_um_max": float(model.max()),
        "residual_um_p50": float(np.median(res)), "residual_um_p95": float(np.percentile(res, 95)),
    })
    return out


def maps(reg: pd.DataFrame, ref0: pd.DataFrame, shape, tile, out: Path, label: str):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    n_rows, n_cols = math.ceil(shape[0] / tile), math.ceil(shape[1] / tile)
    base = np.full((n_rows, n_cols), np.nan)
    base[ref0["row"].to_numpy(), ref0["col"].to_numpy()] = ref0["dapi_sharp"].to_numpy()
    specs = [
        ("shift_um", "raw_vs_raw0", 0, 3, "magma", "how far the slide moved: raw DAPI vs cycle 0 raw DAPI, um"),
        ("shift_um", "int_vs_raw", 0, 6, "magma", "what the instrument's registration applied: intermediate vs its own raw DAPI, um"),
        ("shift_um", "int_vs_int1", 0, 6, "magma", "misregistration left between processed cycles: intermediate vs first intermediate DAPI, um"),
        ("dx_um", "int_vs_int1", -6, 6, "coolwarm", "the same, x component, um"),
        ("dy_um", "int_vs_int1", -6, 6, "coolwarm", "the same, y component, um"),
        ("ncc", "int_vs_int1", 0, 1, "viridis", "DAPI correlation: intermediate vs first intermediate"),
        ("ncc", "raw_vs_raw0", 0, 1, "viridis", "DAPI correlation: raw vs cycle 0 raw"),
        ("sharp_log10", "raw_vs_raw0", None, None, "cividis", "raw DAPI relative sharpness (log10)"),
        ("sharp_ratio_log2", "raw_vs_raw0", -2, 2, "PuOr", "raw DAPI sharpness relative to cycle 0 in the same tile (log2)"),
    ]
    for col, comp, vmin, vmax, cmap, title in specs:
        d = reg[reg["comparison"] == comp]
        cycles = sorted(d["cycle"].unique())
        if not cycles:
            continue
        ncol = 8
        nrow = math.ceil(len(cycles) / ncol)
        fig, axes = plt.subplots(nrow, ncol, figsize=(ncol * 1.9, nrow * 2.9 + 0.6), squeeze=False)
        grids = {}
        for c in cycles:
            s = d[d["cycle"] == c]
            g = np.full((n_rows, n_cols), np.nan)
            if col == "sharp_log10":
                v = np.log10(s["dapi_sharp"].clip(lower=1e-9))
            elif col == "sharp_ratio_log2":
                v = np.log2(s["dapi_sharp"].clip(lower=1e-9) / np.clip(base[s["row"].to_numpy(), s["col"].to_numpy()], 1e-9, None))
            else:
                v = s[col]
            g[s["row"].to_numpy(), s["col"].to_numpy()] = np.asarray(v)
            grids[c] = g
        allv = np.concatenate([g[np.isfinite(g)] for g in grids.values()])
        lo = float(np.percentile(allv, 1)) if vmin is None else vmin
        hi = float(np.percentile(allv, 99)) if vmax is None else vmax
        im = None
        for ax in axes.ravel():
            ax.axis("off")
        for ax, c in zip(axes.ravel(), cycles):
            im = ax.imshow(grids[c], vmin=lo, vmax=hi, cmap=cmap, interpolation="nearest")
            ax.set_title(f"cycle {c}", fontsize=8)
        fig.suptitle(f"{label}: {title}\n(blank = no tissue; one tile = {tile * PX_UM:.0f} um)", fontsize=9)
        fig.colorbar(im, ax=axes.ravel().tolist(), shrink=0.6)
        fig.savefig(out / f"map_{comp}_{col}.png", dpi=130)
        plt.close(fig)


def final_dapi_source(final: Path, files: dict, out: Path, level: int = 3):
    """Which cycle's DAPI the final image carries: correlation of the final's DAPI with every cycle's
    intermediate (and raw) DAPI at a low pyramid level."""
    with tifffile.TiffFile(final) as tf:
        fin = tf.series[0].levels[level].pages[0].asarray().astype(np.float32)
    rows = []
    for c, kinds in sorted(files.items()):
        for kind, f in kinds.items():
            if not f:
                continue
            p = read_plane(Path(f), 0, level=level).astype(np.float32)
            if p.shape != fin.shape:
                rows.append((c, kind, float("nan"), f"shape {p.shape} != final {fin.shape}"))
                continue
            a, b = fin - fin.mean(), p - p.mean()
            rows.append((c, kind, float((a * b).sum() / (np.linalg.norm(a) * np.linalg.norm(b))), ""))
    pd.DataFrame(rows, columns=["cycle", "kind", "ncc_with_final_dapi", "note"]).to_csv(
        out / "final_dapi_source.tsv", sep="\t", index=False)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--temp", type=Path, required=True)
    ap.add_argument("--prefix", required=True, help="e.g. 6457_Scan1")
    ap.add_argument("--month", required=True, help='acquisition month as the TIFF DateTime starts, e.g. "2026:04"')
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--tile", type=int, default=512, help="tile edge in full-resolution pixels")
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--cycles", default="", help="comma list; default every cycle found")
    ap.add_argument("--final", type=Path, default=None, help="the final combined qptiff (optional)")
    ap.add_argument("--raw-reference-cycle", type=int, default=0,
                    help="compare every raw DAPI against THIS cycle's raw DAPI instead of cycle 0's, and nothing "
                         "else (no intermediates, no reporter channels): how well the raw cycles agree with one "
                         "another when cycle 0 is the odd one out. The tissue tiles are still found on cycle 0.")
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    started = time.time()

    raw0_file = cycle_file(args.temp, args.prefix, 0, "raw", args.month)
    if raw0_file is None:
        raise SystemExit(f"no cycle-0 raw file acquired in {args.month} under {args.temp}")
    wanted = [int(c) for c in args.cycles.split(",") if c] or list(range(1, 40))
    found = {}
    for c in range(0, 40):
        kinds = {k: cycle_file(args.temp, args.prefix, c, k, args.month) for k in ("raw", "intermediate")}
        if any(kinds.values()):
            found[c] = {k: (str(v) if v else None) for k, v in kinds.items()}
    int_cycles = [c for c in sorted(found) if found[c]["intermediate"]]
    int1_cycle = int_cycles[0] if int_cycles else None
    int1_file = Path(found[int1_cycle]["intermediate"]) if int1_cycle is not None else None
    raw_only = args.raw_reference_cycle != 0
    run = [c for c in sorted(found) if c in wanted and c != 0 and c != args.raw_reference_cycle]

    ref = read_plane(raw0_file, 0)
    tissue = tissue_grid(ref, args.tile)
    ref_sharp = [(i, j, rel_sharpness(ref[ys, xs]), float(ref[ys, xs].mean()))
                 for i, j, ys, xs in tiles_of(ref.shape, args.tile) if tissue[i, j]]
    shape = ref.shape
    del ref
    print(f"reference {raw0_file.name}: {shape}, {int(tissue.sum())} tissue tiles of {tissue.size}; "
          f"first intermediate: cycle {int1_cycle}; cycles {run}", flush=True)
    raw_ref_file = raw0_file
    if raw_only:
        raw_ref_file = cycle_file(args.temp, args.prefix, args.raw_reference_cycle, "raw", args.month)
        if raw_ref_file is None:
            raise SystemExit(f"no raw file for reference cycle {args.raw_reference_cycle}")
        print(f"raw reference: {raw_ref_file.name} (raw DAPI only)", flush=True)

    from joblib import Parallel, delayed
    results = Parallel(n_jobs=args.workers, verbose=5)(
        delayed(one_cycle)(args.temp, args.prefix, args.month, c, raw_ref_file, int1_file, int1_cycle, tissue,
                           args.tile, raw_only)
        for c in run)
    reg = pd.concat([r for r, _ in results if r is not None], ignore_index=True)
    marks = [m for _, m in results if m is not None]
    mark = pd.concat(marks, ignore_index=True) if marks else pd.DataFrame(
        columns=["cycle", "kind", "channel", "row", "col", "sharp_l1", "mean_l1"])
    ref0 = pd.DataFrame(ref_sharp, columns=["row", "col", "dapi_sharp", "dapi_mean"])
    reg.to_parquet(args.out / "registration_tiles.parquet", index=False)
    mark.to_parquet(args.out / "channel_sharpness_tiles.parquet", index=False)
    ref0.to_parquet(args.out / "cycle0_dapi_tiles.parquet", index=False)

    # ---- per cycle and comparison: plain numbers, no verdict ----
    def q(s, p):
        return float(np.nanpercentile(s, p)) if s.notna().any() else float("nan")
    rows = []
    for (cycle, comp), d in reg.groupby(["cycle", "comparison"]):
        row = {"cycle": cycle, "comparison": comp, "tissue_tiles": len(d),
               "ncc_median": float(d["ncc"].median()), "ncc_p5": q(d["ncc"], 5),
               "shift_um_median": float(d["shift_um"].median()), "shift_um_p95": q(d["shift_um"], 95),
               "tiles_shift_over_1um": float((d["shift_um"] > 1.0).mean()),
               "tiles_shift_over_3um": float((d["shift_um"] > 3.0).mean())}
        row.update(affine_fit(d, args.tile))
        rows.append(row)
    summary = pd.DataFrame(rows).sort_values(["comparison", "cycle"])
    summary.to_csv(args.out / "cycle_summary.tsv", sep="\t", index=False)

    # ---- sharpness: DAPI per cycle against cycle 0, and the reporters where DAPI is least / most sharp ----
    raw = reg[reg["comparison"] == "raw_vs_raw0"].merge(ref0, on=["row", "col"], suffixes=("", "_c0"))
    rows = []
    for cycle, d in raw.groupby("cycle"):
        ratio = d["dapi_sharp"] / d["dapi_sharp_c0"]
        rows.append({"cycle": cycle, "dapi_sharp_median": float(d["dapi_sharp"].median()),
                     "cycle0_dapi_sharp_median": float(d["dapi_sharp_c0"].median()),
                     "ratio_to_cycle0_median": float(ratio.median()),
                     "tiles_below_half_of_cycle0": float((ratio < 0.5).mean())})
    pd.DataFrame(rows).to_csv(args.out / "dapi_sharpness_by_cycle.tsv", sep="\t", index=False)
    m = mark[(mark["kind"] == "raw")].merge(raw[["cycle", "row", "col", "dapi_sharp"]], on=["cycle", "row", "col"])
    rows = []
    for (cycle, channel), d in m.groupby(["cycle", "channel"]):
        lo, hi = d["dapi_sharp"].quantile([0.25, 0.75])
        a, b = d.loc[d["dapi_sharp"] <= lo, "sharp_l1"], d.loc[d["dapi_sharp"] >= hi, "sharp_l1"]
        rows.append({"cycle": cycle, "channel": channel, "tiles": len(d),
                     "sharp_where_dapi_least_sharp_quarter": float(a.median()),
                     "sharp_where_dapi_sharpest_quarter": float(b.median()),
                     "ratio": float(a.median() / b.median()) if b.median() > 0 else float("nan"),
                     "spearman_with_dapi_sharpness": float(d["sharp_l1"].corr(d["dapi_sharp"], method="spearman"))})
    pd.DataFrame(rows).to_csv(args.out / "reporter_sharpness_vs_dapi.tsv", sep="\t", index=False)

    (args.out / "files.json").write_text(json.dumps(
        {"reference_raw": str(raw0_file), "raw_reference_cycle": args.raw_reference_cycle,
         "first_intermediate_cycle": int1_cycle, "month": args.month, "cycles": found,
         "shape": shape, "tile_px": args.tile, "pixel_um": PX_UM, "seconds": round(time.time() - started, 1)}, indent=2))
    maps(reg, ref0, shape, args.tile, args.out, args.prefix)
    if args.final is not None and args.final.exists():
        final_dapi_source(args.final, found, args.out)
    pd.set_option("display.width", 250)
    pd.set_option("display.max_columns", 30)
    print(summary.round(3).to_string(index=False))
    print(f"-> {args.out}  ({time.time() - started:.0f} s)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
