#!/usr/bin/env python3
"""eval_6457_report.py -- the numbers of eval_6457_cycles.py, per cycle, on INTERIOR tissue tiles.

A tile on the edge of the tissue is mostly glass; its cross-correlation peak is noise (shifts of
30-40 um in every cycle, raw or not), so the per-cycle numbers here use only tiles that are covered
by tissue (see `interior`). Reads the tile tables; touches no image. No verdict is written.

    eval_6457_report.py <eval output dir> [<control output dir>] --out report.md
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

import eval_6457_cycles as ev


def interior(reg: pd.DataFrame, c0: pd.DataFrame) -> pd.DataFrame:
    """Tiles that are tissue through and through: cycle 0's DAPI mean is at least half that of a
    fully covered tile (the 75th percentile of the tissue tiles), and no neighbour is off the tissue.

    Measured on 6457: the tiles that give a 30 um "shift" in every raw cycle are the partly covered
    ones (cycle-0 DAPI mean 15 against 75) -- glass, where the correlation peak is noise.
    """
    full = c0["dapi_mean"].quantile(0.75)
    covered = set(map(tuple, c0.loc[c0["dapi_mean"] >= 0.5 * full, ["row", "col"]].to_numpy()))
    tiles = set(map(tuple, c0[["row", "col"]].to_numpy()))
    keep = {(r, c) for r, c in covered
            if all((r + dr, c + dc) in tiles for dr in (-1, 0, 1) for dc in (-1, 0, 1))}
    mask = [(r, c) in keep for r, c in zip(reg["row"], reg["col"])]
    return reg[mask]


def per_cycle(reg: pd.DataFrame, tile: int) -> pd.DataFrame:
    rows = []
    for (comp, cycle), d in reg.groupby(["comparison", "cycle"]):
        fit = ev.affine_fit(d, tile)
        rows.append({
            "comparison": comp, "cycle": int(cycle), "tiles": len(d),
            "shift_um_median": d["shift_um"].median(), "shift_um_p95": d["shift_um"].quantile(0.95),
            "shift_um_max": d["shift_um"].max(),
            "tiles_over_1um": (d["shift_um"] > 1).mean(), "tiles_over_3um": (d["shift_um"] > 3).mean(),
            "ncc_median": d["ncc"].median(),
            "affine_max_um": fit.get("affine_shift_um_max"), "scale_ppm": fit.get("scale_ppm"),
            "stretch_ppm": fit.get("stretch_ppm"), "rotation_mdeg": fit.get("rotation_mdeg"),
            "translation_um": fit.get("translation_um"),
            "left_after_affine_um_median": fit.get("residual_um_p50"),
            "left_after_affine_um_p95": fit.get("residual_um_p95"),
        })
    return pd.DataFrame(rows).sort_values(["comparison", "cycle"])


def sharpness(d: Path, reg: pd.DataFrame) -> pd.DataFrame:
    c0 = pd.read_parquet(d / "cycle0_dapi_tiles.parquet")
    raw = reg[reg["comparison"] == "raw_vs_raw0"].merge(c0, on=["row", "col"], suffixes=("", "_c0"))
    raw["ratio"] = raw["dapi_sharp"] / raw["dapi_sharp_c0"]
    out = raw.groupby("cycle").agg(dapi_sharp_median=("dapi_sharp", "median"),
                                   ratio_to_cycle0_median=("ratio", "median"),
                                   tiles_below_half_of_cycle0=("ratio", lambda s: float((s < 0.5).mean()))).reset_index()
    out.attrs["cycle0"] = c0
    return out


def table(df: pd.DataFrame, cols: list[str], digits: int = 2) -> str:
    head = "| " + " | ".join(cols) + " |\n|" + "---|" * len(cols) + "\n"
    body = ""
    for _, r in df.iterrows():
        cells = []
        for c in cols:
            v = r[c]
            cells.append(str(int(v)) if c in ("cycle", "tiles") else ("" if pd.isna(v) else f"{v:.{digits}f}"))
        body += "| " + " | ".join(cells) + " |\n"
    return head + body


def describe(d: Path, label: str) -> tuple[str, dict]:
    meta = json.loads((d / "files.json").read_text())
    tile = int(meta["tile_px"])
    reg_all = pd.read_parquet(d / "registration_tiles.parquet")
    reg = interior(reg_all, pd.read_parquet(d / "cycle0_dapi_tiles.parquet"))
    n_all = reg_all[["row", "col"]].drop_duplicates().shape[0]
    n_int = reg[["row", "col"]].drop_duplicates().shape[0]
    pc = per_cycle(reg, tile)
    pc.to_csv(d / "cycle_summary_interior.tsv", sep="\t", index=False)
    sh = sharpness(d, reg)
    c0 = sh.attrs["cycle0"]
    c0_int = c0.merge(reg[["row", "col"]].drop_duplicates(), on=["row", "col"])
    text = f"## {label}\n\n"
    text += (f"Files: `{meta['reference_raw']}` (cycle 0 raw) and {len(meta['cycles']) - 1} later cycles, acquired "
             f"{meta['month'].replace(':', '-')}. {n_all:,} tissue tiles of {tile * ev.PX_UM:.0f} um; {n_int:,} fully covered tiles used below.\n\n")
    names = {"raw_vs_raw0": "Raw DAPI of each cycle against cycle 0's raw DAPI (how far the slide really moved)",
             "int_vs_raw": "Intermediate DAPI against the same cycle's raw DAPI (what the instrument's registration applied)",
             "int_vs_int1": f"Intermediate DAPI against the first intermediate, cycle {meta['first_intermediate_cycle']} (left between processed cycles)"}
    cols = ["cycle", "shift_um_median", "shift_um_p95", "tiles_over_1um", "tiles_over_3um", "affine_max_um",
            "scale_ppm", "stretch_ppm", "rotation_mdeg", "left_after_affine_um_median", "left_after_affine_um_p95", "ncc_median"]
    stats = {}
    for comp, title in names.items():
        sub = pc[pc["comparison"] == comp]
        if sub.empty:
            continue
        text += f"### {title}\n\n" + table(sub, cols) + "\n"
        stats[comp] = {"median_of_cycle_medians_um": float(sub["shift_um_median"].median()),
                       "worst_cycle_median_um": float(sub["shift_um_median"].max()),
                       "worst_cycle": int(sub.loc[sub["shift_um_median"].idxmax(), "cycle"]),
                       "cycles_median_over_1um": int((sub["shift_um_median"] > 1).sum()), "cycles": int(len(sub)),
                       "worst_left_after_affine_median_um": float(sub["left_after_affine_um_median"].max())}
    text += "### DAPI sharpness (Laplacian variance / mean^2; raw files)\n\n"
    text += (f"Cycle 0: median {c0_int['dapi_sharp'].median():.3f}, 5th-95th percentile of tiles "
             f"{c0_int['dapi_sharp'].quantile(0.05):.3f}-{c0_int['dapi_sharp'].quantile(0.95):.3f}.\n\n")
    text += table(sh, ["cycle", "dapi_sharp_median", "ratio_to_cycle0_median", "tiles_below_half_of_cycle0"], 3) + "\n"
    stats["cycle0_sharp_median"] = float(c0_int["dapi_sharp"].median())
    stats["cycle0_sharp_p5"] = float(c0_int["dapi_sharp"].quantile(0.05))
    stats["later_cycles_ratio_median"] = float(sh["ratio_to_cycle0_median"].median())
    rep = d / "reporter_sharpness_vs_dapi.tsv"
    if rep.exists():
        r = pd.read_csv(rep, sep="\t")
        r = r[r["channel"] != "DAPI"]
        text += ("### Reporter channels where DAPI is least sharp vs sharpest (same cycle, raw files)\n\n"
                 f"Ratio of the reporter's median tile sharpness in the least-sharp DAPI quarter to the sharpest quarter, "
                 f"over {len(r)} cycle x channel pairs: median {r['ratio'].median():.2f}, "
                 f"5th-95th percentile {r['ratio'].quantile(0.05):.2f}-{r['ratio'].quantile(0.95):.2f}; "
                 f"Spearman correlation with DAPI sharpness: median {r['spearman_with_dapi_sharpness'].median():.2f}.\n\n")
        stats["reporter_ratio_median"] = float(r["ratio"].median())
        stats["reporter_spearman_median"] = float(r["spearman_with_dapi_sharpness"].median())
    fin = d / "final_dapi_source.tsv"
    if fin.exists():
        f = pd.read_csv(fin, sep="\t").sort_values("ncc_with_final_dapi", ascending=False)
        top = f.iloc[0]
        text += (f"### Which DAPI the final image carries\n\nCorrelation of the final image's DAPI with each cycle file's DAPI "
                 f"(pyramid level 3): highest is cycle {int(top['cycle'])} {top['kind']} at {top['ncc_with_final_dapi']:.4f}; "
                 f"next {f.iloc[1]['ncc_with_final_dapi']:.4f} (cycle {int(f.iloc[1]['cycle'])} {f.iloc[1]['kind']}).\n\n")
        stats["final_dapi"] = f"cycle {int(top['cycle'])} {top['kind']} ({top['ncc_with_final_dapi']:.4f})"
    return text, stats


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("dirs", nargs="+", type=Path)
    ap.add_argument("--labels", nargs="*", default=None)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    labels = args.labels or [d.name for d in args.dirs]
    parts, allstats = [], {}
    for d, label in zip(args.dirs, labels):
        text, stats = describe(d, label)
        parts.append(text)
        allstats[label] = stats
    args.out.write_text("# Cycle-file evaluation: registration and DAPI sharpness\n\n"
                        "Shifts are in micrometres; ppm = micrometres per metre (100 ppm over 25 mm = 2.5 um). "
                        "Tiles fully covered by tissue only. Tables and maps, no verdict.\n\n" + "\n".join(parts))
    args.out.with_suffix(".json").write_text(json.dumps(allstats, indent=2))
    print(json.dumps(allstats, indent=2))
    print("->", args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
