#!/usr/bin/env python3
"""eval_6457_figure.py -- one page from the cycle-file evaluation: the slide in question beside a control.

Reads only the tables eval_6457_cycles.py / eval_6457_focus.py / eval_6457_report.py wrote.

  A  median DAPI shift per cycle (tiles fully covered by tissue), log scale:
       raw vs cycle 0 raw            how far the slide really moved between cycles
       intermediate vs its own raw   what the instrument's registration applied
     drawn for both slides, so the control shows what an ordinary run looks like
  B  DAPI contrast per raw cycle, median and 5th percentile of tiles, both slides
  C  DAPI correlation of each raw cycle with the reference cycle's raw DAPI:
     both slides against cycle 0, and (if measured) the slide in question against cycle 1

    eval_6457_figure.py --case <dir> --control <dir> [--case-rawref <dir>] --out fig.png
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

import eval_6457_report as rep  # noqa: E402


def summary(d: Path) -> pd.DataFrame:
    f = d / "cycle_summary_interior.tsv"
    if not f.exists():
        meta = json.loads((d / "files.json").read_text())
        reg = rep.interior(pd.read_parquet(d / "registration_tiles.parquet"), pd.read_parquet(d / "cycle0_dapi_tiles.parquet"))
        rep.per_cycle(reg, int(meta["tile_px"])).to_csv(f, sep="\t", index=False)
    return pd.read_csv(f, sep="\t")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--case", type=Path, required=True)
    ap.add_argument("--control", type=Path, required=True)
    ap.add_argument("--case-rawref", type=Path, default=None)
    ap.add_argument("--case-label", default="6457")
    ap.add_argument("--control-label", default="6539 (control)")
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()

    case, control = summary(args.case), summary(args.control)
    fig, axes = plt.subplots(1, 3, figsize=(19, 5.6))
    colors = {"case": "#b2182b", "control": "#2166ac"}

    ax = axes[0]
    for label, s, key, comp, style in (
            (f"{args.case_label}: applied by the instrument's registration", case, "case", "int_vs_raw", "-"),
            (f"{args.case_label}: real movement (raw vs cycle 0 raw)", case, "case", "raw_vs_raw0", ":"),
            (f"{args.control_label}: applied by the instrument's registration", control, "control", "int_vs_raw", "-"),
            (f"{args.control_label}: real movement (raw vs cycle 0 raw)", control, "control", "raw_vs_raw0", ":")):
        d = s[s["comparison"] == comp].sort_values("cycle")
        ax.plot(d["cycle"], d["shift_um_median"].clip(lower=0.01), style, marker="o", ms=4, color=colors[key], label=label)
        if comp == "int_vs_raw":
            ax.plot(d["cycle"], d["shift_um_p95"].clip(lower=0.01), style, lw=0.8, alpha=0.5, color=colors[key])
    ax.set_yscale("log")
    ax.axhline(0.5078, color="grey", lw=0.8, ls="--")
    ax.text(0.2, 0.5078 * 1.08, "one pixel", color="grey", fontsize=8)
    ax.set_xlabel("cycle")
    ax.set_ylabel("DAPI shift, um (median over tissue tiles; thin line: 95th percentile)")
    ax.set_title("A  How far each cycle moved, and what the instrument applied")
    ax.legend(fontsize=7.5, loc="upper left")

    ax = axes[1]
    for label, d, key in ((args.case_label, args.case, "case"), (args.control_label, args.control, "control")):
        f = pd.read_csv(d / "focus" / "dapi_focus_by_cycle.tsv", sep="\t").sort_values("cycle")
        ax.plot(f["cycle"], f["contrast_median"], "-", marker="o", ms=4, color=colors[key], label=f"{label}: median")
        ax.plot(f["cycle"], f["contrast_p5"], ":", marker=".", ms=3, color=colors[key], label=f"{label}: 5th percentile of tiles")
    ax.set_ylim(0, 1)
    ax.set_xlabel("cycle (0 is the one whose DAPI the final image carries)")
    ax.set_ylabel("DAPI contrast, (p95 - p5) / (p95 + p5) per tile")
    ax.set_title("B  DAPI contrast in each raw cycle")
    ax.legend(fontsize=7.5, loc="lower right")

    ax = axes[2]
    for label, s, key, style in ((f"{args.case_label} vs its cycle 0", case, "case", "-"),
                                 (f"{args.control_label} vs its cycle 0", control, "control", "-")):
        d = s[s["comparison"] == "raw_vs_raw0"].sort_values("cycle")
        ax.plot(d["cycle"], d["ncc_median"], style, marker="o", ms=4, color=colors[key], label=label)
    if args.case_rawref is not None:
        meta = json.loads((args.case_rawref / "files.json").read_text())
        r = summary(args.case_rawref)
        d = r[r["comparison"] == "raw_vs_raw0"].sort_values("cycle")
        ax.plot(d["cycle"], d["ncc_median"], "--", marker="s", ms=4, color="#ef8a62",
                label=f"{args.case_label} vs its cycle {meta.get('raw_reference_cycle')}")
    ax.set_ylim(0, 1.02)
    ax.set_xlabel("cycle")
    ax.set_ylabel("correlation of raw DAPI with the reference cycle's (median over tiles)")
    ax.set_title("C  Does each cycle's DAPI look like the reference cycle's?")
    ax.legend(fontsize=7.5, loc="lower right")

    fig.suptitle("Cycle files of the Fusion run: registration and DAPI, tiles of 260 um fully covered by tissue", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(args.out, dpi=120)
    print("->", args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
