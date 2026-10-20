#!/usr/bin/env python3
"""uniform_bright_end.py -- does UniFORM's one factor per donor bring the bright cells into line, or only the background?

UniFORM moves each donor's whole distribution of a marker by ONE factor, found from the histogram as a whole -- in
practice from the background peak, because that is where the cells are. That is the right correction when a donor
differs from the others by a gain: background and bright cells move together. It is not when a donor's background
alone sits lower or higher (an offset): the factor that lines up the background then moves the bright cells apart.

For every marker and donor this reads the values the alignment read (the pancreas block's eligible cells,
``run_step2.marker_values``) and the factor the alignment applied (``shifts.json``), and writes, before and after
the factor:

  mode       the background peak of the log histogram (the lane's own landmark)
  p99.9      the 99.9th percentile of the donor's cells; p99.99 beside it
  top/mode   p99.9 over the background peak -- a ratio the factor cannot change

and, per marker, the four 70-channel donors beside the other 21: the same quantities as medians and ranges; and over
all donors of a marker, how far their bright ends lie apart before and after the factor (SD of the log), with the
slope of the bright end on the background: near 1 the donors differ by a gain and the factor brings both ends into
line, near 0 only their backgrounds differ and the factor moves the bright ends apart.

For one marker (--level-for, default CD8) it also writes what the typing's level means in stored counts for each
donor: the immune definitions ask a scaled value of 0.10, i.e. floor + 0.10 x span of the cohort scale in aligned
units; divided by the donor's factor that is a count in its own image.

Nothing here reads a cell type or a disease label, and no percentile is used as a threshold: they describe where a
donor's brightest cells are. A donor with few positive cells has a lower p99.9 for that reason alone, so one donor
proves nothing; the table is for reading a GROUP of donors against another.

    senior_v3.sh --cd core/phenocycler python <this file> --shifts <shifts.json> --out <results>/typing/uniform_levers
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from uniform_levers import BINS, FLOOR, histograms, mode_bin  # noqa: E402

FOUR = ("0112", "6347", "6533", "6580")       # the 70-channel same-slide donors (cohort table: acquisition group)
QUANTILES = (50.0, 99.0, 99.9, 99.99)


def donor_rows(marker: str, vals: dict, scale: dict, provenance: str) -> list[dict]:
    counts, lo, hi = histograms(vals, FLOOR, BINS)
    width = (hi - lo) / BINS
    rows = []
    for d, v in sorted(vals.items()):
        if d not in counts:
            continue
        v = v[np.isfinite(v)]
        q = np.percentile(v, QUANTILES)
        mode = float(np.exp(lo + (mode_bin(counts[d]) + 0.5) * width))
        f = float(scale.get(d, np.nan))
        rows.append({"marker": marker, "donor_id": d, "four": d in FOUR, "provenance": provenance, "n_cells": int(v.size),
                     "under_1": float(np.mean(v < FLOOR)), "factor": f, "mode": mode,
                     "p50": q[0], "p99": q[1], "p99.9": q[2], "p99.99": q[3],
                     "mode_aligned": mode * f, "p99.9_aligned": q[2] * f, "p99.99_aligned": q[3] * f,
                     "top_over_mode": q[2] / mode if mode > 0 else np.nan})
    return rows


def rng(x: pd.Series, fmt: str = "{:.3g}") -> str:
    x = x.dropna()
    return "" if x.empty else f"{fmt.format(x.min())}-{fmt.format(x.max())}"


def group_table(per: pd.DataFrame) -> pd.DataFrame:
    out = []
    for m, g in per.groupby("marker", sort=True):
        a, b = g[g["four"]], g[~g["four"]]
        if a.empty or b.empty:
            continue
        med = lambda col: (float(a[col].median()), float(b[col].median()))  # noqa: E731
        (ma, mb), (ta, tb), (fa, fb) = med("mode"), med("p99.9"), med("factor")
        (maa, mba), (taa, tba) = med("mode_aligned"), med("p99.9_aligned")
        (ua, ub) = med("p99.99_aligned")
        out.append({
            "marker": m, "provenance": g["provenance"].iloc[0], "n_four": len(a), "n_others": len(b),
            "factor_four": rng(a["factor"], "{:.2f}"), "factor_others": rng(b["factor"], "{:.2f}"),
            "factor_four_over_others": fa / fb,
            "mode_four_over_others_before": ma / mb, "mode_four_over_others_after": maa / mba,
            "top_four_over_others_before": ta / tb, "top_four_over_others_after": taa / tba,
            "top9999_four_over_others_after": ua / ub,
            "top_over_mode_four": rng(a["top_over_mode"], "{:.0f}"), "top_over_mode_others": rng(b["top_over_mode"], "{:.0f}"),
            "top_over_mode_four_over_others": float(a["top_over_mode"].median() / b["top_over_mode"].median()),
            "four_tops_all_above_others": bool(a["p99.9_aligned"].min() > b["p99.9_aligned"].max()),
            "four_tops_all_below_others": bool(a["p99.9_aligned"].max() < b["p99.9_aligned"].min()),
            "four_tops_before_all_above": bool(a["p99.9"].min() > b["p99.9"].max()),
            "four_tops_before_all_below": bool(a["p99.9"].max() < b["p99.9"].min()),
        })
    return pd.DataFrame(out)


def spread_table(per: pd.DataFrame, min_donors: int = 8) -> pd.DataFrame:
    """Per marker, over its donors: the spread of the bright end before and after the factor."""
    out = []
    for m, g in per.groupby("marker", sort=True):
        if len(g) < min_donors:
            continue
        bg, top, top_a = np.log(g["mode"]), np.log(g["p99.99"]), np.log(g["p99.99_aligned"])
        varies = float(bg.std(ddof=1)) > 0
        out.append({
            "marker": m, "provenance": g["provenance"].iloc[0], "n_donors": len(g),
            "background_fold": float(np.exp(bg.max() - bg.min())),
            "top_fold_before": float(np.exp(top.max() - top.min())), "top_fold_after": float(np.exp(top_a.max() - top_a.min())),
            "top_sd_before": float(top.std(ddof=1)), "top_sd_after": float(top_a.std(ddof=1)),
            "top_sd_after_over_before": float(top_a.std(ddof=1) / top.std(ddof=1)) if top.std(ddof=1) > 0 else np.nan,
            "slope_top_on_background": float(np.polyfit(bg, top, 1)[0]) if varies else np.nan,
            "r_top_background": float(np.corrcoef(bg, top)[0, 1]) if varies else np.nan,
        })
    return pd.DataFrame(out)


def level_table(per: pd.DataFrame, marker: str, floor: float, span: float, level: float) -> pd.DataFrame:
    """What `level` of the cohort scale is in each donor's stored counts, as aligned and if the donor were not shifted."""
    g = per[per["marker"] == marker].copy()
    aligned = floor + level * span
    g["count_at_level"] = aligned / g["factor"]
    g["count_at_level_over_background"] = g["count_at_level"] - g["mode"]
    g["count_at_level_share_of_p99.99"] = g["count_at_level"] / g["p99.99"]
    g["unshifted_over_background"] = aligned - g["mode"]
    return g[["marker", "donor_id", "four", "factor", "mode", "p99.99", "count_at_level", "count_at_level_over_background",
              "count_at_level_share_of_p99.99", "unshifted_over_background"]].sort_values(["four", "factor"])


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--shifts", type=Path, required=True, help="the alignment's shifts.json")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--markers", nargs="*", default=None)
    ap.add_argument("--scale", type=Path, default=None, help="the cohort scale (floor, span per marker) for --level-for")
    ap.add_argument("--level-for", nargs="*", default=["CD8"], help="markers to write the typing's level in stored counts for")
    ap.add_argument("--level", type=float, default=0.10, help="the scaled value the immune definitions ask (subpanel.IMMUNE_REQUIRE_MIN)")
    args = ap.parse_args(argv)

    sys.path.insert(0, str(Path.cwd() / "scripts" / "celltype"))
    import panel  # noqa: PLC0415
    import run_step2 as S2  # noqa: PLC0415

    record = json.loads(args.shifts.read_text())
    donors_all = panel.cohort()["donor_id"].tolist()
    markers = args.markers or sorted(record["shifts"])
    rows = []
    for m in markers:
        entry = record["shifts"][m]
        vals = S2.marker_values(m, donors_all, {})
        if len(vals) < 2:
            continue
        rows += donor_rows(m, vals, entry.get("scale") or {}, entry.get("provenance", ""))
        print(f"{m:<18} {len(vals)} donors", flush=True)
        del vals
    per = pd.DataFrame(rows)
    args.out.mkdir(parents=True, exist_ok=True)
    per.to_csv(args.out / "bright_end_by_marker_donor.tsv", sep="\t", index=False, float_format="%.5g")
    group_table(per).to_csv(args.out / "bright_end_four_vs_21.tsv", sep="\t", index=False, float_format="%.4g")
    spread_table(per).to_csv(args.out / "bright_end_spread_by_marker.tsv", sep="\t", index=False, float_format="%.4g")
    if args.scale is not None:
        scale = pd.read_parquet(args.scale).set_index("marker")
        levels = [level_table(per, m, float(scale.loc[m, "floor"]), float(scale.loc[m, "span"]), args.level)
                  for m in args.level_for if m in scale.index and (per["marker"] == m).any()]
        if levels:
            pd.concat(levels).to_csv(args.out / "level_in_stored_counts.tsv", sep="\t", index=False, float_format="%.4g")
    print(f"-> {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
