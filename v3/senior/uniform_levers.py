#!/usr/bin/env python3
"""uniform_levers.py -- what each setting of the UniFORM step does, marker by marker and donor by donor.

The maintainer, 2026-10-02: "UniFORM should be able to handle the different pixel intensity distributions, check if
there are levers that need to be adjusted."

The lane's UniFORM step (core/phenocycler/scripts/celltype/run_step2.py) has these settings:

  refusal      a marker whose median differs between acquisition groups by 2x the spread within a group is NOT
               aligned at all (usability.BARS['acq_ratio_flag']). The lane's own addition; not part of UniFORM.
  floor        cell values under 1.0 are dropped before the log histogram (UniFORM's log_transform_intensities)
  reference    the donor every other donor is shifted onto: UniFORM picks the histogram closest to the mean one
  counts       histograms enter the cross-correlation as raw counts (a donor with more cells weighs more in the mean)
  bins         1,024 bins over the cohort's log range
  method       cross-correlation of whole histograms (UniFORM's automatic registration); the alternative UniFORM
               offers is a landmark per donor

For every marker this measures, with each setting moved one at a time, what UniFORM's own registration gives:

  mode_sd      spread (SD, log units) of the donors' background modes after the shift -- the lane's 0.10 bar
  edge_sd      spread of the donors' lower half-maximum edge of the background peak after the shift: a second
               landmark of the same peak, not used by any setting, so it can judge all of them
  shape_r      correlation of each donor's shifted histogram with the reference's (median over donors)
  per donor    the residual of its background mode from the cohort median after the shift, and the share of its
               cells dropped by the floor

Nothing here reads a cell type, a frequency or a disease label. It writes tables; it installs nothing.

    senior_v3.sh --cd core/phenocycler python <this file> --out <results>/typing/uniform_levers [--markers ...]
"""
from __future__ import annotations

import argparse
import contextlib
import io
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

BINS = 1024
FLOOR = 1.0


def log_values(values: np.ndarray, floor: float) -> np.ndarray:
    v = values[np.isfinite(values) & (values >= floor)]
    return np.log(v)


def histograms(vals: dict, floor: float, bins: int):
    logs = {d: log_values(v, floor) for d, v in vals.items()}
    logs = {d: lv for d, lv in logs.items() if len(lv)}
    lo = min(float(lv.min()) for lv in logs.values()); hi = max(float(lv.max()) for lv in logs.values())
    counts = {d: np.histogram(lv, bins=bins, range=(lo, hi))[0].astype(float) for d, lv in logs.items()}
    return counts, lo, hi


def mode_bin(c: np.ndarray) -> int:
    return int(np.argmax(c[1:])) + 1                      # skip bin 0: the floor


def lower_edge_bin(c: np.ndarray) -> float:
    """Where the background peak rises through half its height, on its low side (interpolated bin position)."""
    k = mode_bin(c)
    smooth = np.convolve(c, np.ones(9) / 9.0, mode="same")
    half = smooth[k] / 2.0
    j = k
    while j > 1 and smooth[j] > half:
        j -= 1
    if smooth[j + 1] == smooth[j]:
        return float(j)
    return float(j + (half - smooth[j]) / (smooth[j + 1] - smooth[j]))


def register(reg, counts: dict, donors: list, bins: int, *, reference: str | None = None, density: bool = False, landmarks: list | None = None):
    mat = [counts[d] / counts[d].sum() if density else counts[d] for d in donors]
    hist = {"m": {"counts": mat, "bin_edges": [None] * len(donors)}}
    with contextlib.redirect_stdout(io.StringIO()):
        shift_map, refs, _ = reg.automatic_registration(hist, ["m"], ["m"], donors, [reference] if reference else None,
                                                        {"m": landmarks} if landmarks is not None else None, bins)
    import matplotlib.pyplot as plt  # noqa: PLC0415
    plt.close("all")
    return dict(zip(donors, shift_map["m"])), refs[0]


def judge(counts: dict, donors: list, shifts: dict, ref: str, lo: float, hi: float, bins: int) -> dict:
    binw = (hi - lo) / (bins - 1)
    mode = np.array([(mode_bin(counts[d]) - shifts[d]) * binw for d in donors])
    edge = np.array([(lower_edge_bin(counts[d]) - shifts[d]) * binw for d in donors])
    ref_c = counts[ref] / counts[ref].sum()
    shape = []
    for d in donors:
        c = counts[d] / counts[d].sum()
        s = int(shifts[d])
        moved = np.zeros_like(c)
        if s >= 0:
            moved[:bins - s] = c[s:]
        else:
            moved[-s:] = c[:bins + s]
        shape.append(float(np.corrcoef(moved, ref_c)[0, 1]))
    return {"mode_sd": float(mode.std()), "edge_sd": float(edge.std()), "shape_r": float(np.median(shape)),
            "mode_resid": dict(zip(donors, (mode - np.median(mode)).tolist())),
            "edge_resid": dict(zip(donors, (edge - np.median(edge)).tolist())),
            "shape_by_donor": dict(zip(donors, shape)),
            "scale": {d: float(np.exp(-shifts[d] * binw)) for d in donors}}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--markers", nargs="*", default=None)
    args = ap.parse_args(argv)

    sys.path.insert(0, str(Path.cwd() / "scripts" / "celltype"))
    import panel  # noqa: PLC0415
    import run_step2 as S2  # noqa: PLC0415
    import uniform_adapter as UA  # noqa: PLC0415
    reg, _, _ = UA._import_uniform()

    cohort = panel.cohort()
    donors_all = cohort["donor_id"].tolist()
    group = dict(zip(cohort["donor_id"], cohort["acq_group"]))
    usab = pd.read_parquet(panel.OUT / "marker_donor_usability.parquet").groupby("marker").agg(acq_ratio=("acq_ratio", "first"))
    markers = args.markers or sorted(usab.index)
    args.out.mkdir(parents=True, exist_ok=True)

    settings = {
        "as_run":        dict(floor=FLOOR, bins=BINS, density=False, reference=None),
        "density":       dict(floor=FLOOR, bins=BINS, density=True, reference=None),
        "floor_0.25":    dict(floor=0.25, bins=BINS, density=False, reference=None),
        "floor_0.25_density": dict(floor=0.25, bins=BINS, density=True, reference=None),
        "bins_512":      dict(floor=FLOOR, bins=512, density=False, reference=None),
        "bins_2048":     dict(floor=FLOOR, bins=2048, density=False, reference=None),
    }
    rows, per_donor = [], []
    for m in markers:
        vals = S2.marker_values(m, donors_all, {})
        donors = sorted(vals)
        if len(donors) < 2:
            continue
        dropped = {d: float(np.mean(~(np.isfinite(vals[d]) & (vals[d] >= FLOOR)))) for d in donors}
        zero = {d: float(np.mean(vals[d][np.isfinite(vals[d])] <= 0)) for d in donors}
        for name, s in settings.items():
            counts, lo, hi = histograms(vals, s["floor"], s["bins"])
            have = [d for d in donors if d in counts]
            shifts, ref = register(reg, counts, have, s["bins"], reference=s["reference"], density=s["density"])
            j = judge(counts, have, shifts, ref, lo, hi, s["bins"])
            rows.append({"marker": m, "setting": name, "n_donors": len(have), "reference": ref, "acq_ratio": float(usab.loc[m, "acq_ratio"]) if m in usab.index else np.nan,
                         "mode_sd": j["mode_sd"], "edge_sd": j["edge_sd"], "shape_r": j["shape_r"]})
            for d in have:
                per_donor.append({"marker": m, "setting": name, "donor_id": d, "group": group[d], "scale": j["scale"][d],
                                  "mode_resid": j["mode_resid"][d], "edge_resid": j["edge_resid"][d], "shape_r": j["shape_by_donor"][d],
                                  "dropped_by_floor_1": dropped[d], "zero_or_less": zero[d], "n_cells": int(len(vals[d]))})
        # a landmark per donor, the way UniFORM offers it: here the background mode itself (then mode_sd is 0 by
        # construction and only edge_sd and shape_r can judge it)
        counts, lo, hi = histograms(vals, FLOOR, BINS)
        have = [d for d in donors if d in counts]
        lm = [float(mode_bin(counts[d])) for d in have]
        shifts, ref = register(reg, counts, have, BINS, landmarks=lm)
        j = judge(counts, have, shifts, ref, lo, hi, BINS)
        rows.append({"marker": m, "setting": "landmark_mode", "n_donors": len(have), "reference": ref, "acq_ratio": float(usab.loc[m, "acq_ratio"]) if m in usab.index else np.nan,
                     "mode_sd": j["mode_sd"], "edge_sd": j["edge_sd"], "shape_r": j["shape_r"]})
        for d in have:
            per_donor.append({"marker": m, "setting": "landmark_mode", "donor_id": d, "group": group[d], "scale": j["scale"][d],
                              "mode_resid": j["mode_resid"][d], "edge_resid": j["edge_resid"][d], "shape_r": j["shape_by_donor"][d],
                              "dropped_by_floor_1": dropped[d], "zero_or_less": zero[d], "n_cells": int(len(vals[d]))})
        r = [x for x in rows if x["marker"] == m]
        print(f"{m:<18} " + " | ".join(f"{x['setting']}: mode {x['mode_sd']:.3f} edge {x['edge_sd']:.3f} r {x['shape_r']:.3f}" for x in r), flush=True)
        del vals
    pd.DataFrame(rows).to_csv(args.out / "levers_by_marker.tsv", sep="\t", index=False, float_format="%.5g")
    pd.DataFrame(per_donor).to_csv(args.out / "levers_by_marker_donor.tsv", sep="\t", index=False, float_format="%.5g")
    print(f"-> {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
