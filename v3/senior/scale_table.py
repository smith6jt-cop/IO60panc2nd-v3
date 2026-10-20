#!/usr/bin/env python3
"""scale_table.py -- the cohort scale on the v3 aligned values, beside the frozen one.

For the Gate-2 scale ruling (plan: "the frozen cosine_scale.parquet -- recompute with the existing
run_step3.build_scale on the new aligned values, keep, or type under both"). This gives the numbers the ruling
is made on and no verdict:

  frozen   data/development/broadtype/scratch/cosine_scale.parquet of the v2 checkout -- what the typing lane reads
           today through subpanel_io.load_scale(): floor (aligned background mode) and top (p99.9) per marker,
           built 2026-09-05 over 20 donors
  v3       the same builder with the same arguments (run_step3.build_scale, 200,000 cells per donor, seed 0, no
           DAPI ratio) over the v3 aligned values of the 25 donors' pancreas blocks -- Tissue plus same-block
           node cells, the sample UniFORM itself used (ruling of 2026-10-02)
  v3, Tissue cells only
           the same again with the same-block node cells left out of the SCALE sample (they stay in the
           alignment). Shown because a lymph node's cells are bright on lymphocyte markers and the top is the
           99.9th percentile of the pooled sample: this column says how much of a top's change is theirs.
  v3, Tissue cells only, the frozen file's own donors
           the same over the donors the frozen file was built on (read from its metadata): what remains
           is the change the re-segmentation and the redrawn Tissue outlines made on the same donors.

A floor and a top are in the aligned UNIT, and the unit of a marker is its reference donor's: UniFORM shifts every
donor onto one of them. So the table carries each marker's reference donor in both runs -- where it changed, the
two scales are in different units and their ratio says nothing about the tissue.

Run inside the v3 worktree:
    scripts/v3/senior/senior_v3.sh --cd core/phenocycler python <this file> \\
        --frozen <v2 checkout>/data/development/broadtype/scratch/cosine_scale.parquet \\
        --frozen-shifts <v2 checkout>/data/development/celltype/shifts.json \\
        --candidate ../../data/development/celltype/scale_v3_candidate.parquet --out <results>/scale_table.tsv

The candidate is written beside the aligned values and is NOT installed: the typing lane keeps reading the frozen
file until the maintainer rules.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

N_PER_DONOR, SEED, RATIO = 200_000, 0, False        # the frozen file's own arguments (its parquet metadata)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--frozen", type=Path, required=True)
    ap.add_argument("--frozen-shifts", type=Path, default=None)
    ap.add_argument("--candidate", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()

    sys.path.insert(0, str(Path.cwd() / "scripts" / "celltype"))
    import panel  # noqa: PLC0415
    import run_step3 as S3  # noqa: PLC0415

    donors = panel.cohort()["donor_id"].tolist()
    universe = sorted(panel.markers(include_variant=True))
    new = S3.build_scale(donors, universe, n_per_donor=N_PER_DONOR, seed=SEED, ratio=RATIO)
    args.candidate.parent.mkdir(parents=True, exist_ok=True)
    new.reset_index().to_parquet(args.candidate, index=False)
    # the same builder over Tissue cells only: for this table, the block rule is narrowed to the Tissue region --
    # every cell outside a pLN region. Since pass J a cell in an islet object reads core / cluster / single, not
    # `tissue`; it is a Tissue cell all the same, so the test is "not pLN" (2026-10-02; it was `== "tissue"`).
    block_population = S3.AB.population
    try:
        S3.AB.population = lambda donor, object_id, cell_region, root=None: np.asarray(cell_region, dtype=object) != "pLN"
        tissue = S3.build_scale(donors, universe, n_per_donor=N_PER_DONOR, seed=SEED, ratio=RATIO)
    finally:
        S3.AB.population = block_population
    tissue.reset_index().to_parquet(args.candidate.with_name(args.candidate.stem + "_tissue_only.parquet"), index=False)
    import pyarrow.parquet as pq  # noqa: PLC0415
    frozen_meta = json.loads((pq.read_schema(args.frozen).metadata or {}).get(b"broadtype_cosine_scale", b"{}"))
    frozen_donors = [d for d in donors if str(d).lstrip("0") in {str(x).lstrip("0") for x in frozen_meta.get("donors", [])}]
    same_donors = None
    if frozen_donors:
        try:
            S3.AB.population = lambda donor, object_id, cell_region, root=None: np.asarray(cell_region, dtype=object) != "pLN"
            same_donors = S3.build_scale(frozen_donors, universe, n_per_donor=N_PER_DONOR, seed=SEED, ratio=RATIO)
        finally:
            S3.AB.population = block_population

    old = pd.read_parquet(args.frozen).set_index("marker")
    shifts_new = json.loads((panel.OUT / "shifts.json").read_text())["shifts"]
    shifts_old = json.loads(args.frozen_shifts.read_text())["shifts"] if args.frozen_shifts and args.frozen_shifts.is_file() else {}

    rows = []
    for m in sorted(set(old.index) | set(new.index)):
        o = old.loc[m] if m in old.index else None
        n = new.loc[m] if m in new.index else None
        sn, so = shifts_new.get(m, {}), shifts_old.get(m, {})
        ratio = lambda a, b: (float(a) / float(b)) if a is not None and b is not None and np.isfinite(a) and np.isfinite(b) and b > 0 else np.nan  # noqa: E731
        rows.append({
            "marker": m,
            "frozen_floor": None if o is None else float(o["floor"]), "frozen_top": None if o is None else float(o["top"]),
            "v3_floor": None if n is None else float(n["floor"]), "v3_top": None if n is None else float(n["top"]),
            "v3_tissue_only_floor": float(tissue.loc[m, "floor"]) if m in tissue.index else None,
            "v3_tissue_only_top": float(tissue.loc[m, "top"]) if m in tissue.index else None,
            "top_tissue_only_over_frozen": ratio(tissue.loc[m, "top"] if m in tissue.index else None, None if o is None else o["top"]),
            "v3_tissue_only_frozen_donors_top": float(same_donors.loc[m, "top"]) if same_donors is not None and m in same_donors.index else None,
            "top_tissue_only_frozen_donors_over_frozen": ratio(same_donors.loc[m, "top"] if same_donors is not None and m in same_donors.index else None,
                                                               None if o is None else o["top"]),
            "floor_v3_over_frozen": ratio(None if n is None else n["floor"], None if o is None else o["floor"]),
            "top_v3_over_frozen": ratio(None if n is None else n["top"], None if o is None else o["top"]),
            "span_v3_over_frozen": ratio(None if n is None else n["span"], None if o is None else o["span"]),
            "reference_donor_frozen": so.get("reference_donor", ""), "reference_donor_v3": sn.get("reference_donor", ""),
            "same_reference_donor": bool(so) and str(so.get("reference_donor", "")).lstrip("0") == str(sn.get("reference_donor", "")).lstrip("0"),
            "v3_alignment": sn.get("provenance", ""), "v3_donors_aligned": sn.get("n_donors", ""),
            "v3_mode_sd_before": sn.get("mode_sd_before", ""), "v3_mode_sd_after": sn.get("mode_sd_after", ""),
        })
    table = pd.DataFrame(rows)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(args.out, sep="\t", index=False, float_format="%.4g")

    both = table.dropna(subset=["floor_v3_over_frozen", "top_v3_over_frozen"])
    same = both[both["same_reference_donor"]]
    print(f"markers: frozen {len(old)}, v3 {len(new)}, in both {len(both)}; only in v3: {sorted(set(new.index) - set(old.index))}; "
          f"only frozen: {sorted(set(old.index) - set(new.index))}")
    print(f"same reference donor in both runs: {len(same)} of {len(both)} markers")
    for name, sub in (("all markers in both", both), ("same reference donor", same), ("reference donor changed", both[~both["same_reference_donor"]])):
        if len(sub):
            q = lambda c: np.percentile(sub[c], [10, 50, 90])  # noqa: E731
            f, t = q("floor_v3_over_frozen"), q("top_v3_over_frozen")
            print(f"  {name} ({len(sub)}): floor v3/frozen p10 {f[0]:.2f}, median {f[1]:.2f}, p90 {f[2]:.2f} | top v3/frozen p10 {t[0]:.2f}, median {t[1]:.2f}, p90 {t[2]:.2f}")
    print(f"donors: v3 {len(donors)}; the frozen file's own: {len(frozen_donors)} of {len(frozen_meta.get('donors', []))} are in v3")
    print(f"-> {args.out}\n-> {args.candidate} (not installed)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
