#!/usr/bin/env python3
"""Build an isolated all-raw QuPath expression tree for the typing calibration comparison.

The registry-selected expression tree uses REDSEA whole-cell measurements for the markers that
take REDSEA correction and raw QuPath nucleus means for the nuclear markers. This producer builds
the declared all-raw comparator: replace every non-nuclear marker by its raw QuPath whole-cell
mean, preserve the same raw nuclear measurements, and retain identical object/QC/geometry rows.
It is never written over the registry-selected inputs.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq


ROOT = Path("/home/smith6jt/IO60panc2nd/Islet-Explorer-Senior-v3")
PKG = ROOT / "core/phenocycler"
DEFAULT_CELLS = PKG / "data/cells"
DEFAULT_EXPRESSION = PKG / "data/expression"
NUCLEAR = {"DAPI", "FOXP3", "Ki67", "PCNA", "SOX2", "TCF_1", "TOX", "TP63"}


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cells", type=Path, default=DEFAULT_CELLS)
    ap.add_argument("--expression", type=Path, default=DEFAULT_EXPRESSION)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--donors", nargs="+", required=True)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    manifest = []
    for donor in args.donors:
        ep = args.expression / f"donor_id={donor}/data.parquet"
        cp = args.cells / f"donor_id={donor}/data_0.parquet"
        if not ep.is_file():
            # A small number of partitions use the alternate name data.parquet.
            cp2 = args.cells / f"donor_id={donor}/data.parquet"
            if cp2.is_file():
                cp = cp2
        if not ep.is_file() or not cp.is_file():
            raise FileNotFoundError(f"{donor}: missing expression or cell export: {ep}, {cp}")

        expr_file = pq.ParquetFile(ep)
        cell_file = pq.ParquetFile(cp)
        expression_names = expr_file.schema_arrow.names
        cells_names = set(cell_file.schema_arrow.names)
        metadata = {"donor_id", "object_id", "qc_analysis_eligible", "qc_estimation_eligible",
                    "image", "X_centroid", "Y_centroid", "cell_region", "islet_num",
                    "in_alignment_block"}
        markers = [c for c in expression_names if c not in metadata]
        def raw_name(marker: str) -> str:
            if marker in cells_names:
                return marker
            # The raw QuPath exporter uses a stricter sanitization than the registry for the
            # phospho-Histone channel: `Histone_H3_(p_Ser28)` -> `Histone_H3_p_Ser28`.
            return re.sub(r"[^0-9A-Za-z]+", "_", marker).strip("_")

        raw_columns = {m: raw_name(m) for m in markers if m not in NUCLEAR}
        replaced = [m for m, c in raw_columns.items() if c in cells_names]
        absent = [m for m, c in raw_columns.items() if c not in cells_names]
        if absent:
            raise KeyError(f"{donor}: raw whole-cell QuPath values absent for {absent}")

        # Read only the joined measurement columns, retaining expression row order as the source
        # of eligibility/geometry. The sorted-key merge is explicit and validated one-to-one.
        expr = expr_file.read().to_pandas()
        raw_source_columns = list(dict.fromkeys(raw_columns[m] for m in replaced))
        raw = cell_file.read(columns=["object_id", *raw_source_columns]).to_pandas()
        if expr.object_id.duplicated().any() or raw.object_id.duplicated().any():
            raise ValueError(f"{donor}: duplicate object_id in one of the input partitions")
        aligned = expr[["object_id"]].merge(raw, on="object_id", how="left", sort=False,
                                             validate="one_to_one")
        if len(aligned) != len(expr) or aligned[raw_source_columns].isna().all(axis=None):
            raise ValueError(f"{donor}: raw-to-expression join failed or returned all missing")
        missing_ids = int(aligned[raw_source_columns].isna().all(axis=1).sum())
        if missing_ids:
            raise ValueError(f"{donor}: {missing_ids:,} expression objects lack raw measurements")
        for m in replaced:
            expr[m] = aligned[raw_columns[m]].to_numpy()

        target = args.out / f"donor_id={donor}/data.parquet"
        target.parent.mkdir(parents=True, exist_ok=True)
        table = pa.Table.from_pandas(expr, preserve_index=False)
        metadata_bytes = dict(table.schema.metadata or {})
        metadata_bytes[b"typing_input"] = b"all_raw_qupath_cell_mean_except_raw_nucleus"
        metadata_bytes[b"raw_cells_sha256"] = digest(cp).encode()
        metadata_bytes[b"registry_expression_sha256"] = digest(ep).encode()
        table = table.replace_schema_metadata(metadata_bytes)
        pq.write_table(table, target, compression="zstd", row_group_size=500_000)
        manifest.append({"donor_id": donor, "rows": len(expr), "raw_cell_markers": replaced,
                         "retained_raw_nuclear_markers": sorted(set(markers) & NUCLEAR),
                         "missing_raw_markers": absent, "missing_join_rows": missing_ids,
                         "cells_source": str(cp), "cells_sha256": digest(cp),
                         "registry_expression_source": str(ep), "expression_sha256": digest(ep),
                         "raw_output": str(target)})
        print(f"{donor}: {len(expr):,} rows; replaced {len(replaced)} whole-cell markers; "
              f"retained {len(set(markers) & NUCLEAR)} raw nuclear measurements", flush=True)

    out_manifest = args.out / "raw_input_manifest.json"
    out_manifest.write_text(json.dumps({"method": "all_raw_qupath_cell_mean_except_raw_nucleus",
                                        "donors": manifest}, indent=2) + "\n")
    print(f"-> {out_manifest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
