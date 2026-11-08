#!/usr/bin/env python3
"""Inspect the GSE296577 Seurat RDS metadata with pyreadr (no R on this VM).
Prints: object structure, assays, meta.data columns + levels, cell counts
per patient/timepoint, and whether any survival/response field exists.
Usage: inspect_rds.py  (run after the RDS file is complete)
"""
import json
import os

import pyreadr

RDS = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                   "GSE296577_GBM_oHSV_seurat_Full_20260107.rds")
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "rds_meta_summary.json")

print("file size:", os.path.getsize(RDS), "expected: 3843442210")
res = pyreadr.read_r(RDS)
print("top-level keys:", list(res.keys()))
summary = {"keys": list(res.keys())}

for key, obj in res.items():
    print(f"\n=== {key} === type={type(obj).__name__}")
    if hasattr(obj, "shape"):
        print("  shape:", obj.shape)
        summary[key] = {"type": str(type(obj)), "shape": list(obj.shape)}
    if hasattr(obj, "columns"):
        cols = list(obj.columns)
        print("  columns:", cols)
        summary[key] = {"columns": cols}
        # value counts for low-cardinality categorical columns
        for c in cols:
            try:
                vc = obj[c].value_counts(dropna=False)
                if len(vc) <= 30:
                    print(f"    {c}: {vc.to_dict()}")
                    summary[key][c] = {str(k): int(v) for k, v in vc.items()}
            except Exception as e:
                print(f"    {c}: (unreadable: {e})")

with open(OUT, "w") as f:
    json.dump(summary, f, indent=2, default=str)
print("\nsummary written to", OUT)
