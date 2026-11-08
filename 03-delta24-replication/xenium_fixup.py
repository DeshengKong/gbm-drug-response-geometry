#!/usr/bin/env python3
"""GSE330794 Xenium fixup: re-summarize R-Q2 with mean-based metrics
(median is 0 in both groups under zero-inflation) and add a sensitive
T-cell definition (CD3D or CD3E count > 0) to cross-check fractions.
"""
import glob, os, re
import numpy as np
import pandas as pd
import h5py
from scipy import sparse, stats

D = "/home/hatch/workspace/bioinfo_downloads/gse330794"
OUT = os.path.join(D, "xenium_analysis")
MM2_NO_VEGFA = ["MET", "CAV1", "LOX", "ICAM1", "CXCL10", "CD44"]
T_MARKERS = ["CD3D", "CD3E", "CD8A", "CD4"]
TUMOR_MARKERS = ["OLIG2", "SOX2", "EGFR"]

def load_section(h5_path):
    with h5py.File(h5_path, "r") as f:
        names = [n.decode() for n in f["matrix/features/name"][:]]
        data = f["matrix/data"][:]
        indices = f["matrix/indices"][:]
        indptr = f["matrix/indptr"][:]
        shape = tuple(f["matrix/shape"][:])
    csc = sparse.csc_matrix((data, indices, indptr), shape=shape)
    return names, csc.tocsr()

def log_norm(csr):
    tot = np.asarray(csr.sum(axis=0)).ravel()
    tot[tot == 0] = 1
    norm = (csr @ sparse.diags(10000.0 / tot)).tocsr()
    norm.data = np.log1p(norm.data)
    return norm

def score(norm, names, genes):
    idx = [names.index(g) for g in genes if g in names]
    if not idx:
        return np.zeros(norm.shape[1])
    return np.asarray(norm[idx, :].mean(axis=0)).ravel()

rows = []
for h5 in sorted(glob.glob(os.path.join(D, "GSM97313*_section_*_cell_feature_matrix.h5"))):
    m = re.search(r"(P\d+)-(pre|post)_section_(\d+)", os.path.basename(h5))
    pat, tp, sec = m.group(1), m.group(2), int(m.group(3))
    names, csr = load_section(h5)
    n = csr.shape[1]
    norm = log_norm(csr)
    # sensitive T: any of CD3D/CD3E/CD8A/CD4 count > 0 (raw counts)
    t_idx = [names.index(g) for g in T_MARKERS]
    t_any = np.asarray(csr[t_idx, :].sum(axis=0)).ravel() > 0
    # tumor: strict rule as before
    sc_tumor = score(norm, names, TUMOR_MARKERS)
    sc_t = score(norm, names, T_MARKERS)
    sc_m = score(norm, names, ["CD68", "CD163", "CSF1R", "AIF1"])
    sc_b = score(norm, names, ["MS4A1", "CD79A"])
    sc_nk = score(norm, names, ["NCAM1", "KLRD1"])
    sc_s = score(norm, names, ["ACTA2", "FN1", "LUM"])
    S = np.vstack([sc_t, sc_m, sc_b, sc_nk, sc_tumor, sc_s])
    order = np.argsort(-S, axis=0)
    best = S[order[0], np.arange(n)]
    second = S[order[1], np.arange(n)]
    types = np.array(["T", "Myeloid", "B", "NK", "Tumor", "Stromal"])[order[0]]
    assign = np.where((best > 0.5) & ((best - second) >= 0.2), types, "Other")
    is_tumor = assign == "Tumor"
    vegfa_idx = names.index("VEGFA")
    vegfa = np.asarray(norm[vegfa_idx, :].todense()).ravel() > 0
    mm2 = score(norm, names, MM2_NO_VEGFA)
    # Q2 mean-based metrics within tumor cells
    if is_tumor.sum() > 0 and vegfa[is_tumor].sum() > 10 and (~vegfa[is_tumor]).sum() > 10:
        a = mm2[is_tumor & vegfa]
        b = mm2[is_tumor & ~vegfa]
        rows.append(dict(
            patient=pat, timepoint=tp, section=sec,
            n_cells=n, n_tumor=int(is_tumor.sum()),
            frac_T_strict=float((assign == "T").mean()),
            frac_T_sensitive=float(t_any.mean()),
            vegfa_pos_tumor=float(vegfa[is_tumor].mean()),
            mm2_mean_pos=float(a.mean()), mm2_mean_neg=float(b.mean()),
            mm2_mean_diff=float(a.mean() - b.mean()),
            mm2_pos_frac_pos=float((a > 0).mean()),
            mm2_pos_frac_neg=float((b > 0).mean()),
            mw_p=float(stats.mannwhitneyu(a, b).pvalue)))
    print(f"{pat}-{tp} s{sec}: T strict={float((assign=='T').mean()):.4f} "
          f"T sens={float(t_any.mean()):.4f} tumor={float(is_tumor.mean()):.3f}")

df = pd.DataFrame(rows)
df.to_csv(os.path.join(OUT, "xenium_vegfa_mm2_fixup.csv"), index=False)
print("\n=== Q2 (mean-based) ===")
print(df[["patient", "timepoint", "section", "vegfa_pos_tumor",
          "mm2_mean_diff", "mm2_pos_frac_pos", "mm2_pos_frac_neg", "mw_p"]]
      .to_string(index=False))
print("\nmedian mm2 mean-diff:", float(df["mm2_mean_diff"].median()),
      "| positive in", int((df["mm2_mean_diff"] > 0).sum()), "/", len(df))
w = stats.wilcoxon(df["mm2_mean_diff"])
print("Wilcoxon on per-section mean diffs: p =", w.pvalue)
print("median MM2-positive frac: VEGFA+",
      float(df["mm2_pos_frac_pos"].median()), "vs VEGFA-",
      float(df["mm2_pos_frac_neg"].median()))
print("\n=== T fraction (sensitive) paired P17/P20 ===")
for pat in ["P17", "P20"]:
    sub = df[df.patient == pat].groupby("timepoint")[["frac_T_strict", "frac_T_sensitive"]].mean()
    print(pat)
    print(sub.to_string())
    print("post - pre:", (sub.loc["post"] - sub.loc["pre"]).to_dict())
print("\nDONE")
