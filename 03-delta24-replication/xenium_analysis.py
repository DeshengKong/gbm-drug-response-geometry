#!/usr/bin/env python3
"""GSE330794 Xenium analysis (pre-registered, replication_design.md).
Marker-score cell typing per section; R-Q1 (T fraction, P17/P20 descriptive),
R-Q2 (VEGFA vs MM2-hypoxia score in tumor cells), R-Q4 (myeloid, descriptive).
MM2 score uses the 6 non-VEGFA MM2 genes on the panel
(MET, CAV1, LOX, ICAM1, CXCL10, CD44) to avoid circularity with VEGFA.
"""
import glob, os, re
import numpy as np
import pandas as pd
import h5py
from scipy import sparse, stats

D = "/home/hatch/workspace/bioinfo_downloads/gse330794"
OUT = os.path.join(D, "xenium_analysis")
os.makedirs(OUT, exist_ok=True)

MARKERS = {
    "T": ["CD3D", "CD3E", "CD8A", "CD4"],
    "Myeloid": ["CD68", "CD163", "CSF1R", "AIF1"],
    "B": ["MS4A1", "CD79A"],
    "NK": ["NCAM1", "KLRD1"],
    "Tumor": ["OLIG2", "SOX2", "EGFR"],
    "Stromal": ["ACTA2", "FN1", "LUM"],
}
MM2_NO_VEGFA = ["MET", "CAV1", "LOX", "ICAM1", "CXCL10", "CD44"]

def load_section(h5_path):
    with h5py.File(h5_path, "r") as f:
        names = [n.decode() for n in f["matrix/features/name"][:]]
        data = f["matrix/data"][:]
        indices = f["matrix/indices"][:]
        indptr = f["matrix/indptr"][:]
        shape = tuple(f["matrix/shape"][:])
    # CSC (genes x cells), indptr length = n_cells+1 -> convert to CSR (gene-major)
    csc = sparse.csc_matrix((data, indices, indptr), shape=shape)
    csr = csc.tocsr()
    return names, csr

def log_norm(csr):
    tot = np.asarray(csr.sum(axis=0)).ravel()
    tot[tot == 0] = 1
    # scale each cell (column) to 10k then log1p; keep sparse
    norm = (csr @ sparse.diags(10000.0 / tot)).tocsr()
    norm.data = np.log1p(norm.data)
    return norm

def score(norm, names, genes):
    idx = [names.index(g) for g in genes if g in names]
    if not idx:
        return np.zeros(norm.shape[1])
    sub = norm[idx, :]
    return np.asarray(sub.mean(axis=0)).ravel()

rows = []
tumor_vegfa = []   # per-section VEGFA vs MM2 comparison
for h5 in sorted(glob.glob(os.path.join(D, "GSM97313*_section_*_cell_feature_matrix.h5"))):
    base = os.path.basename(h5)
    m = re.search(r"(P\d+)-(pre|post)_section_(\d+)", base)
    pat, tp, sec = m.group(1), m.group(2), int(m.group(3))
    names, csr = load_section(h5)
    n_cells = csr.shape[1]
    norm = log_norm(csr)
    sc = {k: score(norm, names, v) for k, v in MARKERS.items()}
    # assign: highest score wins if > 0.5 and beats runner-up by >= 0.2
    S = np.vstack([sc[k] for k in MARKERS])          # (6, n)
    order = np.argsort(-S, axis=0)
    best = S[order[0], np.arange(n_cells)]
    second = S[order[1], np.arange(n_cells)]
    types = np.array(list(MARKERS.keys()))[order[0]]
    assign = np.where((best > 0.5) & ((best - second) >= 0.2), types, "Other")
    frac = {k: float((assign == k).mean()) for k in list(MARKERS.keys()) + ["Other"]}
    # VEGFA positivity overall and in tumor cells
    vegfa_idx = names.index("VEGFA")
    vegfa = np.asarray(norm[vegfa_idx, :].todense()).ravel() > 0
    mm2 = score(norm, names, MM2_NO_VEGFA)
    is_tumor = assign == "Tumor"
    n_tumor = int(is_tumor.sum())
    vegfa_in_tumor = float(vegfa[is_tumor].mean()) if n_tumor else np.nan
    # VEGFA+ vs VEGFA- tumor cells: median MM2 difference
    if n_tumor and vegfa[is_tumor].sum() > 10 and (~vegfa[is_tumor]).sum() > 10:
        mm2_pos = np.median(mm2[is_tumor & vegfa])
        mm2_neg = np.median(mm2[is_tumor & ~vegfa])
        w = stats.mannwhitneyu(mm2[is_tumor & vegfa], mm2[is_tumor & ~vegfa])
        tumor_vegfa.append(dict(patient=pat, timepoint=tp, section=sec,
                                n_tumor=n_tumor, vegfa_pos_frac=vegfa_in_tumor,
                                mm2_median_pos=mm2_pos, mm2_median_neg=mm2_neg,
                                mm2_diff=mm2_pos - mm2_neg, mw_p=w.pvalue))
    rows.append(dict(patient=pat, timepoint=tp, section=sec, n_cells=n_cells,
                     **{f"frac_{k}": v for k, v in frac.items()},
                     vegfa_pos_all=float(vegfa.mean()),
                     vegfa_pos_tumor=vegfa_in_tumor,
                     n_tumor=n_tumor))
    print(f"{pat}-{tp} s{sec}: n={n_cells} T={frac['T']:.3f} "
          f"Myeloid={frac['Myeloid']:.3f} Tumor={frac['Tumor']:.3f} "
          f"Other={frac['Other']:.3f} VEGFA(tumor)={vegfa_in_tumor:.3f}")

tab = pd.DataFrame(rows)
tab.to_csv(os.path.join(OUT, "xenium_fractions.csv"), index=False)
tv = pd.DataFrame(tumor_vegfa)
tv.to_csv(os.path.join(OUT, "xenium_vegfa_mm2.csv"), index=False)

print("\n=== R-Q1 / R-Q4: paired sections (P17, P20) ===")
for pat in ["P17", "P20"]:
    sub = tab[tab.patient == pat].groupby("timepoint")[["frac_T", "frac_Myeloid", "frac_Tumor"]].mean()
    print(pat)
    print(sub.to_string())
    print("post - pre:", (sub.loc["post"] - sub.loc["pre"]).to_dict())

print("\n=== R-Q1 cross-section: T fraction by timepoint (all sections) ===")
print(tab.groupby("timepoint")[["frac_T", "frac_Myeloid", "frac_Tumor"]].agg(["mean", "median", "count"]).to_string())

print("\n=== R-Q2: VEGFA+ vs VEGFA- tumor MM2 (per-section, 15 sections) ===")
print(tv[["patient", "timepoint", "section", "vegfa_pos_frac", "mm2_diff", "mw_p"]].to_string(index=False))
print("\nmedian VEGFA+ fraction in tumor cells:", float(tv["vegfa_pos_frac"].median()))
print("MM2 diff (pos - neg): median", float(tv["mm2_diff"].median()),
      "positive in", int((tv["mm2_diff"] > 0).sum()), "/", len(tv))
w = stats.wilcoxon(tv["mm2_diff"])
print("Wilcoxon on per-section MM2 diffs: p =", w.pvalue)
print("\nDONE")
