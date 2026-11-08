"""GBmap core analysis: trajectory validation + Point 3 ligand replication.

Design (see analysis_design.md):
- Exclude author == 'Zhao2020' (= GSE148842, PW032 only) for independence.
- X is log-normalized (float32, values ~0.3-6.5). Module score = mean of the
  module genes' log-expression (simple mean; no scanpy reference subtraction
  - noted as a method difference vs the canonical scanpy score_genes TRAJ).
- TRAJ = mean(LATE MM2/3/5) - mean(EARLY MM1/4/6/7/8/10/11).
- Process X in row chunks (CSR) to fit 7 GB RAM; only needed genes are kept.
"""
import h5py, json
import numpy as np
import scipy.sparse as sp
import pandas as pd
from scipy.stats import spearmanr, wilcoxon

H5 = "/home/hatch/workspace/bioinfo_downloads/gbmap/gbmap_core.h5ad"
S10_JSON = "/home/hatch/workspace/bioinfo_downloads/gbm-point1/a-side/trajectory_metamodules_s10.json"
OUT = "/home/hatch/workspace/bioinfo_downloads/gbmap/gbmap_analysis"
import os
os.makedirs(OUT, exist_ok=True)

EARLY = [1, 4, 6, 7, 8, 10, 11]
LATE = [2, 3, 5]
MYELOID_TYPES = {"macrophage", "microglial cell", "monocyte", "dendritic cell"}
TUMOR_TYPE = "malignant cell"

def read_col(obs, name):
    d = obs[name]
    if isinstance(d, h5py.Group):
        cats = [c.decode() if isinstance(c, bytes) else c for c in d["categories"][()]]
        codes = d["codes"][()]
        return np.array([cats[c] if c >= 0 else None for c in codes], dtype=object)
    v = d[()]
    if v.dtype.kind in ("S", "O"):
        return np.array([x.decode() if isinstance(x, bytes) else x for x in v], dtype=object)
    return v

with h5py.File(H5, "r") as f:
    obs = f["obs"]
    cell_type = read_col(obs, "cell_type")
    author = read_col(obs, "author")
    donor = read_col(obs, "donor_id")
    # gene symbols
    fn = f["var"]["feature_name"]
    cats = [c.decode() if isinstance(c, bytes) else str(c) for c in fn["categories"][:]]
    codes = fn["codes"][:]
    syms = np.array([cats[c] if c >= 0 else "" for c in codes])
    sym_to_idx = {}
    for i, s in enumerate(syms):
        if s and s not in sym_to_idx:
            sym_to_idx[s] = i
    n_cells = len(cell_type)
    n_genes = len(syms)
    print(f"cells={n_cells} genes={n_genes}")

    s10 = json.load(open(S10_JSON))
    # module gene indices (present genes only)
    mod_genes = {}
    for mi, m in enumerate(s10):
        idxs = [sym_to_idx[g] for g in m["genes"] if g in sym_to_idx]
        mod_genes[mi] = idxs
        print(f"MM{mi} {m['annotation'][:28]:28s} present {len(idxs)}/{len(m['genes'])}")
    lig_genes = {g: sym_to_idx[g] for g in ["VEGFA", "VCAN", "ANXA1"] if g in sym_to_idx}
    print("ligands:", lig_genes)

    # needed gene columns (union of module genes + ligands)
    needed = sorted(set(sum(mod_genes.values(), [])) | set(lig_genes.values()))
    needed_set = set(needed)
    print("needed genes:", len(needed))
    # weight matrix: needed gene -> module (1/n_present_in_module)
    W = np.zeros((len(needed), 12), dtype=np.float32)
    pos_of = {g: j for j, g in enumerate(needed)}
    for mi, idxs in mod_genes.items():
        for g in idxs:
            W[pos_of[g], mi] = 1.0 / len(idxs)
    # ligand positions within needed
    lig_pos = {g: pos_of[i] for g, i in lig_genes.items()}

    keep = author != "Zhao2020"
    print("excluding Zhao2020 (GSE148842):", int((~keep).sum()), "cells")
    is_tumor = (cell_type == TUMOR_TYPE) & keep
    is_myeloid = np.isin(cell_type, list(MYELOID_TYPES)) & keep
    print("tumor cells (excl GSE148842):", int(is_tumor.sum()))
    print("myeloid cells (excl GSE148842):", int(is_myeloid.sum()))

    # chunked CSR read
    X = f["X"]
    data = X["data"]
    indices = X["indices"]
    indptr = X["indptr"]
    CH = 20000
    mod_scores = np.zeros((n_cells, 12), dtype=np.float32)
    lig_vals = np.zeros((n_cells, len(lig_genes)), dtype=np.float32)
    lig_names = list(lig_genes.keys())
    for start in range(0, n_cells, CH):
        end = min(start + CH, n_cells)
        ip0, ip1 = int(indptr[start]), int(indptr[end])
        indptr_c = indptr[start:end + 1] - ip0
        data_c = data[ip0:ip1]
        idx_c = indices[ip0:ip1]
        Xc = sp.csr_matrix((data_c, idx_c, indptr_c), shape=(end - start, n_genes))
        Xsub = Xc[:, needed]  # sparse, only needed genes
        # module scores = Xsub @ W (mean over present module genes)
        mod_scores[start:end] = Xsub @ W
        # ligand values
        Xd = Xsub.toarray()
        for li, g in enumerate(lig_names):
            lig_vals[start:end, li] = Xd[:, lig_pos[g]]
        if (start // CH) % 5 == 0:
            print(f"  chunk {start}-{end}")
        del Xc, Xsub, Xd

# TRAJ
traj = mod_scores[:, LATE].mean(axis=1) - mod_scores[:, EARLY].mean(axis=1)
mm2 = mod_scores[:, 2]
np.savez(os.path.join(OUT, "gbmap_per_cell.npz"),
         traj=traj, mod_scores=mod_scores, lig_vals=lig_vals,
         lig_names=np.array(lig_names),
         is_tumor=is_tumor, is_myeloid=is_myeloid, donor=donor, author=author,
         cell_type=cell_type)

# --- Check 1: TRAJ vs PC1 of the 12 module scores (tumor cells) ---
from sklearn.decomposition import PCA
tm = is_tumor
ms_t = mod_scores[tm]
pca = PCA(n_components=2).fit(ms_t)
pc1 = pca.transform(ms_t)[:, 0]
rho_pc, p_pc = spearmanr(traj[tm], pc1)
# sign of PC1 is arbitrary; report |rho|
print(f"\n[Check 1] TRAJ vs PC1 (12 module scores, tumor n={tm.sum()}): "
      f"Spearman rho={rho_pc:.3f} (abs {abs(rho_pc):.3f}), p={p_pc:.2e}")
print(f"  PC1 explains {pca.explained_variance_ratio_[0]*100:.1f}% of module-score variance")

# --- Check 2: Point 3 replication, per-donor ---
donors = pd.DataFrame({
    "donor": donor, "traj": traj, "is_tumor": is_tumor, "is_myeloid": is_myeloid,
    "VEGFA": lig_vals[:, lig_names.index("VEGFA")],
    "VCAN": lig_vals[:, lig_names.index("VCAN")],
    "ANXA1": lig_vals[:, lig_names.index("ANXA1")],
})
rows = []
for d, g in donors.groupby("donor"):
    t = g[g["is_tumor"]]
    m = g[g["is_myeloid"]]
    if len(t) >= 50 and len(m) >= 20:
        rows.append({
            "donor": d, "n_tumor": len(t), "n_myeloid": len(m),
            "traj_mean": t["traj"].mean(),
            "VEGFA_myeloid": m["VEGFA"].mean(),
            "VCAN_myeloid": m["VCAN"].mean(),
            "ANXA1_myeloid": m["ANXA1"].mean(),
        })
per_donor = pd.DataFrame(rows)
per_donor.to_csv(os.path.join(OUT, "gbmap_per_donor.csv"), index=False)
print(f"\n[Check 2] donors with >=50 tumor & >=20 myeloid cells: {len(per_donor)}")
for lig in ["VEGFA", "VCAN", "ANXA1"]:
    r, p = spearmanr(per_donor["traj_mean"], per_donor[f"{lig}_myeloid"])
    print(f"  tumor TRAJ vs myeloid {lig}: rho={r:.3f}, p={p:.2e} (n={len(per_donor)})")

# --- Check 3: VEGFA in MM2-high vs MM2-low tumor cells ---
vegfa_t = lig_vals[tm, lig_names.index("VEGFA")]
mm2_t = mm2[tm]
hi = mm2_t >= np.median(mm2_t)
from scipy.stats import mannwhitneyu
u, p_mw = mannwhitneyu(vegfa_t[hi], vegfa_t[~hi])
print(f"\n[Check 3] tumor VEGFA, MM2-high (n={hi.sum()}) mean={vegfa_t[hi].mean():.3f} "
      f"vs MM2-low mean={vegfa_t[~hi].mean():.3f}, MW p={p_mw:.2e}")
print(f"  VEGFA+ fraction MM2-high={(vegfa_t[hi]>0).mean():.3f} vs low={(vegfa_t[~hi]>0).mean():.3f}")
# also VCAN/ANXA1 for context
for lig in ["VCAN", "ANXA1"]:
    v = lig_vals[tm, lig_names.index(lig)]
    print(f"  {lig}: MM2-high mean={v[hi].mean():.3f} vs low={v[~hi].mean():.3f}")
print("DONE")
