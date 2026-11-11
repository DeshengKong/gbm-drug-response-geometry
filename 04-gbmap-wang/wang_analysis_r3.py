#!/usr/bin/env python3
"""Wang 2022 primary-vs-recurrent trajectory analysis (pre-registered).
Per-sample: load MTX (batch1+batch2 combined), log-normalize, score the 12 S10
modules with scanpy score_genes, TRAJ = mean(LATE)-mean(EARLY).
Tumor = tumor-score > 0.2 and tumor-score > PTPRC. Myeloid = CD68/CD163 high.
T cells = CD3D/CD3E high. Per-sample summaries only (no per-cell dump).
"""
import gzip, os, re
import numpy as np
import pandas as pd
import scipy.io, scipy.sparse
import scanpy as sc
from scipy import stats

D = "/home/hatch/workspace/bioinfo_downloads/wang2022"
OUT = os.path.join(D, "wang_analysis")
os.makedirs(OUT, exist_ok=True)
import json
with open("/home/hatch/workspace/bioinfo_downloads/gbm-point1/a-side/trajectory_metamodules_s10.json") as f:
    MODS = json.load(f)
MM_ANNOT = {m["metamodule"]: m.get("annotation", "") for m in MODS}
EARLY = [1, 4, 6, 7, 8, 10, 11]
LATE = [2, 3, 5]
rng = np.random.default_rng(0)

# --- pair list (from SOFT parse) ---
PAIRS = {
 "#1": ("GSM5319575", "GSM5319534"),
 "#9": ("GSM5319503", "GSM5319536"),
 "#10": ("GSM5319525", "GSM6199409"), "#11": ("GSM5319526", "GSM6199405"),
 "#12": ("GSM5319509", "GSM5319565"), "#13": ("GSM5319541", "GSM5319561"),
 "#14": ("GSM5319522", "GSM5319576"), "#15": ("GSM5319527", "GSM5319555"),
 "#16": ("GSM5319530", "GSM5319568"), "#17": ("GSM5319528", "GSM5319557"),
 "#19": ("GSM5319504", "GSM5319539"), "#20": ("GSM5319531", "GSM5319535"),
 "#21": ("GSM5319524", "GSM5319552"), "#22": ("GSM5319523", "GSM5319572"),
 "#23": ("GSM5319508", "GSM5319573"), "#24": ("GSM5319513", "GSM5319546"),
 "#25": ("GSM5319518", "GSM5319548"), "#26": ("GSM5319517", "GSM5319547"),
 "#27": ("GSM5319519", "GSM5319549"), "#28": ("GSM5319506", "GSM5319562"),
 "#29": ("GSM5319521", "GSM5319551"), "#30": ("GSM5319511", "GSM5319543"),
 "#31": ("GSM5319510", "GSM5319560"), "#32": ("GSM5319520", "GSM5319550"),
 "#35": ("GSM5319505", "GSM5319542"),
 "#37": ("GSM5319559", "GSM5319532"),
 "#39": ("GSM6199404", "GSM5319574"),
 "#40": ("GSM5319507", "GSM5319545"), "#41": ("GSM5319514", "GSM5319569"),
}  # 29 pairs; #33 dropped (GSM5319529 MTX truncated, see design doc); #2/#8 lack snRNA on one side; #13 corrected; #37 primary/recurrent swapped per Diaz 2026-10-09 (SF9798 primary, SF9494 recurrent; GEO progression annotations reversed); NA pair excluded
SAMPLES = {}  # gsm -> (pair, timepoint)
for pair, (pri, rec) in PAIRS.items():
    SAMPLES[pri] = (pair, "Primary")
    SAMPLES[rec] = (pair, "Recurrent")

def load_sample(gsm):
    """Load a sample's MTX (batch1 + batch2 if present), dedupe gene symbols."""
    pref = None
    for f in os.listdir(D):
        if f.startswith(gsm + "_") and f.endswith("_matrix.mtx.gz") \
                and "batch2" not in f and "snATAC" not in f:
            pref = f[: -len("_matrix.mtx.gz")]
            break
    if pref is None:
        raise FileNotFoundError(gsm)
    mats = []
    for suffix in ["", "_batch2"]:
        mtx = os.path.join(D, pref + suffix + "_matrix.mtx.gz")
        feat = os.path.join(D, pref + suffix + "_features.tsv.gz")
        if not os.path.exists(mtx):
            continue
        m = scipy.io.mmread(mtx).tocsr()  # genes x cells
        with gzip.open(feat, "rt") as fh:
            genes = [ln.strip().split("\t")[0] for ln in fh]
        # dedupe gene symbols (first occurrence)
        seen, keep_idx, keep_genes = set(), [], []
        for i, g in enumerate(genes):
            if g not in seen:
                seen.add(g); keep_idx.append(i); keep_genes.append(g)
        mats.append((m[keep_idx, :], keep_genes))
    # combine batches (cells concatenated); gene lists should match
    base_genes = mats[0][1]
    assert all(g == base_genes for _, g in mats), f"gene mismatch {gsm}"
    X = scipy.sparse.hstack([m for m, _ in mats]).tocsr()  # genes x cells
    return X.T.tocsr(), base_genes  # cells x genes

def mean_expr(adata, genes):
    idx = [adata.var_names.get_loc(g) for g in genes if g in adata.var_names]
    if not idx:
        return np.zeros(adata.n_obs)
    return np.asarray(adata.X[:, idx].mean(axis=1)).ravel()

rows = []
CSV_PATH = os.path.join(OUT, "wang_per_sample_r3.csv")
done_gsm = set()
if os.path.exists(CSV_PATH):
    done_gsm = set(pd.read_csv(CSV_PATH)["gsm"])
    rows = pd.read_csv(CSV_PATH).to_dict("records")
    print(f"resuming: {len(done_gsm)} samples already done")
for gsm, (pair, tp) in sorted(SAMPLES.items(), key=lambda kv: kv[1][0]):
    if gsm in done_gsm:
        continue
    X, genes = load_sample(gsm)
    ad = sc.AnnData(X)
    ad.var_names = genes
    ad.obs_names = [f"{gsm}_{i}" for i in range(ad.n_obs)]
    sc.pp.normalize_total(ad, target_sum=1e4)
    sc.pp.log1p(ad)
    # S10 module scores
    for m in MODS:
        gl = [g for g in m["genes"] if g in ad.var_names]
        if len(gl) < 5:
            ad.obs[f"MM{m['metamodule']}"] = np.nan
            continue
        sc.tl.score_genes(ad, gl, score_name=f"MM{m['metamodule']}",
                          random_state=0, use_raw=False)
    mm = np.vstack([ad.obs[f"MM{i}"].values for i in range(12)]).T  # cells x 12
    traj = mm[:, LATE].mean(axis=1) - mm[:, EARLY].mean(axis=1)
    # tumor = high total S10 tumor-program expression (sum of 12 modules)
    # + no immune marker; SOX2/OLIG2/EGFR/PDGFRA are too broad (98% of cells)
    total_prog = mm.mean(axis=1)
    ptprc = mean_expr(ad, ["PTPRC"])
    is_tumor = (total_prog > 0) & (ptprc < 0.5)
    tumor_sc = mean_expr(ad, ["SOX2", "OLIG2", "EGFR", "PDGFRA"])
    myel_sc = mean_expr(ad, ["CD68", "CD163"])
    is_myel = (myel_sc > 0.5) & (myel_sc > tumor_sc)
    t_sc = mean_expr(ad, ["CD3D", "CD3E"])
    is_t = (t_sc > 0.5)
    lig = {g: float(mean_expr(ad, [g])[is_myel].mean()) if is_myel.sum() > 10 else np.nan
           for g in ["VEGFA", "VCAN", "ANXA1"]}
    row = dict(pair=pair, timepoint=tp, gsm=gsm, n_cells=int(ad.n_obs),
               n_tumor=int(is_tumor.sum()), frac_tumor=float(is_tumor.mean()),
               n_myeloid=int(is_myel.sum()), frac_t=float(is_t.mean()),
               traj_median_tumor=float(np.median(traj[is_tumor])) if is_tumor.sum() else np.nan,
               traj_mean_tumor=float(np.mean(traj[is_tumor])) if is_tumor.sum() else np.nan,
               traj_median_all=float(np.median(traj)),
               **{f"MM{i}_tumor": float(np.nanmean(mm[is_tumor, i])) if is_tumor.sum() else np.nan
                  for i in range(12)},
               **{f"lig_{k}": v for k, v in lig.items()})
    rows.append(row)
    pd.DataFrame(rows).to_csv(CSV_PATH, index=False)  # incremental save
    print(f"{pair} {tp} {gsm}: n={ad.n_obs} tumor={is_tumor.mean():.2f} "
          f"TRAJ_med(tumor)={row['traj_median_tumor']:.3f} T={is_t.mean():.3f}")
    del ad, X

tab = pd.DataFrame(rows)

# --- main test: paired Wilcoxon on median TRAJ (tumor) ---
piv = tab.pivot(index="pair", columns="timepoint", values="traj_median_tumor").dropna()
print("\npaired pairs (median TRAJ):", len(piv))
d = piv["Recurrent"] - piv["Primary"]
print("median Recurrent-Primary TRAJ:", float(d.median()))
print("up (recurrent higher):", int((d > 0).sum()), "down:", int((d < 0).sum()))
w = stats.wilcoxon(piv["Recurrent"], piv["Primary"])
print("Wilcoxon p =", w.pvalue)
# permutation null
pr, rc = piv["Primary"].values, piv["Recurrent"].values
obs = float(np.median(rc - pr))
nulls = []
for _ in range(10000):
    flip = rng.random(len(pr)) < 0.5
    dd = np.where(flip, pr - rc, rc - pr)
    nulls.append(np.median(dd))
nulls = np.array(nulls)
print(f"perm: obs {obs:+.4f}, null 95% CI [{np.percentile(nulls,2.5):+.4f}, {np.percentile(nulls,97.5):+.4f}], "
      f"empirical p = {(np.abs(nulls) >= abs(obs)).mean():.4f}")
# mean TRAJ version
pivm = tab.pivot(index="pair", columns="timepoint", values="traj_mean_tumor").dropna()
dm = pivm["Recurrent"] - pivm["Primary"]
wm = stats.wilcoxon(pivm["Recurrent"], pivm["Primary"])
print("mean TRAJ: median diff", float(dm.median()), "Wilcoxon p =", wm.pvalue)
# per-module paired
print("\nper-module paired (tumor, Recurrent-Primary median diff, Wilcoxon p):")
for i in range(12):
    pm = tab.pivot(index="pair", columns="timepoint", values=f"MM{i}_tumor").dropna()
    dd = pm["Recurrent"] - pm["Primary"]
    ww = stats.wilcoxon(pm["Recurrent"], pm["Primary"])
    print(f"MM{i} {MM_ANNOT.get(i,'')}: median diff {float(dd.median()):+.4f}, p={ww.pvalue:.4g}")
# T-cell fraction paired
pt = tab.pivot(index="pair", columns="timepoint", values="frac_t").dropna()
dt = pt["Recurrent"] - pt["Primary"]
wt = stats.wilcoxon(pt["Recurrent"], pt["Primary"])
print("\nT-cell fraction: median diff", float(dt.median()), "up", int((dt>0).sum()), "down", int((dt<0).sum()), "Wilcoxon p =", wt.pvalue)
# Point 3 replication: myeloid ligand vs tumor TRAJ (per-sample)
for lig in ["VEGFA", "VCAN", "ANXA1"]:
    sub = tab.dropna(subset=[f"lig_{lig}", "traj_median_tumor"])
    r, p = stats.spearmanr(sub[f"lig_{lig}"], sub["traj_median_tumor"])
    print(f"ligand {lig} vs tumor TRAJ (n={len(sub)}): rho={r:.3f} p={p:.4g}")
print("\nDONE")
