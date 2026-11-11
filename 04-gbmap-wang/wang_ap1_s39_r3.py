"""Wang 2022 AP-1 test: is the AP-1 program elevated in recurrent tumors?

Rationale: Wang et al. claim recurrent GBMs shift to mesenchymal via AP-1.
Our S10 analysis found recurrent tumors enriched for EARLY modules with LATE
modules flat -- apparently contradicting the MES-shift. But S10 modules may
not capture AP-1 activity. This test asks directly whether AP-1 is up in
recurrent tumor cells. If AP-1 is up but S10-LATE is flat, the contradiction
is about how 'mesenchymal' is defined. If AP-1 is flat too, the MES-shift
claim is directly contradicted in this dataset.

Method mirrors wang_analysis.py: same 29 pairs, same tumor definition
(S10 total score > 0 & PTPRC < 0.5), same normalization. AP-1 score =
mean log-expression of AP-1 genes in tumor cells.
"""
import os, glob, gzip
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.io
from scipy.stats import wilcoxon, spearmanr

BASE = "/home/hatch/workspace/bioinfo_downloads/wang2022"
S10_JSON = "/home/hatch/workspace/bioinfo_downloads/gbm-point1/a-side/trajectory_metamodules_s10.json"
import json
s10 = json.load(open(S10_JSON))
EARLY = [1, 4, 6, 7, 8, 10, 11]
LATE = [2, 3, 5]

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
 "#39": ("GSM6199404", "GSM5319563"),
 "#40": ("GSM5319507", "GSM5319545"), "#41": ("GSM5319514", "GSM5319569"),
}

AP1_CORE = ["FOS", "FOSB", "FOSL1", "FOSL2", "JUN", "JUNB", "JUND", "ATF3"]
AP1_EXT = AP1_CORE + ["ATF4", "BATF", "BATF3", "MAF", "MAFB", "MAFF", "MAFG", "MAFK"]

def load_sample(gsm):
    pats = sorted(glob.glob(os.path.join(BASE, f"{gsm}_*_matrix.mtx.gz")))
    pats = [p for p in pats if "snATAC" not in p]
    b1 = [p for p in pats if "batch2" not in p]
    b2 = [p for p in pats if "batch2" in p]
    def load_one(mtx):
        feat = mtx.replace("_matrix.mtx.gz", "_features.tsv.gz")
        with gzip.open(feat, "rt") as fh:
            genes = [ln.strip().split("\t")[0] for ln in fh]
        m = scipy.io.mmread(mtx).tocsr()
        ad = sc.AnnData(m.T)
        ad.var_names = genes
        ad.var_names_make_unique()
        return ad
    ad = load_one(b1[0])
    if b2:
        ad2 = load_one(b2[0])
        ad = sc.concat([ad, ad2])
    return ad

rows = []
for pair in sorted(PAIRS, key=lambda p: (0, int(p[1:])) if p[1:].isdigit() else (1, p)):
    for tp, gsm in (("Primary", PAIRS[pair][0]), ("Recurrent", PAIRS[pair][1])):
        ad = load_sample(gsm)
        sc.pp.normalize_total(ad, target_sum=1e4)
        sc.pp.log1p(ad)
        # S10 module scores for tumor mask
        for mi, m in enumerate(s10):
            g = [x for x in m["genes"] if x in ad.var_names]
            if len(g) >= 10:
                sc.tl.score_genes(ad, g, score_name=f"MM{mi}")
        mod = np.array([ad.obs[f"MM{i}"].values for i in range(12)])
        total = mod.mean(axis=0)
        ptprc = ad[:, "PTPRC"].X.toarray().ravel() if "PTPRC" in ad.var_names else np.zeros(ad.n_obs)
        is_tumor = (total > 0) & (ptprc < 0.5)
        # AP-1 scores (tumor cells only)
        ap1 = {}
        for name, glist in (("AP1_core", AP1_CORE), ("AP1_ext", AP1_EXT)):
            g = [x for x in glist if x in ad.var_names]
            ap1[name] = float(ad[is_tumor, g].X.mean()) if g else np.nan
            ap1[name + "_n"] = len(g)
        # also AP-1 in all cells, and FOS/JUN individually in tumor
        ap1["AP1_core_all"] = float(ad[:, [x for x in AP1_CORE if x in ad.var_names]].X.mean())
        for g in ["FOS", "JUN", "FOSL1", "ATF3"]:
            ap1[f"{g}_tumor"] = float(ad[is_tumor, g].X.mean()) if g in ad.var_names else np.nan
        # TRAJ for cross-check
        traj = mod[LATE].mean(axis=0) - mod[EARLY].mean(axis=0)
        rows.append({"pair": pair, "tp": tp, "gsm": gsm, "n": ad.n_obs,
                     "n_tumor": int(is_tumor.sum()),
                     "traj_med_tumor": float(np.median(traj[is_tumor])),
                     **ap1})
        print(f"{pair} {tp} {gsm}: n_tumor={is_tumor.sum()} AP1_core(tumor)={ap1['AP1_core']:.3f} FOS={ap1['FOS_tumor']:.3f} JUN={ap1['JUN_tumor']:.3f}")
        del ad

tab = pd.DataFrame(rows)
tab.to_csv(os.path.join(BASE, "wang_analysis", "wang_ap1_s39_r3.csv"), index=False)
piv = tab.pivot(index="pair", columns="tp")
for col in ["AP1_core", "AP1_ext", "AP1_core_all", "FOS_tumor", "JUN_tumor", "FOSL1_tumor", "ATF3_tumor", "traj_med_tumor"]:
    d = piv[(col, "Recurrent")] - piv[(col, "Primary")]
    d = d.dropna()
    w, p = wilcoxon(d)
    print(f"{col}: median Recurrent-Primary diff={d.median():+.4f}, up={(d>0).sum()}/{(d<0).sum()} down, Wilcoxon p={p:.4f} (n={len(d)})")
# AP-1 vs TRAJ correlation (per-sample)
for col in ["AP1_core"]:
    r, p = spearmanr(tab[tab["tp"]=="Primary"][col], tab[tab["tp"]=="Primary"]["traj_med_tumor"])
    print(f"Primary AP1 vs TRAJ: rho={r:.3f} p={p:.4f}")
print("DONE")
