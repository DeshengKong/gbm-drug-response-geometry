"""Follow-up: T-cell exclusion test using the AUTHOR's mes2 score and the
discrete tier2 hypoxia state, from the saved per-tumor-cell table.
This avoids re-running the spatial queries (already done) and stays
memory-light. Primary: POST sections, per-section median split + Wilcoxon.
"""
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr, wilcoxon

OUT = "/home/hatch/workspace/bioinfo_downloads/gse296577/meylan_analysis"

tum = pd.read_csv(OUT + "/mm2_tcell_tumor_cells.csv.gz", index_col=0)
print("tumor cells:", tum.shape)
print("cols:", list(tum.columns))
post = tum[tum["timepoint"] == "post"].copy()
pre = tum[tum["timepoint"] == "pre"].copy()
print(f"post n={len(post)}, pre n={len(pre)}")

# ---------- A. Author mes2, POST ----------
print("\n===== A. Author mes2 (POST) =====")
med = post.groupby("orig.ident")["mes2"].transform("median")
post["mes2_hi"] = post["mes2"] > med
hi = post[post["mes2_hi"]]; lo = post[~post["mes2_hi"]]
u, p = mannwhitneyu(hi["dist_T"], lo["dist_T"])
print(f"pooled MW dist: hi (n={len(hi)}) median={hi['dist_T'].median():.1f} um "
      f"vs lo (n={len(lo)}) median={lo['dist_T'].median():.1f} um; p={p:.2e}")
print(f"  mean dist hi={hi['dist_T'].mean():.1f} vs lo={lo['dist_T'].mean():.1f}")
print(f"  frac<=50um hi={hi['has_T_50'].mean():.3f} vs lo={lo['has_T_50'].mean():.3f}")
rho, pr = spearmanr(post["mes2"], post["dist_T"])
print(f"  per-cell Spearman mes2 vs dist_T: rho={rho:.3f}, p={pr:.2e}")

rows = []
for sec, g in post.groupby("orig.ident"):
    m = g["mes2"].median()
    h = g[g["mes2"] > m]; l = g[g["mes2"] <= m]
    rows.append({"section": sec, "n_hi": len(h), "n_lo": len(l),
                 "dist_hi": h["dist_T"].median(), "dist_lo": l["dist_T"].median(),
                 "dist_diff": h["dist_T"].median() - l["dist_T"].median(),
                 "f50_hi": h["has_T_50"].mean(), "f50_lo": l["has_T_50"].mean(),
                 "f50_diff": h["has_T_50"].mean() - l["has_T_50"].mean()})
sec = pd.DataFrame(rows)
print(sec.to_string(index=False))
sec.to_csv(OUT + "/mm2_tcell_per_section_mes2_post.csv", index=False)
w, pw = wilcoxon(sec["dist_hi"], sec["dist_lo"])
print(f"per-section Wilcoxon dist (n={len(sec)}): median diff={sec['dist_diff'].median():.1f} um, "
      f"positive={int((sec['dist_diff']>0).sum())}/{len(sec)}, p={pw:.4f}")
w2, pw2 = wilcoxon(sec["f50_hi"], sec["f50_lo"])
print(f"per-section Wilcoxon frac50 (n={len(sec)}): median diff={sec['f50_diff'].median():.3f}, p={pw2:.4f}")
# tertile on mes2
rows3 = []
for sec_, g in post.groupby("orig.ident"):
    q1, q2 = g["mes2"].quantile([1/3, 2/3])
    h3 = g[g["mes2"] > q2]; l3 = g[g["mes2"] <= q1]
    rows3.append({"section": sec_, "dist_hi3": h3["dist_T"].median(),
                  "dist_lo3": l3["dist_T"].median(),
                  "f50_hi3": h3["has_T_50"].mean(), "f50_lo3": l3["has_T_50"].mean()})
t3 = pd.DataFrame(rows3)
print(t3.to_string(index=False))
w3, p3 = wilcoxon(t3["dist_hi3"], t3["dist_lo3"])
print(f"tertile Wilcoxon dist: p={p3:.4f}, median diff={(t3['dist_hi3']-t3['dist_lo3']).median():.1f} um")
w3b, p3b = wilcoxon(t3["f50_hi3"], t3["f50_lo3"])
print(f"tertile Wilcoxon frac50: p={p3b:.4f}, median diff={(t3['f50_hi3']-t3['f50_lo3']).median():.3f}")

# ---------- B. Discrete tier2 hypoxia state, POST ----------
print("\n===== B. tier2 MES-like 2 (hypoxia) vs other tumor states (POST) =====")
post["is_hyp"] = post["tier2"] == "Tumor MES-like 2 (hypoxia)"
hyp = post[post["is_hyp"]]; oth = post[~post["is_hyp"]]
u, p = mannwhitneyu(hyp["dist_T"], oth["dist_T"])
print(f"pooled MW dist: hypoxia (n={len(hyp)}) median={hyp['dist_T'].median():.1f} um "
      f"vs other (n={len(oth)}) median={oth['dist_T'].median():.1f} um; p={p:.2e}")
print(f"  mean dist hyp={hyp['dist_T'].mean():.1f} vs other={oth['dist_T'].mean():.1f}")
print(f"  frac<=50um hyp={hyp['has_T_50'].mean():.3f} vs other={oth['has_T_50'].mean():.3f}")
rows = []
for sec_, g in post.groupby("orig.ident"):
    h = g[g["is_hyp"]]; o = g[~g["is_hyp"]]
    if len(h) < 20 or len(o) < 20:
        continue
    rows.append({"section": sec_, "n_hyp": len(h), "n_oth": len(o),
                 "dist_hyp": h["dist_T"].median(), "dist_oth": o["dist_T"].median(),
                 "dist_diff": h["dist_T"].median() - o["dist_T"].median(),
                 "f50_hyp": h["has_T_50"].mean(), "f50_oth": o["has_T_50"].mean(),
                 "f50_diff": h["has_T_50"].mean() - o["has_T_50"].mean()})
secb = pd.DataFrame(rows)
print(secb.to_string(index=False))
secb.to_csv(OUT + "/mm2_tcell_per_section_tier2_post.csv", index=False)
wb, pb = wilcoxon(secb["dist_hyp"], secb["dist_oth"])
print(f"per-section Wilcoxon dist (n={len(secb)}): median diff={secb['dist_diff'].median():.1f} um, "
      f"positive={int((secb['dist_diff']>0).sum())}/{len(secb)}, p={pb:.4f}")
wb2, pb2 = wilcoxon(secb["f50_hyp"], secb["f50_oth"])
print(f"per-section Wilcoxon frac50 (n={len(secb)}): median diff={secb['f50_diff'].median():.3f}, p={pb2:.4f}")

# ---------- C. Author mes2, PRE (secondary) ----------
print("\n===== C. Author mes2 (PRE, secondary) =====")
medp = pre.groupby("orig.ident")["mes2"].transform("median")
pre["mes2_hi"] = pre["mes2"] > medp
hip = pre[pre["mes2_hi"]]; lop = pre[~pre["mes2_hi"]]
up, pp = mannwhitneyu(hip["dist_T"], lop["dist_T"])
print(f"pooled MW dist: hi (n={len(hip)}) median={hip['dist_T'].median():.1f} um "
      f"vs lo (n={len(lop)}) median={lop['dist_T'].median():.1f} um; p={pp:.2e}")
print(f"  frac<=50um hi={hip['has_T_50'].mean():.3f} vs lo={lop['has_T_50'].mean():.3f}")
rowsp = []
for sec_, g in pre.groupby("orig.ident"):
    m = g["mes2"].median()
    h = g[g["mes2"] > m]; l = g[g["mes2"] <= m]
    rowsp.append({"section": sec_, "dist_hi": h["dist_T"].median(),
                  "dist_lo": l["dist_T"].median(),
                  "dist_diff": h["dist_T"].median() - l["dist_T"].median(),
                  "f50_hi": h["has_T_50"].mean(), "f50_lo": l["has_T_50"].mean(),
                  "f50_diff": h["has_T_50"].mean() - l["has_T_50"].mean()})
secp = pd.DataFrame(rowsp)
print(secp.to_string(index=False))
wp, ppw = wilcoxon(secp["dist_hi"], secp["dist_lo"])
print(f"per-section Wilcoxon dist (n={len(secp)}): median diff={secp['dist_diff'].median():.1f} um, p={ppw:.4f}")
print("DONE")
