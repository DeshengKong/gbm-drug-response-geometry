"""Meylan MM2-hypoxia vs T-cell spatial exclusion analysis.

- MM2 score = mean of log1p(counts) for the 4 MM2 genes available on the
  panel (MET, ICAM1, CXCL10, CD44). CAV1 and LOX are NOT on the panel,
  so the score is a 4-gene version (validated against the author's mes2).
- Tumor cells: tier1 == 'Tumor cells'. T cells: tier1 == 'T cells'.
- Spatial distances computed WITHIN each orig.ident section (x, y in um)
  with scipy cKDTree. Primary focus: POST sections (avoids pre/post
  composition confounding); pre reported as a secondary check.
- MM2 split: per-section median (controls for section-level MM2/T-density
  differences). Tertiles as a sensitivity check.
"""
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import mannwhitneyu, spearmanr, wilcoxon

BASE = "/home/hatch/workspace/bioinfo_downloads/gse296577"
OUT = BASE + "/meylan_analysis"

MM2_GENES = ["MET", "ICAM1", "CXCL10", "CD44"]

print("loading meta ...")
meta = pd.read_csv(BASE + "/meylan_meta.csv.gz", index_col=0)
print("meta shape:", meta.shape)
print("meta cols:", [c for c in ["tier1", "tier2", "timepoint", "orig.ident", "x", "y", "mes2", "tier1_T"] if c in meta.columns])

print("loading MM2 genes (long) ...")
long = pd.read_csv(OUT + "/mm2_genes_long.csv.gz")
print("long rows:", len(long), "genes:", sorted(long["gene"].unique()))
print("counts per gene:\n", long.groupby("gene")["count"].agg(["count", "sum", "max"]))

# pivot to per-cell mean log1p
long["logc"] = np.log1p(long["count"])
mm2 = long.groupby(["cell", "gene"])["logc"].mean().unstack(fill_value=0.0)
# ensure all 4 genes present as columns
for g in MM2_GENES:
    if g not in mm2.columns:
        mm2[g] = 0.0
mm2_score = mm2[MM2_GENES].mean(axis=1)
mm2_score.name = "MM2_4gene"
print("cells with any MM2 gene:", len(mm2_score))

# join to meta (index = cell barcode)
meta = meta.join(mm2_score, how="left")
meta["MM2_4gene"] = meta["MM2_4gene"].fillna(0.0)
print("joined. MM2 describe (tumor cells):")
tum_all = meta[meta["tier1"] == "Tumor cells"]
print(tum_all["MM2_4gene"].describe())
print("mes2 describe (tumor cells):")
print(tum_all["mes2"].describe())

# --- Validation: 4-gene MM2 vs author's mes2 (tumor cells) ---
rho_v, p_v = spearmanr(tum_all["MM2_4gene"], tum_all["mes2"])
print(f"\n[Validation] 4-gene MM2 vs author mes2 (tumor, n={len(tum_all)}): "
      f"Spearman rho={rho_v:.3f}, p={p_v:.2e}")
# MM2 by tier2 tumor state
print("\nMM2_4gene mean by tier2 tumor state:")
print(tum_all.groupby("tier2")["MM2_4gene"].agg(["mean", "median", "count"]))
print("\nmes2 mean by tier2 tumor state:")
print(tum_all.groupby("tier2")["mes2"].agg(["mean", "median", "count"]))
# MM2-high enrichment in MES-like 2 (hypoxia)
mes2_state = tum_all["tier2"] == "Tumor MES-like 2 (hypoxia)"
print(f"\nMM2_4gene: MES-like 2 mean={tum_all.loc[mes2_state,'MM2_4gene'].mean():.3f} "
      f"vs other tumor states={tum_all.loc[~mes2_state,'MM2_4gene'].mean():.3f}")

# --- Spatial: nearest-T distance per tumor cell, within each section ---
meta["dist_T"] = np.nan
meta["has_T_50"] = np.nan
for sec, g in meta.groupby("orig.ident"):
    t_idx = g.index[g["tier1"] == "Tumor cells"]
    tc_idx = g.index[g["tier1"] == "T cells"]
    if len(tc_idx) == 0 or len(t_idx) == 0:
        continue
    tree = cKDTree(g.loc[tc_idx, ["x", "y"]].values)
    d, _ = tree.query(g.loc[t_idx, ["x", "y"]].values, k=1)
    meta.loc[t_idx, "dist_T"] = d
    meta.loc[t_idx, "has_T_50"] = (d <= 50).astype(float)
print("\nspatial done. tumor cells with dist_T:", meta["dist_T"].notna().sum())

# save per-tumor-cell table (post + pre) for reuse
tum = meta[meta["tier1"] == "Tumor cells"].copy()
tum_out = tum[["orig.ident", "timepoint", "tier2", "x", "y", "mes2",
               "MM2_4gene", "dist_T", "has_T_50"]]
tum_out.to_csv(OUT + "/mm2_tcell_tumor_cells.csv.gz", index=True)
print("saved mm2_tcell_tumor_cells.csv.gz:", tum_out.shape)

def section_stats(df, label):
    """Per-section median-split MM2 high/low comparison."""
    rows = []
    for sec, g in df.groupby("orig.ident"):
        if len(g) < 100:
            continue
        med = g["MM2_4gene"].median()
        hi = g[g["MM2_4gene"] > med]
        lo = g[g["MM2_4gene"] <= med]
        if len(hi) < 20 or len(lo) < 20:
            continue
        rows.append({
            "section": sec,
            "timepoint": g["timepoint"].iloc[0],
            "n_hi": len(hi), "n_lo": len(lo),
            "dist_hi_med": hi["dist_T"].median(),
            "dist_lo_med": lo["dist_T"].median(),
            "dist_diff": hi["dist_T"].median() - lo["dist_T"].median(),
            "frac50_hi": hi["has_T_50"].mean(),
            "frac50_lo": lo["has_T_50"].mean(),
            "frac50_diff": hi["has_T_50"].mean() - lo["has_T_50"].mean(),
            "mm2_hi_med": hi["MM2_4gene"].median(),
            "mm2_lo_med": lo["MM2_4gene"].median(),
        })
    return pd.DataFrame(rows)

# ---------- Primary: POST sections ----------
post = tum[tum["timepoint"] == "post"].copy()
print(f"\n===== POST tumor cells: {len(post)} =====")
# pooled Mann-Whitney (per-section median split)
post_med = post.groupby("orig.ident")["MM2_4gene"].transform("median")
post["MM2_hi"] = post["MM2_4gene"] > post_med
hi = post[post["MM2_hi"]]; lo = post[~post["MM2_hi"]]
u, p_mw = mannwhitneyu(hi["dist_T"], lo["dist_T"])
print(f"[a] POST pooled MW: MM2-high (n={len(hi)}) median dist={hi['dist_T'].median():.1f} um "
      f"vs MM2-low (n={len(lo)}) median={lo['dist_T'].median():.1f} um; "
      f"MW p={p_mw:.2e}")
print(f"    mean dist hi={hi['dist_T'].mean():.1f} vs lo={lo['dist_T'].mean():.1f}")
u2, p_mw2 = mannwhitneyu(hi["has_T_50"], lo["has_T_50"])
print(f"[b] POST pooled MW (has T<=50um): hi frac={hi['has_T_50'].mean():.3f} "
      f"vs lo frac={lo['has_T_50'].mean():.3f}; MW p={p_mw2:.2e}")

# per-section paired Wilcoxon (post sections, n=8)
sec = section_stats(post, "post")
print("\nper-section (POST) table:")
print(sec.to_string(index=False))
sec.to_csv(OUT + "/mm2_tcell_per_section_post.csv", index=False)
if len(sec) >= 3:
    w, p_w = wilcoxon(sec["dist_hi_med"], sec["dist_lo_med"])
    print(f"\n[a] POST per-section Wilcoxon (n={len(sec)} sections): "
          f"median dist_diff={sec['dist_diff'].median():.1f} um, "
          f"all-positive={int((sec['dist_diff']>0).sum())}/{len(sec)}, "
          f"Wilcoxon p={p_w:.4f}")
    w2, p_w2 = wilcoxon(sec["frac50_hi"], sec["frac50_lo"])
    print(f"[b] POST per-section Wilcoxon frac<=50um (n={len(sec)}): "
          f"median frac50_diff={sec['frac50_diff'].median():.3f}, "
          f"Wilcoxon p={p_w2:.4f}")

# Spearman MM2 vs distance (per-cell pooled post; per-section median)
rho_c, p_c = spearmanr(post["MM2_4gene"], post["dist_T"])
print(f"\n[c] POST per-cell Spearman MM2 vs dist_T (n={len(post)}): "
      f"rho={rho_c:.3f}, p={p_c:.2e}")
sec_med = post.groupby("orig.ident").agg(mm2_med=("MM2_4gene", "median"),
                                         dist_med=("dist_T", "median"))
rho_s, p_s = spearmanr(sec_med["mm2_med"], sec_med["dist_med"])
print(f"[c] POST per-section-median Spearman (n={len(sec_med)}): "
      f"rho={rho_s:.3f}, p={p_s:.4f}")

# tertile sensitivity (post, per-section tertiles)
def tertile_stats(df):
    rows = []
    for sec_, g in df.groupby("orig.ident"):
        if len(g) < 150:
            continue
        q1, q2 = g["MM2_4gene"].quantile([1/3, 2/3])
        hi3 = g[g["MM2_4gene"] > q2]
        lo3 = g[g["MM2_4gene"] <= q1]
        rows.append({"section": sec_,
                     "dist_hi3": hi3["dist_T"].median(),
                     "dist_lo3": lo3["dist_T"].median(),
                     "frac50_hi3": hi3["has_T_50"].mean(),
                     "frac50_lo3": lo3["has_T_50"].mean()})
    return pd.DataFrame(rows)
t3 = tertile_stats(post)
print("\n[tertile sensitivity, POST] top vs bottom tertile:")
print(t3.to_string(index=False))
if len(t3) >= 3:
    w3, p3 = wilcoxon(t3["dist_hi3"], t3["dist_lo3"])
    print(f"  Wilcoxon dist (n={len(t3)}): p={p3:.4f}, "
          f"median diff={(t3['dist_hi3']-t3['dist_lo3']).median():.1f} um")
    w3b, p3b = wilcoxon(t3["frac50_hi3"], t3["frac50_lo3"])
    print(f"  Wilcoxon frac50 (n={len(t3)}): p={p3b:.4f}, "
          f"median diff={(t3['frac50_hi3']-t3['frac50_lo3']).median():.3f}")

# ---------- Secondary: PRE sections ----------
pre = tum[tum["timepoint"] == "pre"].copy()
print(f"\n===== PRE tumor cells: {len(pre)} (secondary) =====")
pre_med = pre.groupby("orig.ident")["MM2_4gene"].transform("median")
pre["MM2_hi"] = pre["MM2_4gene"] > pre_med
hi_p = pre[pre["MM2_hi"]]; lo_p = pre[~pre["MM2_hi"]]
u_p, p_mp = mannwhitneyu(hi_p["dist_T"], lo_p["dist_T"])
print(f"[a] PRE pooled MW: hi (n={len(hi_p)}) median={hi_p['dist_T'].median():.1f} um "
      f"vs lo (n={len(lo_p)}) median={lo_p['dist_T'].median():.1f} um; p={p_mp:.2e}")
print(f"[b] PRE pooled frac<=50um: hi={hi_p['has_T_50'].mean():.3f} vs lo={lo_p['has_T_50'].mean():.3f}")
sec_pre = section_stats(pre, "pre")
print("\nper-section (PRE) table:")
print(sec_pre.to_string(index=False))
sec_pre.to_csv(OUT + "/mm2_tcell_per_section_pre.csv", index=False)
if len(sec_pre) >= 3:
    w_p, p_wp = wilcoxon(sec_pre["dist_hi_med"], sec_pre["dist_lo_med"])
    print(f"\n[a] PRE per-section Wilcoxon (n={len(sec_pre)}): "
          f"median dist_diff={sec_pre['dist_diff'].median():.1f} um, p={p_wp:.4f}")

# ---------- Sensitivity: author's mes2 instead of 4-gene MM2 (POST) ----------
post2 = tum[tum["timepoint"] == "post"].copy()
post2_med = post2.groupby("orig.ident")["mes2"].transform("median")
post2["mes2_hi"] = post2["mes2"] > post2_med
hi2 = post2[post2["mes2_hi"]]; lo2 = post2[~post2["mes2_hi"]]
u_s, p_sens = mannwhitneyu(hi2["dist_T"], lo2["dist_T"])
print(f"\n[Sensitivity, author mes2, POST] hi (n={len(hi2)}) median={hi2['dist_T'].median():.1f} um "
      f"vs lo (n={len(lo2)}) median={lo2['dist_T'].median():.1f} um; MW p={p_sens:.2e}")
print(f"  frac<=50um hi={hi2['has_T_50'].mean():.3f} vs lo={lo2['has_T_50'].mean():.3f}")
print("DONE")
