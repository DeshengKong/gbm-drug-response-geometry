#!/usr/bin/env python3
"""Independent verification + delta table for the #37 correction (_r3).

Recomputes EVERY Wang headline statistic from the per-sample CSVs using the
identical methods as wang_analysis.py / wang_ap1.py (scipy wilcoxon defaults,
10k within-pair label permutation with np.random.default_rng(0), Spearman for
ligands), for OLD (backup, pre-correction) vs NEW (_r3) files, both #39
variants. Prints a delta table. Also verifies the #37 rows carry the corrected
labels and that all other rows reproduce.
"""
import numpy as np
import pandas as pd
from scipy import stats

OLD = "/home/hatch/workspace/bioinfo_downloads/wang2022/wang_analysis/backup_pre_r3_20261109"
NEW = "/home/hatch/workspace/bioinfo_downloads/wang2022/wang_analysis"
REPO = "/home/hatch/workspace/gbm-drug-geometry-code/04-gbmap-wang/wang_analysis"

MM_ANNOT = {0: "Cycling", 1: "NPC-neuronal-like 1", 2: "Gliosis-hypoxia",
            3: "AC-gliosis-like 2", 4: "NPC-OPC-like 2", 5: "AC-gliosis-like 1",
            6: "NPC-neuronal-like 2", 7: "NPC-OPC-like 1", 8: "AC-NPC-OPC-like 1",
            9: "AC-progenitor-like", 10: "NPC-neuronal-like 3", 11: "OPC-like"}


def main_stats(tab):
    out = {}
    piv = tab.pivot(index="pair", columns="timepoint", values="traj_median_tumor").dropna()
    d = piv["Recurrent"] - piv["Primary"]
    out["n_pairs"] = len(piv)
    out["traj_med_diff"] = float(d.median())
    out["up"] = int((d > 0).sum()); out["down"] = int((d < 0).sum())
    out["wilcoxon_p"] = float(stats.wilcoxon(piv["Recurrent"], piv["Primary"]).pvalue)
    pr, rc = piv["Primary"].values, piv["Recurrent"].values
    rng = np.random.default_rng(0)
    obs = float(np.median(rc - pr))
    nulls = []
    for _ in range(10000):
        flip = rng.random(len(pr)) < 0.5
        nulls.append(np.median(np.where(flip, pr - rc, rc - pr)))
    nulls = np.array(nulls)
    out["perm_p"] = float((np.abs(nulls) >= abs(obs)).mean())
    out["perm_ci"] = (float(np.percentile(nulls, 2.5)), float(np.percentile(nulls, 97.5)))
    pivm = tab.pivot(index="pair", columns="timepoint", values="traj_mean_tumor").dropna()
    dm = pivm["Recurrent"] - pivm["Primary"]
    out["mean_med_diff"] = float(dm.median())
    out["mean_wilcoxon_p"] = float(stats.wilcoxon(pivm["Recurrent"], pivm["Primary"]).pvalue)
    out["modules"] = {}
    for i in range(12):
        pm = tab.pivot(index="pair", columns="timepoint", values=f"MM{i}_tumor").dropna()
        dd = pm["Recurrent"] - pm["Primary"]
        out["modules"][i] = (float(dd.median()),
                             float(stats.wilcoxon(pm["Recurrent"], pm["Primary"]).pvalue))
    pt = tab.pivot(index="pair", columns="timepoint", values="frac_t").dropna()
    dt = pt["Recurrent"] - pt["Primary"]
    out["t_med_diff"] = float(dt.median())
    out["t_up"] = int((dt > 0).sum()); out["t_down"] = int((dt < 0).sum())
    out["t_p"] = float(stats.wilcoxon(pt["Recurrent"], pt["Primary"]).pvalue)
    out["ligands"] = {}
    for lig in ["VEGFA", "VCAN", "ANXA1"]:
        sub = tab.dropna(subset=[f"lig_{lig}", "traj_median_tumor"])
        r, p = stats.spearmanr(sub[f"lig_{lig}"], sub["traj_median_tumor"])
        out["ligands"][lig] = (len(sub), float(r), float(p))
    out["pair37_diff"] = float(d.loc["#37"]) if "#37" in d.index else None
    return out


def ap1_stats(tab):
    out = {}
    piv = tab.pivot(index="pair", columns="tp")
    for col in ["AP1_core", "AP1_ext", "AP1_core_all", "FOS_tumor", "JUN_tumor",
                "FOSL1_tumor", "ATF3_tumor", "traj_med_tumor"]:
        d = (piv[(col, "Recurrent")] - piv[(col, "Primary")]).dropna()
        w, p = stats.wilcoxon(d)
        out[col] = (float(d.median()), int((d > 0).sum()), int((d < 0).sum()),
                    float(p), len(d))
    r, p = stats.spearmanr(tab[tab["tp"] == "Primary"]["AP1_core"],
                           tab[tab["tp"] == "Primary"]["traj_med_tumor"])
    out["AP1_vs_TRAJ_primary"] = (float(r), float(p))
    out["pair37_ap1_diff"] = float((piv[("AP1_core", "Recurrent")] -
                                    piv[("AP1_core", "Primary")]).loc["#37"])
    return out


def check_labels(tab, tag, tp_col="timepoint"):
    r32 = tab[tab["gsm"] == "GSM5319532"].iloc[0]
    r59 = tab[tab["gsm"] == "GSM5319559"].iloc[0]
    ok = (r32["pair"] == "#37" and r32[tp_col] == "Recurrent"
          and r59["pair"] == "#37" and r59[tp_col] == "Primary")
    print(f"[{tag}] #37 labels: GSM5319532={r32[tp_col]} "
          f"GSM5319559={r59[tp_col]} -> {'OK' if ok else 'STALE!'}")
    return ok


def compare_rows(old, new, tag, tp_col="timepoint"):
    m = old.merge(new, on="gsm", suffixes=("_old", "_new"))
    num = [c for c in old.columns if old[c].dtype.kind in "fi" and c != "gsm"]
    diffs = {}
    for c in num:
        d = (m[c + "_old"] - m[c + "_new"]).abs()
        diffs[c] = float(d.max())
    changed = {c: v for c, v in diffs.items() if v > 0}
    lab = (m["pair_old"] != m["pair_new"]) | (m[tp_col + "_old"] != m[tp_col + "_new"])
    print(f"[{tag}] rows old={len(old)} new={len(new)}; numeric cols with any "
          f"change: {changed if changed else 'NONE'}; label-changed rows: "
          f"{m[lab][['gsm', 'pair_old', tp_col + '_old', 'pair_new', tp_col + '_new']].to_dict('records') if lab.any() else 'NONE'}")
    return changed


def show(name, old, new, fmt="{:.4f}"):
    print(f"  {name:38s} old={fmt.format(old) if isinstance(old, float) else old}  "
          f"new={fmt.format(new) if isinstance(new, float) else new}")


for variant, old_main, new_main, old_ap1, new_ap1 in [
        ("BASE KEEP (#39=GSM5319574)",
         f"{OLD}/wang_per_sample.csv", f"{NEW}/wang_per_sample_r3.csv",
         f"{REPO}/wang_ap1.csv", f"{NEW}/wang_ap1_r3.csv"),
        ("SWAP (#39=GSM5319563)",
         f"{OLD}/wang_per_sample_s39.csv", f"{NEW}/wang_per_sample_s39_r3.csv",
         f"{REPO}/wang_ap1_s39.csv", f"{NEW}/wang_ap1_s39_r3.csv")]:
    print(f"\n===== {variant} =====")
    ot = pd.read_csv(old_main); nt = pd.read_csv(new_main)
    oa = pd.read_csv(old_ap1); na = pd.read_csv(new_ap1)
    check_labels(nt, variant + " main")
    check_labels(na, variant + " ap1", tp_col="tp")
    compare_rows(ot, nt, variant + " main")
    compare_rows(oa, na, variant + " ap1", tp_col="tp")
    o, n = main_stats(ot), main_stats(nt)
    print("  -- main TRAJ --")
    show("n pairs", o["n_pairs"], n["n_pairs"], fmt="{}")
    show("TRAJ median diff", o["traj_med_diff"], n["traj_med_diff"])
    show("up / down", f"{o['up']}/{o['down']}", f"{n['up']}/{n['down']}", fmt="{}")
    show("Wilcoxon p", o["wilcoxon_p"], n["wilcoxon_p"], fmt="{:.4g}")
    show("permutation p", o["perm_p"], n["perm_p"], fmt="{:.4f}")
    show("perm null 95% CI", o["perm_ci"], n["perm_ci"], fmt="{}")
    show("#37 pair diff", o["pair37_diff"], n["pair37_diff"])
    show("mean TRAJ median diff", o["mean_med_diff"], n["mean_med_diff"])
    show("mean Wilcoxon p", o["mean_wilcoxon_p"], n["mean_wilcoxon_p"], fmt="{:.4g}")
    print("  -- per-module (median diff / Wilcoxon p; Bonferroni 0.05/12=0.00417) --")
    for i in range(12):
        od, op = o["modules"][i]; nd, np_ = n["modules"][i]
        flag_o = "BONF" if op < 0.05 / 12 else ("nom" if op < 0.05 else "ns")
        flag_n = "BONF" if np_ < 0.05 / 12 else ("nom" if np_ < 0.05 else "ns")
        print(f"  MM{i:2d} {MM_ANNOT[i]:22s} old={od:+.4f} / {op:.4g} [{flag_o}]  "
              f"new={nd:+.4f} / {np_:.4g} [{flag_n}]")
    print("  -- T-cell --")
    show("T median diff", o["t_med_diff"], n["t_med_diff"])
    show("T up / down", f"{o['t_up']}/{o['t_down']}", f"{n['t_up']}/{n['t_down']}", fmt="{}")
    show("T Wilcoxon p", o["t_p"], n["t_p"], fmt="{:.4g}")
    print("  -- ligands (Spearman vs tumor TRAJ) --")
    for lig in ["VEGFA", "VCAN", "ANXA1"]:
        on, orr, op = o["ligands"][lig]; nn, nrr, npp = n["ligands"][lig]
        print(f"  {lig:6s} old n={on} rho={orr:.3f} p={op:.4g}  new n={nn} rho={nrr:.3f} p={npp:.4g}")
    oa_s, na_s = ap1_stats(oa), ap1_stats(na)
    print("  -- AP-1 (median diff / up-down / Wilcoxon p) --")
    for col in ["AP1_core", "AP1_ext", "AP1_core_all", "FOS_tumor", "JUN_tumor",
                "FOSL1_tumor", "ATF3_tumor", "traj_med_tumor"]:
        od, ou, odwn, op, on = oa_s[col]; nd, nu, ndwn, npp, nn = na_s[col]
        print(f"  {col:14s} old={od:+.4f} {ou}/{odwn} p={op:.4g} (n={on})  "
              f"new={nd:+.4f} {nu}/{ndwn} p={npp:.4g} (n={nn})")
    print(f"  Primary AP1 vs TRAJ: old rho={oa_s['AP1_vs_TRAJ_primary'][0]:.3f} "
          f"p={oa_s['AP1_vs_TRAJ_primary'][1]:.4g}  new rho={na_s['AP1_vs_TRAJ_primary'][0]:.3f} "
          f"p={na_s['AP1_vs_TRAJ_primary'][1]:.4g}")
    print(f"  #37 AP-1 core diff: old={oa_s['pair37_ap1_diff']:+.4f}  "
          f"new={na_s['pair37_ap1_diff']:+.4f}")
print("\nDONE")
