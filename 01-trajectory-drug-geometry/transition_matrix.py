"""B-side compositional transition analysis: which trajectory bins lose/gain cells?

Point 1 found trajectory-bin displacement is 83% compositional (redistribution
across bins, not within-bin movement). This analysis characterizes the flow:
for each of the 43 treatment-vs-vehicle comparisons, compute the fraction of
tumor cells in each TRAJ bin (quintiles) in vehicle vs treated, and aggregate
the changes. This turns the net 'drugs push toward early' into a flow map.

Uses only obs from adata_tumor_selected_full.h5ad (TRAJ, Treatment, TissueID).
No expression matrix needed.
"""
import anndata as ad
import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

H5 = "/home/hatch/workspace/bioinfo_downloads/gbm-point1/b-side/adata_tumor_selected_full.h5ad"
OUT = "/home/hatch/workspace/bioinfo_downloads/gbm-point1/colab_outputs"
import os
os.makedirs(OUT, exist_ok=True)

a = ad.read_h5ad(H5, backed="r")
obs = a.obs[["TissueID", "Treatment", "TRAJ"]].copy()
print("cells:", len(obs))
# bins: quintiles of TRAJ over ALL tumor cells (matches Point 1 five-bin approach)
qs = np.quantile(obs["TRAJ"], [0, 0.2, 0.4, 0.6, 0.8, 1.0])
qs[0] -= 1e-9
qs[-1] += 1e-9
print("bin edges:", np.round(qs, 3))
obs["bin"] = pd.cut(obs["TRAJ"], bins=qs, labels=[0, 1, 2, 3, 4], include_lowest=True).astype(int)
# bin 0 = most early, bin 4 = most late
print("overall bin fractions:", obs["bin"].value_counts(normalize=True).sort_index().round(3).to_dict())

# per-comparison bin fractions
rows = []
for (tid, drug), g in obs[obs["Treatment"] != "vehicle"].groupby(["TissueID", "Treatment"]):
    veh = obs[(obs["TissueID"] == tid) & (obs["Treatment"] == "vehicle")]
    if len(veh) < 100 or len(g) < 100:
        continue
    fv = veh["bin"].value_counts(normalize=True).sort_index()
    ft = g["bin"].value_counts(normalize=True).sort_index()
    for b in range(5):
        rows.append({"tissue": tid, "drug": drug, "bin": b,
                     "frac_vehicle": fv.get(b, 0.0), "frac_treated": ft.get(b, 0.0),
                     "n_vehicle": len(veh), "n_treated": len(g)})
    # net early-ward flow: fraction in bins 0-1 (early) vs bins 3-4 (late)
    rows.append({"tissue": tid, "drug": drug, "bin": 99,
                 "frac_vehicle": fv.get(0, 0) + fv.get(1, 0),
                 "frac_treated": ft.get(0, 0) + ft.get(1, 0),
                 "n_vehicle": len(veh), "n_treated": len(g)})
    rows.append({"tissue": tid, "drug": drug, "bin": 98,
                 "frac_vehicle": fv.get(3, 0) + fv.get(4, 0),
                 "frac_treated": ft.get(3, 0) + ft.get(4, 0),
                 "n_vehicle": len(veh), "n_treated": len(g)})

tab = pd.DataFrame(rows)
tab["delta"] = tab["frac_treated"] - tab["frac_vehicle"]
tab.to_csv(os.path.join(OUT, "transition_matrix.csv"), index=False)
print(f"\ncomparisons: {tab['tissue'].nunique()} tissues, {len(tab[tab['bin']<5]['drug'].unique())} drugs, "
      f"{tab[tab['bin']<5].groupby(['tissue','drug']).ngroups} pairs")

print("\nMean change in bin fraction (treated - vehicle), across comparisons:")
for b in range(5):
    d = tab[tab["bin"] == b]["delta"]
    w, p = wilcoxon(d)
    print(f"  bin {b} ({'early' if b==0 else 'late' if b==4 else 'mid'}): "
          f"mean {d.mean():+.4f}, median {d.median():+.4f}, "
          f"up {(d>0).sum()}/{(d<0).sum()} down, Wilcoxon p={p:.4f} (n={len(d)})")
for b, label in ((99, "early bins 0-1"), (98, "late bins 3-4")):
    d = tab[tab["bin"] == b]["delta"]
    w, p = wilcoxon(d)
    print(f"  {label}: mean {d.mean():+.4f}, median {d.median():+.4f}, "
          f"up {(d>0).sum()}/{(d<0).sum()} down, Wilcoxon p={p:.4f} (n={len(d)})")

print("\nPer-drug mean delta by bin (rows=drug, cols=bin 0..4):")
piv = tab[tab["bin"] < 5].pivot_table(index="drug", columns="bin", values="delta", aggfunc="mean")
print(piv.round(3).to_string())
piv.to_csv(os.path.join(OUT, "transition_per_drug.csv"))

# net flow: total mass moving from late to early (sum of positive deltas in early bins)
print("\nNet mass gained by early bins 0-1 (mean across comparisons):",
      round(float(tab[tab['bin']==99]['delta'].mean()), 4))
print("Net mass lost by late bins 3-4 (mean across comparisons):",
      round(float(tab[tab['bin']==98]['delta'].mean()), 4))
print("DONE")
