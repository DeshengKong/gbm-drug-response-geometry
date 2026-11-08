#!/usr/bin/env python3
"""Fix-up for the Meylan analysis: correct the paired-test bug.
Bug: orig.ident is 'P05_post'/'P05_pre', so pivoting on it misaligns
pre/post (each 'patient' appears in only one column -> NaN Wilcoxon).
Fix: derive patient_base by stripping the _pre/_post suffix.
Also adds a more sensitive spatial metric: distance from each tumor cell
to the nearest T cell, and from each T cell to the nearest tumor cell.
Saves corrected tables + a summary JSON.
"""
import json
import numpy as np
import pandas as pd
from scipy import stats, spatial

rng = np.random.default_rng(7)

def base(p):
    return p.rsplit('_', 1)[0]

# ---- corrected paired tests on saved tables ----
comp = pd.read_csv('meylan_analysis/q1_composition.csv')
comp['patient_base'] = comp['patient'].apply(base)
nb = pd.read_csv('meylan_analysis/q1_q4_neighborhood.csv')
nb['patient_base'] = nb['patient'].apply(base)

out = {}
for df, cols in [(comp, ['T_frac', 'myeloid_frac', 'tumor_frac', 'T_per_tumor']),
                 (nb, ['T_frac_near', 'myeloid_frac_near'])]:
    for col in cols:
        piv = df.pivot(index='patient_base', columns='timepoint', values=col)
        d = piv['post'] - piv['pre']
        w = stats.wilcoxon(piv['post'], piv['pre'])
        out[col] = dict(n=int(len(piv)), median_diff=float(d.median()),
                        mean_diff=float(d.mean()),
                        n_up=int((d > 0).sum()), n_down=int((d < 0).sum()),
                        wilcoxon_stat=float(w.statistic), wilcoxon_p=float(w.pvalue),
                        per_patient={k: float(v) for k, v in d.items()})
        print(f"{col}: n={len(piv)} median_diff={d.median():.4f} "
              f"up={int((d>0).sum())}/down={int((d<0).sum())} "
              f"Wilcoxon p={w.pvalue:.4g}")

# ---- N1 permutation null (corrected): permute pre/post within patient ----
piv = comp.pivot(index='patient_base', columns='timepoint', values='T_frac')
obs = (piv['post'] - piv['pre']).median()
nulls = []
for _ in range(10000):
    flip = rng.random(len(piv)) < 0.5
    d = np.where(flip, piv['pre'].values - piv['post'].values,
                 piv['post'].values - piv['pre'].values)
    nulls.append(np.median(d))
nulls = np.array(nulls)
p_emp = float((np.abs(nulls) >= abs(obs)).mean())
out['N1_T_frac_perm'] = dict(obs=float(obs),
    null_lo=float(np.percentile(nulls, 2.5)),
    null_hi=float(np.percentile(nulls, 97.5)), p_emp=p_emp)
print(f"N1 T_frac perm: obs={obs:.4f} null95%=[{np.percentile(nulls,2.5):.4f},"
      f"{np.percentile(nulls,97.5):.4f}] p_emp={p_emp:.4f}")

# ---- sensitive spatial metric: T-cell <-> tumor distances ----
meta = pd.read_csv('meylan_meta.csv.gz', index_col=0,
                   usecols=['orig.ident', 'timepoint', 'tier1', 'tier1_T', 'x', 'y'])
T_STATES = ['CD8 Early activated', 'CD8 Stem-like', 'CD8 Tissue resident',
            'CD8 IFN-responsive', 'CD8 Proliferative', 'CD8 Chemokine producing',
            'CD8 NK-like', 'CD8 Undefined', 'CD4 Tissue resident', 'CD4 T regulatory']
meta['is_T'] = meta['tier1_T'].isin(T_STATES)
meta['is_tumor'] = meta['tier1'] == 'Tumor cells'
rows = []
for fov, g in meta.groupby('orig.ident'):
    xy = g[['x', 'y']].values
    txy = g.loc[g['is_T'], ['x', 'y']].values
    tumxy = g.loc[g['is_tumor'], ['x', 'y']].values
    # tumor -> nearest T cell
    treeT = spatial.cKDTree(txy)
    d_tum_to_T, _ = treeT.query(tumxy, k=1)
    # T cell -> nearest tumor cell
    treeTum = spatial.cKDTree(tumxy)
    d_T_to_tum, _ = treeTum.query(txy, k=1)
    rows.append(dict(patient=fov, timepoint=g['timepoint'].iloc[0],
                     n_tumor=len(tumxy), n_T=len(txy),
                     med_tumor_to_T=float(np.median(d_tum_to_T)),
                     frac_tumor_with_T_50=float((d_tum_to_T <= 50).mean()),
                     frac_tumor_with_T_20=float((d_tum_to_T <= 20).mean()),
                     med_T_to_tumor=float(np.median(d_T_to_tum)),
                     frac_T_near_tumor_20=float((d_T_to_tum <= 20).mean())))
    print(fov, 'med tumor->T %.1f um, T within 20um of tumor %.3f'
          % (np.median(d_tum_to_T), (d_T_to_tum <= 20).mean()))
dist = pd.DataFrame(rows)
dist['patient_base'] = dist['patient'].apply(base)
dist.to_csv('meylan_analysis/q1_distances.csv', index=False)
for col in ['med_tumor_to_T', 'frac_tumor_with_T_50', 'frac_tumor_with_T_20',
            'med_T_to_tumor', 'frac_T_near_tumor_20']:
    piv = dist.pivot(index='patient_base', columns='timepoint', values=col)
    d = piv['post'] - piv['pre']
    w = stats.wilcoxon(piv['post'], piv['pre'])
    out[col] = dict(n=int(len(piv)), median_diff=float(d.median()),
                    n_up=int((d > 0).sum()), n_down=int((d < 0).sum()),
                    wilcoxon_p=float(w.pvalue),
                    per_patient={k: float(v) for k, v in d.items()})
    print(f"{col}: median_diff={d.median():.4f} up={int((d>0).sum())}/"
          f"down={int((d<0).sum())} Wilcoxon p={w.pvalue:.4g}")

with open('meylan_analysis/summary.json', 'w') as f:
    json.dump(out, f, indent=2)
print('SAVED summary.json + q1_distances.csv')
