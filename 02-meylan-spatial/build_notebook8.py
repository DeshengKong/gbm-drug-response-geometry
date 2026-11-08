#!/usr/bin/env python3
"""Build the Meylan Q1/Q2/Q4 analysis notebook (Python, runs in the
gbm-meylan Colab session on the extracted tables).
Q1: T-cell infiltration post-oHSV (composition + tumor-neighborhood).
Q2: ligand zoning (VEGFA/VCAN/ANXA1 vs tier2 state + local myeloid density).
Q4: myeloid spatial rearrangement.
Q3 (TCR) is dead — no usable clonotype mapping in the RDS or GEO.
"""
import nbformat as nbf

NB_PATH = "/home/hatch/workspace/bioinfo_downloads/gse296577/gbm-meylan-analysis.ipynb"

cells = []

cells.append(nbf.v4.new_markdown_cell(
    "# Meylan GSE296577 — Q1 / Q2 / Q4 spatial analysis\n"
    "Data: extracted tables (meta 23 cols, 22-gene counts, x/y in meta).\n"
    "Q3 (TCR) is dead: no usable clonotype mapping in the RDS or GEO.\n"
    "Patient is the unit for paired tests (n=8). Spatial nulls are\n"
    "permutation-based. No trajectory projection (panel coverage 4.8%)."))

cells.append(nbf.v4.new_code_cell(
    "import os, json\n"
    "import numpy as np, pandas as pd\n"
    "from scipy import stats, spatial\n"
    "import matplotlib\n"
    "matplotlib.use('Agg')\n"
    "import matplotlib.pyplot as plt\n"
    "os.makedirs('meylan_analysis', exist_ok=True)\n"
    "rng = np.random.default_rng(7)\n"
    "print('ok')"))

cells.append(nbf.v4.new_code_cell(
    "# Load meta + gene counts (row order matches: verified identical barcodes)\n"
    "meta = pd.read_csv('meylan_meta.csv.gz', index_col=0)\n"
    "genes = pd.read_csv('meylan_genes.csv.gz', index_col=0)\n"
    "assert (meta.index == genes.index).all(), 'row order mismatch'\n"
    "print('meta:', meta.shape, '| genes:', genes.shape)\n"
    "print('patients:', sorted(meta['orig.ident'].unique()))\n"
    "print('timepoints:', meta['timepoint'].value_counts().to_dict())\n"
    "T_STATES = ['CD8 Early activated','CD8 Stem-like','CD8 Tissue resident',\n"
    "            'CD8 IFN-responsive','CD8 Proliferative','CD8 Chemokine producing',\n"
    "            'CD8 NK-like','CD8 Undefined','CD4 Tissue resident','CD4 T regulatory']\n"
    "meta['is_T'] = meta['tier1_T'].isin(T_STATES)\n"
    "meta['is_myeloid'] = meta['tier1'] == 'Macrophages/Monocytes'\n"
    "meta['is_tumor'] = meta['tier1'] == 'Tumor cells'\n"
    "meta.to_pickle('meylan_analysis/meta_flags.pkl')\n"
    "print('T cells:', int(meta['is_T'].sum()), '| myeloid:', int(meta['is_myeloid'].sum()),\n"
    "      '| tumor:', int(meta['is_tumor'].sum()))"))

cells.append(nbf.v4.new_code_cell(
    "# Q1 (part 1): composition per patient x timepoint\n"
    "rows = []\n"
    "for (pat, tp), g in meta.groupby(['orig.ident', 'timepoint']):\n"
    "    rows.append(dict(patient=pat, timepoint=tp, n=len(g),\n"
    "        T_frac=g['is_T'].mean(), myeloid_frac=g['is_myeloid'].mean(),\n"
    "        tumor_frac=g['is_tumor'].mean(),\n"
    "        T_per_tumor=g['is_T'].sum() / max(g['is_tumor'].sum(), 1)))\n"
    "comp = pd.DataFrame(rows)\n"
    "comp.to_csv('meylan_analysis/q1_composition.csv', index=False)\n"
    "print(comp.to_string())\n"
    "piv = comp.pivot(index='patient', columns='timepoint', values='T_frac')\n"
    "print('\\npaired Wilcoxon T_frac post-pre:')\n"
    "print(stats.wilcoxon(piv['post'], piv['pre']))\n"
    "print('median post-pre:', float((piv['post'] - piv['pre']).median()))"))

cells.append(nbf.v4.new_code_cell(
    "# Q1 (part 2) + Q4: distance to nearest tumor cell, per FOV (cKDTree)\n"
    "# Tumor-neighborhood = within 50 um of any tumor cell.\n"
    "meta['dist_tumor'] = np.nan\n"
    "meta['near_tumor'] = False\n"
    "for fov, g in meta.groupby('orig.ident'):\n"
    "    xy = g[['x', 'y']].values\n"
    "    txy = g.loc[g['is_tumor'], ['x', 'y']].values\n"
    "    tree = spatial.cKDTree(txy)\n"
    "    d, _ = tree.query(xy, k=1)\n"
    "    meta.loc[g.index, 'dist_tumor'] = d\n"
    "    meta.loc[g.index, 'near_tumor'] = d <= 50.0\n"
    "    print(fov, 'median dist to tumor: %.1f um, near<=50um: %.3f' % (\n"
    "        np.median(d), (d <= 50).mean()))\n"
    "meta[['dist_tumor', 'near_tumor']].to_pickle('meylan_analysis/dist_tumor.pkl')\n"
    "print('done')"))

cells.append(nbf.v4.new_code_cell(
    "# Q1 (part 3): T and myeloid fractions inside tumor neighborhoods\n"
    "rows = []\n"
    "for (pat, tp), g in meta.groupby(['orig.ident', 'timepoint']):\n"
    "    near = g[g['near_tumor']]\n"
    "    rows.append(dict(patient=pat, timepoint=tp, n_near=len(near),\n"
    "        T_frac_near=near['is_T'].mean(), myeloid_frac_near=near['is_myeloid'].mean(),\n"
    "        T_frac_all=g['is_T'].mean(), myeloid_frac_all=g['is_myeloid'].mean()))\n"
    "nb = pd.DataFrame(rows)\n"
    "nb.to_csv('meylan_analysis/q1_q4_neighborhood.csv', index=False)\n"
    "print(nb.to_string())\n"
    "for col in ['T_frac_near', 'myeloid_frac_near']:\n"
    "    piv = nb.pivot(index='patient', columns='timepoint', values=col)\n"
    "    print(col, 'paired Wilcoxon:', stats.wilcoxon(piv['post'], piv['pre']))\n"
    "    print('  median post-pre:', float((piv['post'] - piv['pre']).median()))"))

cells.append(nbf.v4.new_code_cell(
    "# Q2A: ligand positivity by tier2 tumor state\n"
    "lig = ['VEGFA', 'VCAN', 'ANXA1']\n"
    "tumor = meta[meta['is_tumor']].copy()\n"
    "for L in lig:\n"
    "    tumor[L + '_pos'] = genes.loc[tumor.index, L] > 0\n"
    "tumor['lig_any'] = tumor[[L + '_pos' for L in lig]].any(axis=1)\n"
    "tab = tumor.groupby('tier2')[[L + '_pos' for L in lig] + ['lig_any']].mean()\n"
    "cnt = tumor.groupby('tier2').size()\n"
    "print('ligand positivity by tier2 state (tumor cells only):')\n"
    "print(tab.round(4).to_string())\n"
    "print('counts:'); print(cnt.to_string())\n"
    "tab.to_csv('meylan_analysis/q2a_ligand_by_tier2.csv')\n"
    "# per-patient per-state rates for the hypoxic MES-like 2 vs rest test\n"
    "rows = []\n"
    "for (pat, tp, st), g in tumor.groupby(['orig.ident', 'timepoint', 'tier2']):\n"
    "    rows.append(dict(patient=pat, timepoint=tp, state=st, n=len(g),\n"
    "        VEGFA=g['VEGFA_pos'].mean(), VCAN=g['VCAN_pos'].mean(),\n"
    "        ANXA1=g['ANXA1_pos'].mean(), lig_any=g['lig_any'].mean(),\n"
    "        mes2=g['mes2'].mean()))\n"
    "q2a = pd.DataFrame(rows)\n"
    "q2a.to_csv('meylan_analysis/q2a_patient_state.csv', index=False)\n"
    "hyp = q2a[q2a['state'] == 'Tumor MES-like 2 (hypoxia)']\n"
    "oth = q2a[q2a['state'] != 'Tumor MES-like 2 (hypoxia)']\n"
    "for L in ['VEGFA', 'VCAN', 'ANXA1', 'lig_any']:\n"
    "    print(L, 'hypoxic MES2 mean %.4f vs other states %.4f (Mann-Whitney U p=%.3g)' % (\n"
    "        hyp[L].mean(), oth[L].mean(),\n"
    "        stats.mannwhitneyu(hyp[L], oth[L]).pvalue))"))

cells.append(nbf.v4.new_code_cell(
    "# Q2B: local myeloid density around ligand+ vs ligand- tumor cells\n"
    "# Stratified subsample: up to 8k tumor cells per patient x timepoint.\n"
    "R = 50.0\n"
    "tumor = meta[meta['is_tumor']].copy()\n"
    "for L in ['VEGFA', 'VCAN', 'ANXA1']:\n"
    "    tumor[L + '_pos'] = genes.loc[tumor.index, L] > 0\n"
    "tumor['lig_any'] = tumor[['VEGFA_pos', 'VCAN_pos', 'ANXA1_pos']].any(axis=1)\n"
    "sub = []\n"
    "for (pat, tp), g in tumor.groupby(['orig.ident', 'timepoint']):\n"
    "    take = min(8000, len(g))\n"
    "    sub.append(g.sample(n=take, random_state=7))\n"
    "sub = pd.concat(sub)\n"
    "print('subsample tumor cells:', len(sub))\n"
    "res = []\n"
    "for fov, g in sub.groupby('orig.ident'):\n"
    "    allf = meta[meta['orig.ident'] == fov]\n"
    "    tree = spatial.cKDTree(allf[['x', 'y']].values)\n"
    "    my = allf['is_myeloid'].values\n"
    "    idx = tree.query_ball_point(g[['x', 'y']].values, r=R)\n"
    "    dens = np.array([my[list(ii)].mean() if len(ii) else np.nan for ii in idx])\n"
    "    cnt = np.array([len(ii) for ii in idx])\n"
    "    tmp = g[['tier2', 'lig_any', 'VEGFA_pos', 'VCAN_pos', 'ANXA1_pos']].copy()\n"
    "    tmp['myeloid_dens_50'] = dens\n"
    "    tmp['n_within_50'] = cnt\n"
    "    tmp['patient'] = fov\n"
    "    res.append(tmp)\n"
    "q2b = pd.concat(res)\n"
    "q2b.to_csv('meylan_analysis/q2b_local_myeloid.csv.gz', compression='gzip')\n"
    "print(q2b.groupby('lig_any')['myeloid_dens_50'].agg(['mean', 'median', 'count']))\n"
    "print('overall: ligand+ vs ligand- myeloid density (Mann-Whitney):',\n"
    "      stats.mannwhitneyu(q2b.loc[q2b['lig_any'], 'myeloid_dens_50'],\n"
    "                          q2b.loc[~q2b['lig_any'], 'myeloid_dens_50']))\n"
    "# paired per patient-timepoint: mean density ligand+ minus ligand-\n"
    "pr = q2b.groupby(['patient', 'lig_any'])['myeloid_dens_50'].mean().unstack()\n"
    "pr['diff'] = pr[True] - pr[False]\n"
    "print(pr.round(4).to_string())\n"
    "print('paired diff median: %.4f, Wilcoxon p=%.3g' % (\n"
    "    pr['diff'].median(), stats.wilcoxon(pr['diff']).pvalue))\n"
    "print('n units with diff>0: %d / %d' % ((pr['diff'] > 0).sum(), len(pr)))"))

cells.append(nbf.v4.new_code_cell(
    "# Q2C: ligand positivity vs mes2 (hypoxia) score — continuous\n"
    "tumor = meta[meta['is_tumor']].copy()\n"
    "for L in ['VEGFA', 'VCAN', 'ANXA1']:\n"
    "    tumor[L + '_pos'] = genes.loc[tumor.index, L] > 0\n"
    "# per patient-timepoint: mean mes2 in ligand+ vs ligand- tumor cells\n"
    "rows = []\n"
    "for (pat, tp), g in tumor.groupby(['orig.ident', 'timepoint']):\n"
    "    for L in ['VEGFA', 'VCAN', 'ANXA1']:\n"
    "        pos = g[g[L + '_pos']]['mes2'].mean()\n"
    "        neg = g[~g[L + '_pos']]['mes2'].mean()\n"
    "        rows.append(dict(patient=pat, timepoint=tp, ligand=L,\n"
    "                         mes2_pos=pos, mes2_neg=neg, diff=pos - neg))\n"
    "q2c = pd.DataFrame(rows)\n"
    "q2c.to_csv('meylan_analysis/q2c_ligand_vs_mes2.csv', index=False)\n"
    "for L in ['VEGFA', 'VCAN', 'ANXA1']:\n"
    "    s = q2c[q2c['ligand'] == L]\n"
    "    print(L, 'mes2 ligand+ minus ligand-: median %.4f, Wilcoxon p=%.3g, %d/%d positive' % (\n"
    "        s['diff'].median(), stats.wilcoxon(s['diff']).pvalue,\n"
    "        (s['diff'] > 0).sum(), len(s)))\n"
    "print(q2c.round(4).to_string())"))

cells.append(nbf.v4.new_code_cell(
    "# Permutation nulls (pre-registered):\n"
    "#  N1: Q1 T_frac — permute pre/post labels within patient, recompute\n"
    "#      median paired difference; 10k permutations.\n"
    "#  N2: Q2B — permute lig_any labels within patient-timepoint, recompute\n"
    "#      the paired median density difference; 1k permutations (subsample).\n"
    "comp = pd.read_csv('meylan_analysis/q1_composition.csv')\n"
    "obs = (comp.pivot(index='patient', columns='timepoint', values='T_frac')['post']\n"
    "       - comp.pivot(index='patient', columns='timepoint', values='T_frac')['pre']).median()\n"
    "nulls = []\n"
    "piv = comp.pivot(index='patient', columns='timepoint', values='T_frac')\n"
    "for _ in range(10000):\n"
    "    flip = rng.random(len(piv)) < 0.5\n"
    "    d = np.where(flip, piv['pre'] - piv['post'], piv['post'] - piv['pre'])\n"
    "    nulls.append(np.median(d))\n"
    "nulls = np.array(nulls)\n"
    "print('N1 Q1 T_frac: observed median diff %.4f, null 95%% CI [%.4f, %.4f], '\n"
    "      'empirical p=%.4f' % (obs, np.percentile(nulls, 2.5), np.percentile(nulls, 97.5),\n"
    "      (np.abs(nulls) >= abs(obs)).mean()))\n"
    "q2b = pd.read_csv('meylan_analysis/q2b_local_myeloid.csv.gz')\n"
    "pr = q2b.groupby(['patient', 'lig_any'])['myeloid_dens_50'].mean().unstack()\n"
    "obs2 = (pr[True] - pr[False]).median()\n"
    "nulls2 = []\n"
    "for _ in range(1000):\n"
    "    q = q2b.copy()\n"
    "    q['lig_any'] = q.groupby('patient')['lig_any'].transform(\n"
    "        lambda s: s.sample(frac=1, random_state=int(rng.integers(1e9))).values)\n"
    "    p2 = q.groupby(['patient', 'lig_any'])['myeloid_dens_50'].mean().unstack()\n"
    "    nulls2.append((p2[True] - p2[False]).median())\n"
    "nulls2 = np.array(nulls2)\n"
    "print('N2 Q2B myeloid density: observed median diff %.4f, null 95%% CI [%.4f, %.4f], '\n"
    "      'empirical p=%.4f' % (obs2, np.percentile(nulls2, 2.5), np.percentile(nulls2, 97.5),\n"
    "      (np.abs(nulls2) >= abs(obs2)).mean()))"))

cells.append(nbf.v4.new_code_cell(
    "# Figures\n"
    "comp = pd.read_csv('meylan_analysis/q1_composition.csv')\n"
    "nb = pd.read_csv('meylan_analysis/q1_q4_neighborhood.csv')\n"
    "q2a = pd.read_csv('meylan_analysis/q2a_ligand_by_tier2.csv', index_col=0)\n"
    "fig, axes = plt.subplots(1, 3, figsize=(13, 4))\n"
    "for ax, col, title in [(axes[0], 'T_frac', 'T-cell fraction'),\n"
    "                      (axes[1], 'myeloid_frac', 'Myeloid fraction'),\n"
    "                      (axes[2], 'tumor_frac', 'Tumor fraction')]:\n"
    "    piv = comp.pivot(index='patient', columns='timepoint', values=col)\n"
    "    for pat in piv.index:\n"
    "        ax.plot([0, 1], [piv.loc[pat, 'pre'], piv.loc[pat, 'post']],\n"
    "                marker='o', color='gray', alpha=0.6, lw=1)\n"
    "    ax.set_xticks([0, 1]); ax.set_xticklabels(['pre', 'post'])\n"
    "    ax.set_title(title); ax.set_xlim(-0.2, 1.2)\n"
    "fig.tight_layout(); fig.savefig('meylan_analysis/fig_q1_composition.png', dpi=150)\n"
    "print('fig_q1 saved')\n"
    "fig, ax = plt.subplots(figsize=(8, 4.5))\n"
    "q2a[['VEGFA_pos', 'VCAN_pos', 'ANXA1_pos']].plot(kind='bar', ax=ax)\n"
    "ax.set_ylabel('fraction of tumor cells expressing')\n"
    "ax.set_title('Ligand positivity by tumor state')\n"
    "plt.xticks(rotation=30, ha='right'); fig.tight_layout()\n"
    "fig.savefig('meylan_analysis/fig_q2a_ligand.png', dpi=150)\n"
    "print('fig_q2a saved')\n"
    "print('ALL DONE')"))

nb = nbf.v4.new_notebook(cells=cells, metadata={
    "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
    "language_info": {"name": "python", "version": "3.11.0"},
})
nbf.write(nb, NB_PATH)
print("wrote", NB_PATH)
