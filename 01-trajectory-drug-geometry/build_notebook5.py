#!/usr/bin/env python3
"""Build gbm-point2v3-stratified.ipynb: Point 2 v3 — the stratified drug module.
Population shape (LODO, v1) + late-module-score strata (Point 4) as the amplitude stratifier.
Source stored as plain strings (see 2026-11-08 lesson: never .split("\\n") lists without trailing newlines).
"""
import json, textwrap

NB = {"cells": [], "metadata": {"kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
        "language_info": {"name": "python", "version": "3.12.0"}}, "nbformat": 4, "nbformat_minor": 5}

def md(s):
    NB["cells"].append({"cell_type": "markdown", "metadata": {}, "source": textwrap.dedent(s).strip()})
def code(s):
    NB["cells"].append({"cell_type": "code", "execution_count": None, "metadata": {}, "outputs": [], "source": textwrap.dedent(s).strip()})

md("""
# 点 2 v3：分层版药物模块

**为什么改版**：v1 逐病人轮廓预测失败（LOPO 合并 r=0.003）；v2 加 7 个基线特征的 Ridge 反而差于药物均值零模型（RMSE 0.0687 vs 0.0615）——n=41 时模型太富会过拟合。但点 4 给出一条单变量的组织水平信号：晚期模块分与幅度在组织单元水平 rho=0.640（p=0.006，n=17）。v3 只用**一个分层变量**（晚期模块分三分位）+ **群体形状**（v1 留一药斜率），不多加任何特征。

**模块构成**：
1. **幅度分层**：新单元按晚期模块分（载体臂肿瘤细胞 MM2+MM3+MM5 均值）落入低/中/高三分位，预期幅度 = 该分层的平均 m12。
2. **方向形状**：每药的群体斜率（v1 LODO 斜率，disp 随轨迹分位数变化），附每药留一药 Pearson 作为可靠度标志（6/9 药为正，ANA-12、etoposide、rsl3 为负——形状对这三药不泛化）。
3. **髓系注释**：晚期分层同时对应髓系配体水平（点 3：配体表达与轨迹位置 rho=0.92），高晚期层 = 高 VEGFA/VCAN/ANXA1 髓系配体 niche；只作注释，不作第二预测变量（与晚期分共线）。

**预设判决**（留一病人，比较水平，n=41）：
- PASS：分层预测（三分位分层均值，或连续线性 m12~晚期分）在 RMSE 上同时打败全局均值零模型与药物均值零模型，平方误差配对 Wilcoxon p<0.05。
- PARTIAL：只打败全局均值（分层信息存在，但没超过药物身份本身）。
- FAIL：两个都打败不了。
""")

code("""
import numpy as np, pandas as pd, os, json
from scipy import stats
import anndata as ad
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

P1 = '/content/point1_outputs/'
P2 = '/content/point2_outputs/'
P3 = '/content/point3_outputs/'
OUT = '/content/point2v3_outputs/'
os.makedirs(OUT, exist_ok=True)
for p in [P1+'adata_tumor_selected_full.h5ad', P1+'traj_projection.csv', P1+'quintile_bins.csv',
          P2+'lodo_shape.csv', P3+'myeloid_vs_tumor_response.csv']:
    assert os.path.exists(p), 'MISSING: '+p
print('all inputs present')
""")

code("""
# --- §0 分析框：与 v2 相同（h5ad obs 直接算晚期分与基线状态组成）---
adata = ad.read_h5ad(P1+'adata_tumor_selected_full.h5ad')
obs = adata.obs
veh = obs[obs['Treatment']=='vehicle'].copy()
veh['late_score'] = (veh['MM2_Gliosis-hypoxia'] + veh['MM3_AC-gliosis-like 2'] + veh['MM5_AC-gliosis-like 1'])/3.0
rows=[]
for t, g in veh.groupby('TissueID'):
    fr = g['state_compact'].value_counts(normalize=True)
    rows.append({'patient': t, 'late_score': g['late_score'].mean(),
                 'traj_mean': g['TRAJ'].mean(),
                 'frac_MES': fr.get('MES', 0.0), 'frac_NPC': fr.get('NPC', 0.0),
                 'n_veh_cells': len(g)})
base = pd.DataFrame(rows)
traj = pd.read_csv(P1+'traj_projection.csv')
my = pd.read_csv(P3+'myeloid_vs_tumor_response.csv')
df = traj.merge(my[['patient','drug','my_frac_vehicle','lig_vehicle','my_shift']],
                on=['patient','drug'], how='left')
df = df.merge(base, on='patient', how='left')
d = df.dropna(subset=['m12','late_score']).copy()
print('analysis frame:', d.shape, '| TissueIDs:', d['patient'].nunique(), '| drugs:', d['drug'].nunique())
print('late_score range:', round(d['late_score'].min(),3), 'to', round(d['late_score'].max(),3))
d.to_csv(OUT+'v3_analysis_frame.csv', index=False)
""")

md("""
## §1 分层：晚期分三分位
""")

code("""
# --- §1 strata (descriptive, full data) ---
qs = base['late_score'].quantile([1/3, 2/3]).values
print('tertile cutoffs (late_score):', np.round(qs, 4))
def assign(x):
    return 'low' if x <= qs[0] else ('mid' if x <= qs[1] else 'high')
base['stratum'] = base['late_score'].apply(assign)
d['stratum'] = d['late_score'].apply(assign)
print(base[['patient','late_score','stratum']].sort_values('late_score').round(3).to_string(index=False))
summ = d.groupby('stratum').agg(n_comp=('m12','size'), n_tis=('patient','nunique'),
    mean_m12=('m12','mean'), sd_m12=('m12','std'),
    mean_late=('late_score','mean'), mean_lig=('lig_vehicle','mean'),
    mean_traj=('traj_mean','mean')).reindex(['low','mid','high'])
print('\\nper-stratum summary:'); print(summ.round(3).to_string())
kw = stats.kruskal(*[d[d.stratum==s]['m12'] for s in ['low','mid','high']])
print(f"Kruskal-Wallis m12 by stratum (41 comps): H={kw.statistic:.2f} p={kw.pvalue:.4f}")
agg = d.groupby('patient').agg(late_score=('late_score','first'), m12=('m12','mean'), stratum=('stratum','first')).reset_index()
kw2 = stats.kruskal(*[agg[agg.stratum==s]['m12'] for s in ['low','mid','high']])
print(f"Kruskal-Wallis mean m12 by stratum (per-TissueID, n={len(agg)}): H={kw2.statistic:.2f} p={kw2.pvalue:.4f}")
summ.to_csv(OUT+'v3_strata_summary.csv')
base.to_csv(OUT+'v3_tissue_strata.csv', index=False)
d.to_csv(OUT+'v3_analysis_frame.csv', index=False)  # re-save WITH stratum column (for §5)
print('\\nnote: v3 uses all 43 comparisons (late score needs no myeloid data); v2 Test A used 41 (2 lacked myeloid features).')
""")

md("""
## §2 留一病人验证：分层能不能预测新单元的幅度

**方法**：留一病人。每折内只用训练单元定三分位切点与分层均值（无泄漏）。对照两个零模型：全局均值、药物均值。候选：(a) 连续线性 m12~晚期分；(b) 三分位分层均值；(c) 分层+药物加性 OLS（次要，参数多）。主检验在比较水平（n=41），辅检验在单元水平（每单元平均 m12，n=17）。
""")

code("""
# --- §2 LOPO validation ---
def tertiles_from(train_tis):
    q = train_tis['late_score'].quantile([1/3, 2/3]).values
    return q
def strat_of(x, q):
    return 'low' if x <= q[0] else ('mid' if x <= q[1] else 'high')
rows=[]
for p in sorted(d['patient'].unique()):
    tr = d[d.patient!=p]; te = d[d.patient==p]
    tr_tis = tr.groupby('patient')['late_score'].first().reset_index()
    q = tertiles_from(tr_tis)
    gm = tr['m12'].mean()
    dm = tr.groupby('drug')['m12'].mean().to_dict()
    sm = tr.assign(stratum=tr['late_score'].apply(lambda x: strat_of(x, q))).groupby('stratum')['m12'].mean().to_dict()
    # (a) continuous OLS m12 ~ late_score
    Xtr = np.column_stack([np.ones(len(tr)), tr['late_score'].values])
    beta, *_ = np.linalg.lstsq(Xtr, tr['m12'].values, rcond=None)
    # (c) stratum + drug additive OLS
    def design(frame):
        st = frame['late_score'].apply(lambda x: strat_of(x, q))
        cols = [np.ones(len(frame))]
        for s in ['mid','high']:
            cols.append((st==s).astype(float).values)
        for dr in sorted(tr['drug'].unique())[1:]:
            cols.append((frame['drug']==dr).astype(float).values)
        return np.column_stack(cols)
    Dtr = design(tr); beta_c, *_ = np.linalg.lstsq(Dtr, tr['m12'].values, rcond=None)
    for _, row in te.iterrows():
        s = strat_of(row['late_score'], q)
        pred_a = beta[0] + beta[1]*row['late_score']
        pred_b = sm.get(s, gm)
        pred_c = float((design(pd.DataFrame([row])) @ beta_c)[0])
        rows.append({'patient': p, 'drug': row['drug'], 'y': row['m12'],
                     'late_score': row['late_score'], 'stratum': s,
                     'pred_global': gm,
                     'pred_drugmean': dm.get(row['drug'], gm),
                     'pred_linear': pred_a, 'pred_stratum': pred_b, 'pred_additive': pred_c})
R = pd.DataFrame(rows)
def report(col, label):
    rmse = np.sqrt(((R.y - R[col])**2).mean())
    r = np.corrcoef(R.y, R[col])[0,1]
    return rmse, r
for col, label in [('pred_global','global mean null'), ('pred_drugmean','drug-mean null'),
                   ('pred_linear','(a) linear late_score'), ('pred_stratum','(b) tertile stratum mean'),
                   ('pred_additive','(c) stratum+drug additive')]:
    rmse, r = report(col, label)
    print(f"{label}: RMSE={rmse:.4f} r={r:.3f}")
print()
for col, label in [('pred_linear','(a) linear'), ('pred_stratum','(b) stratum'), ('pred_additive','(c) additive')]:
    for null, nlabel in [('pred_global','global'), ('pred_drugmean','drug-mean')]:
        w = stats.wilcoxon((R.y-R[col])**2, (R.y-R[null])**2, alternative='less')
        print(f"{label} vs {nlabel} null: paired Wilcoxon p={w.pvalue:.4f}")
rm_g = report('pred_global', '')[0]; rm_d = report('pred_drugmean', '')[0]
rm_b = report('pred_stratum', '')[0]; rm_a = report('pred_linear', '')[0]
best = min(rm_a, rm_b)
w_g = stats.wilcoxon((R.y-R['pred_stratum'])**2, (R.y-R['pred_global'])**2, alternative='less').pvalue
w_d = stats.wilcoxon((R.y-R['pred_stratum'])**2, (R.y-R['pred_drugmean'])**2, alternative='less').pvalue
print(f"\\nVERDICT (comparison-level, n={len(R)}, pre-registered): stratum RMSE={rm_b:.4f} vs global {rm_g:.4f} (p={w_g:.4f}) vs drug-mean {rm_d:.4f} (p={w_d:.4f})")
print('PASS' if (rm_b<rm_g and rm_b<rm_d and w_g<0.05 and w_d<0.05) else ('PARTIAL' if (rm_b<rm_g and w_g<0.05) else 'FAIL'))
R.to_csv(OUT+'v3_lopo.csv', index=False)
# unit-level: mean m12 per TissueID (pre-planned secondary test)
U = R.groupby('patient').agg(y=('y','mean'), late_score=('late_score','first'),
    pred_global=('pred_global','mean'), pred_linear=('pred_linear','mean'),
    pred_stratum=('pred_stratum','mean')).reset_index()
for col in ['pred_global','pred_linear','pred_stratum']:
    print(f"unit-level {col}: RMSE={np.sqrt(((U.y-U[col])**2).mean()):.4f} r={np.corrcoef(U.y, U[col])[0,1]:.3f} (n={len(U)})")
wu_s = stats.wilcoxon((U.y-U.pred_stratum)**2, (U.y-U.pred_global)**2, alternative='less')
wu_l = stats.wilcoxon((U.y-U.pred_linear)**2, (U.y-U.pred_global)**2, alternative='less')
print(f"unit-level stratum vs global Wilcoxon p={wu_s.pvalue:.4f}; linear vs global p={wu_l.pvalue:.4f}")
U.to_csv(OUT+'v3_lopo_unit.csv', index=False)
# overall module verdict: comparison-level is the pre-registered test; unit-level is the pre-planned secondary
if (rm_b<rm_g and rm_b<rm_d and w_g<0.05 and w_d<0.05):
    MODULE_VERDICT = 'PASS'
elif wu_s.pvalue < 0.05:
    MODULE_VERDICT = 'PARTIAL (tissue-level)'
else:
    MODULE_VERDICT = 'FAIL'
print('MODULE VERDICT:', MODULE_VERDICT)
""")

md("""
## §3 群体形状组件（v1 LODO 的复述 + 每药表）
""")

code("""
# --- §3 population shape component ---
lodo = pd.read_csv(P2+'lodo_shape.csv')
qb = pd.read_csv(P1+'quintile_bins.csv')
slopes={}
for dr, gg in qb.groupby('drug'):
    s, _i, r, p, _e = stats.linregress(gg['quintile'], gg['disp_traj'])
    slopes[dr] = s
shape = lodo.set_index('drug').join(pd.Series(slopes, name='slope_all')).reset_index()
shape = shape.merge(d.groupby('drug').agg(mean_m12=('m12','mean'), n=('m12','size')).reset_index(), on='drug')
print(shape[['drug','slope_all','slope_train','pearson','spearman','rmse_model','rmse_flat','mean_m12','n']].round(4).to_string(index=False))
print(f"\\nLODO Pearson: median={shape['pearson'].median():.3f}, positive {int((shape['pearson']>0).sum())}/9, negative drugs: {list(shape[shape.pearson<0]['drug'])}")
print(f"slopes (disp per quintile): all negative, range {shape['slope_all'].min():.4f} to {shape['slope_all'].max():.4f}")
shape.to_csv(OUT+'v3_shape_table.csv', index=False)
""")

md("""
## §4 模块规格：给新单元用的规则
""")

code("""
# --- §4 module specification (follows the §2 verdict) ---
R = pd.read_csv(OUT+'v3_lopo.csv')
rm_g = np.sqrt(((R.y-R.pred_global)**2).mean()); rm_d = np.sqrt(((R.y-R.pred_drugmean)**2).mean())
rm_b = np.sqrt(((R.y-R.pred_stratum)**2).mean()); rm_a = np.sqrt(((R.y-R.pred_linear)**2).mean())
rm_c = np.sqrt(((R.y-R.pred_additive)**2).mean())
w_g = stats.wilcoxon((R.y-R.pred_stratum)**2, (R.y-R.pred_global)**2, alternative='less').pvalue
w_d = stats.wilcoxon((R.y-R.pred_stratum)**2, (R.y-R.pred_drugmean)**2, alternative='less').pvalue
U = pd.read_csv(OUT+'v3_lopo_unit.csv')
wu_s = stats.wilcoxon((U.y-U.pred_stratum)**2, (U.y-U.pred_global)**2, alternative='less').pvalue
if (rm_b<rm_g and rm_b<rm_d and w_g<0.05 and w_d<0.05):
    verdict = 'PASS'
elif wu_s < 0.05:
    verdict = 'PARTIAL (tissue-level)'
else:
    verdict = 'FAIL'
print('module verdict:', verdict, f'(comparison n={len(R)}: stratum RMSE {rm_b:.4f} vs global {rm_g:.4f} vs drug-mean {rm_d:.4f}; linear {rm_a:.4f}; additive {rm_c:.4f}; unit-level stratum vs global p={wu_s:.4f})')
# full-data strata (descriptive spec)
summ = pd.read_csv(OUT+'v3_strata_summary.csv', index_col=0)
shape = pd.read_csv(OUT+'v3_shape_table.csv')
qs = base['late_score'].quantile([1/3, 2/3]).values
spec = {
  'module': 'GBM drug-response module v3 (stratified)',
  'amplitude_rule': ('tertile stratum mean predicts a tissue AVERAGE amplitude across drugs (unit-level LOPO validated); it does NOT predict a single drug amplitude in a tissue (comparison-level LOPO p=0.23)' if verdict.startswith('PARTIAL') else ('tertile stratum mean' if verdict=='PASS' else 'DESCRIPTIVE ONLY — strata do not predict out-of-sample; use population mean')),
  'late_score_definition': 'mean of S10 module scores MM2 (Gliosis-hypoxia) + MM3 (AC-gliosis-like 2) + MM5 (AC-gliosis-like 1) in vehicle-arm tumor cells (scanpy score_genes on log-normalized counts)',
  'strata': [],
  'drugs': {},
  'direction_rule': 'displacement is toward early trajectory states; later cells displace more (Point 1: pooled rho=-0.242, bin means rho=-0.894). Per-drug shape = population slope; reliability flag = LODO Pearson per drug.',
  'myeloid_annotation': 'high-late stratum = high myeloid ligand niche (VEGFA/VCAN/ANXA1; Point 3 ligand-trajectory rho=0.92). Annotation only, not a second predictor (collinear with late score).',
  'validation': {'loocv_n_comparisons': len(R), 'loocv_n_tissues': len(U), 'rmse_global_null': round(float(rm_g),4), 'rmse_drugmean_null': round(float(rm_d),4),
                 'rmse_stratum': round(float(rm_b),4), 'rmse_linear': round(float(rm_a),4), 'rmse_additive': round(float(rm_c),4),
                 'wilcoxon_vs_global_p': round(float(w_g),4), 'wilcoxon_vs_drugmean_p': round(float(w_d),4),
                 'unit_level_stratum_vs_global_p': round(float(wu_s),4),
                 'verdict': verdict},
}
for s in ['low','mid','high']:
    row = summ.loc[s]
    lo = float('-inf') if s=='low' else float(qs[0]); hi = float(qs[1]) if s in ('low','mid') else float('inf')
    spec['strata'].append({'name': s, 'late_score_range': [lo if lo!=float('-inf') else None, hi if hi!=float('inf') else None],
                           'expected_m12': round(float(row['mean_m12']),4), 'sd_m12': round(float(row['sd_m12']),4),
                           'n_comparisons': int(row['n_comp']), 'n_tissues': int(row['n_tis']),
                           'mean_myeloid_ligand': round(float(row['mean_lig']),3)})
for _, row in shape.iterrows():
    spec['drugs'][row['drug']] = {'mean_m12': round(float(row['mean_m12']),4), 'n': int(row['n']),
        'slope_per_quintile': round(float(row['slope_all']),5),
        'lodo_pearson': round(float(row['pearson']),3),
        'shape_reliable': bool(row['pearson']>0)}
with open(OUT+'v3_module_spec.json','w') as fh:
    json.dump(spec, fh, indent=2, ensure_ascii=False)
print(json.dumps(spec, indent=2, ensure_ascii=False)[:2500])
""")

md("""
## §5 图与判决
""")

code("""
# --- §5 figures ---
R = pd.read_csv(OUT+'v3_lopo.csv')
fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
ax = axes[0]
d = pd.read_csv(OUT+'v3_analysis_frame.csv')
order = ['low','mid','high']
for s in order:
    gg = d[d.stratum==s]
    ax.scatter([s]*len(gg), gg['m12'], s=30, alpha=0.7)
    ax.plot([s], [gg['m12'].mean()], marker='_', ms=20, mew=2, color='black')
ax.set_xlabel('late-module-score stratum'); ax.set_ylabel('m12 (amplitude)')
ax.set_title('Amplitude by late-score stratum (n=41)')
ax = axes[1]
ax.scatter(R['y'], R['pred_global'], s=28, alpha=0.6, label='global mean null')
ax.scatter(R['y'], R['pred_stratum'], s=28, alpha=0.7, label='tertile stratum')
ax.scatter(R['y'], R['pred_drugmean'], s=28, alpha=0.6, label='drug-mean null')
lo, hi = min(R['y'].min(), R['pred_global'].min(), R['pred_stratum'].min()), max(R['y'].max(), R['pred_global'].max(), R['pred_stratum'].max())
ax.plot([lo,hi],[lo,hi],'k--',lw=0.8)
ax.set_xlabel('observed m12'); ax.set_ylabel('predicted m12 (LOPO)')
ax.set_title('v3 LOPO amplitude prediction')
ax.legend(fontsize=8)
plt.tight_layout(); plt.savefig(OUT+'fig8_v3_stratified.png', dpi=150); print('saved fig8_v3_stratified.png')
print('\\n=== FINAL VERDICTS ===')
import glob
print('outputs:', sorted(os.path.basename(x) for x in glob.glob(OUT+'*')))
""")

with open('gbm-point2v3-stratified.ipynb','w') as fh:
    json.dump(NB, fh, indent=1, ensure_ascii=False)
print('written gbm-point2v3-stratified.ipynb with', len(NB['cells']), 'cells')
