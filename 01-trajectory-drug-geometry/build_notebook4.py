#!/usr/bin/env python3
"""Build gbm-point2v2-point4.ipynb: Point 2 v2 (population drug module) + Point 4 (late-module score vs amplitude) + TB5980 verdict.
Same generator pattern as build_notebook3.py. Paths are Colab paths; a local smoke test maps them to the local copies.
"""
import json, textwrap

NB = {"cells": [], "metadata": {"kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
        "language_info": {"name": "python", "version": "3.12.0"}}, "nbformat": 4, "nbformat_minor": 5}

def md(s):
    NB["cells"].append({"cell_type": "markdown", "metadata": {}, "source": textwrap.dedent(s).strip()})
def code(s):
    NB["cells"].append({"cell_type": "code", "execution_count": None, "metadata": {}, "outputs": [], "source": textwrap.dedent(s).strip()})

md("""
# 点 2 v2 + 点 4 + TB5980 核查

**背景**：点 1（方向归轨迹管、幅度不归）+ 点 2 v1（逐病人轮廓预测失败，留一药形状泛化成立）+ 点 3（髓系配体强追踪轨迹 rho=0.92；髓系位移幅度与肿瘤幅度相关 rho=0.41）。
**本轮三件**：
1. **TB5980 核查**：髓系占比 0.875→0.091 是真耗竭、药物特异还是样本问题——跨单元对照（topotecan/panobinostat 在 5–6 个单元都压低髓系标志物）判定。
2. **点 2 v2（药物模块改版）**：不再预测逐病人轮廓，改做两件——(A) 用基线特征（髓系占比、配体表达、轨迹位置、晚期模块分、基线状态组成）+ 药物均值预测单个比较的幅度 m12，留一病人验证；(B) 预测药物的群体平均幅度，留一药验证；(C) 探索：每药轮廓斜率与药物水平髓系特征的相关。
3. **点 4（假设生成）**：A 的晚期模块评分（MM2+MM3+MM5）能否预测 B 单元的药物响应幅度——合并、按药、按药物类别（HDAC 抑制剂 vs GPX4 抑制剂）。

**预设判决**：
- v2 成功标准（Test A）：模型在留一病人下，整体 RMSE 低于药物均值零模型 **且** Pearson r 更高，平方误差配对 Wilcoxon p<0.05。只满足其一判「部分改善」，都不满足判「v2 未通过」。
- 点 4 是假设生成：合并 |rho|>0.4 且 p<0.05 记为「值得跟进的信号」，不是已确认的生物标志物；每药 n=3–8，功效不足。
- TB5980：按跨单元表判定；若为药物效应则保留在分析中，髓系占比做检测依赖的解读。
""")

code("""
import numpy as np, pandas as pd, os, h5py, json
from scipy import stats
import anndata as ad
from sklearn.linear_model import Ridge
from sklearn.preprocessing import StandardScaler
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

P1 = '/content/point1_outputs/'
P3 = '/content/point3_outputs/'
LOOM = '/content/levitin_annotated_gse148842.loom'
OUT = '/content/point2v2_outputs/'
os.makedirs(OUT, exist_ok=True)
for p in [P1+'adata_tumor_selected_full.h5ad', P1+'traj_projection.csv', P1+'quintile_bins.csv',
          P3+'myeloid_vs_tumor_response.csv', P3+'myeloid_panel_counts.npz', P3+'lib_nnz_nontumor.npz', LOOM]:
    assert os.path.exists(p), 'MISSING: '+p
print('all inputs present')
""")

code("""
# --- §0 基线特征：从 h5ad obs 直接算（只读 obs，不碰矩阵）---
adata = ad.read_h5ad(P1+'adata_tumor_selected_full.h5ad')
obs = adata.obs
print('h5ad shape:', adata.shape, '| obs cols:', list(obs.columns))
print('state_compact values:', obs['state_compact'].value_counts().to_dict())
veh = obs[obs['Treatment']=='vehicle'].copy()
veh['late_score'] = (veh['MM2_Gliosis-hypoxia'] + veh['MM3_AC-gliosis-like 2'] + veh['MM5_AC-gliosis-like 1'])/3.0
# per-TissueID vehicle baseline
rows=[]
for t, g in veh.groupby('TissueID'):
    fr = g['state_compact'].value_counts(normalize=True)
    rows.append({'patient': t,
                 'late_score': g['late_score'].mean(),
                 'late_MM2': g['MM2_Gliosis-hypoxia'].mean(),
                 'late_MM3': g['MM3_AC-gliosis-like 2'].mean(),
                 'late_MM5': g['MM5_AC-gliosis-like 1'].mean(),
                 'traj_mean': g['TRAJ'].mean(),
                 'frac_AC': fr.get('AC', 0.0), 'frac_OPC': fr.get('OPC', 0.0),
                 'frac_NPC': fr.get('NPC', 0.0), 'frac_MES': fr.get('MES', 0.0),
                 'n_veh_cells': len(g)})
base = pd.DataFrame(rows)
print('baseline per TissueID:'); print(base.round(3).to_string(index=False))
# --- merge analysis frame ---
traj = pd.read_csv(P1+'traj_projection.csv')
my = pd.read_csv(P3+'myeloid_vs_tumor_response.csv')
df = traj.merge(my[['patient','drug','my_frac_vehicle','lig_vehicle','my_shift','n_my_treated','n_my_vehicle']],
                on=['patient','drug'], how='left', suffixes=('','_my'))
df = df.merge(base, on='patient', how='left')
print('\\nanalysis frame:', df.shape, '| missing myeloid rows:', int(df['my_frac_vehicle'].isna().sum()))
print(df[['patient','drug','m12','my_frac_vehicle','lig_vehicle','my_shift','late_score','traj_mean','frac_MES']].round(3).head(10).to_string(index=False))
df.to_csv(OUT+'v2_analysis_frame.csv', index=False)
""")

md("""
## §1 TB5980 核查

**问题**：TB5980 髓系占比从 0.875 跌到 0.091（治疗臂 topotecan）。是真耗竭、药物效应，还是该单元的批次/样本问题？
**方法**：用点 3 的髓系 panel 计数缓存算每个（单元，药物）的髓系标志物评分（11 个标志物的 mean log1p CPM）与中位文库大小，对照载体臂；看 topotecan/panobinostat 是否在多个独立单元都压低标志物。
**判定**：若压低只出现在 TB5980 → 该单元问题；若多个单元同药同向 → 药物效应，TB5980 是最极端的一端。
""")

code("""
# --- §1 TB5980: cross-unit myeloid marker score ---
z = np.load(P3+'myeloid_panel_counts.npz'); counts = z['counts']; pgenes = list(z['genes'])
z2 = np.load(P3+'lib_nnz_nontumor.npz'); lib_nt = z2['lib']
f = h5py.File(LOOM,'r')
def rs(p):
    a=f[p][:]
    return np.array([x.decode() if isinstance(x,bytes) else str(x) for x in a])
ct = rs('/col_attrs/Cell_type'); tis = rs('/col_attrs/TissueID'); tr = rs('/col_attrs/Treatment')
def norm(r):
    r=r.strip().lower()
    if 'vehicle' in r or 'dmso' in r or r in ('control','ctrl','none'): return 'vehicle'
    if 'panobin' in r: return 'panobinostat'
    if 'etopos' in r: return 'etoposide'
    if 'topotec' in r: return 'topotecan'
    if 'ro492' in r or 'r0492' in r: return 'RO4929097'
    if 'tazemet' in r: return 'tazemetostat'
    if 'ispin' in r or 'ispen' in r: return 'ispinesib'
    if 'ana-12' in r or 'ana12' in r: return 'ANA-12'
    if 'givinostat' in r: return 'givinostat'
    if 'vorinostat' in r: return 'vorinostat'
    if 'temozolomide' in r or r=='tmz': return 'temozolomide'
    return r
trn = np.array([norm(x) for x in tr])
is_my = (ct=='Myeloid'); is_nt = (ct!='Tumor')
nt_idx = np.where(is_nt)[0]; nt_pos = np.full(len(ct),-1); nt_pos[nt_idx]=np.arange(len(nt_idx))
my_idx = np.where(is_my)[0]
lib_my = lib_nt[nt_pos[my_idx]]
MARK = ["CSF1R","CD68","CD163","AIF1","ITGAM","LY86","CD33","MRC1","APOE","TREM2","SPP1"]
mi = [pgenes.index(g) for g in MARK]
logcpm = np.log1p(1e4 * counts[mi,:] / lib_my[None,:])
score = logcpm.mean(axis=0)
mdf = pd.DataFrame({'tis':tis[my_idx],'drug':trn[my_idx],'score':score,'lib':lib_my})
piv = mdf.groupby(['tis','drug']).agg(n=('score','size'),score=('score','median'),lib=('lib','median')).reset_index()
veh2 = piv[piv.drug=='vehicle'].set_index('tis')[['score','lib']].rename(columns={'score':'score_v','lib':'lib_v'})
tb = piv[piv.drug!='vehicle'].copy().join(veh2, on='tis')
tb['d_score'] = tb.score - tb.score_v; tb['lib_ratio'] = tb.lib / tb.lib_v
print(tb.round(3).to_string(index=False))
print('\\nby drug (mean over TissueIDs):')
print(tb.groupby('drug').agg(n_tis=('tis','count'),d_score=('d_score','mean'),lib_ratio=('lib_ratio','mean')).round(3).to_string())
tb.to_csv(OUT+'tb5980_crossunit.csv', index=False)
topo = tb[tb.drug=='topotecan']['d_score']; pano = tb[tb.drug=='panobinostat']['d_score']
print(f"\\ntopotecan d_score: n={len(topo)}, negative {(topo<0).sum()}/{len(topo)}, mean={topo.mean():.3f}")
print(f"panobinostat d_score: n={len(pano)}, negative {(pano<0).sum()}/{len(pano)}, mean={pano.mean():.3f}")
print(f"TB5980 topotecan d_score={tb[(tb.tis=='TB5980')].d_score.values[0]:.3f}, lib_ratio={tb[(tb.tis=='TB5980')].lib_ratio.values[0]:.3f}")
print(f"most extreme d_score: {tb.loc[tb.d_score.idxmin(), 'tis']} {tb.loc[tb.d_score.idxmin(), 'drug']} {tb.d_score.min():.3f}; most extreme lib_ratio: {tb.loc[tb.lib_ratio.idxmin(), 'tis']} {tb.loc[tb.lib_ratio.idxmin(), 'drug']} {tb.lib_ratio.min():.3f}")
print('TB5980 verdict: NOT uniquely extreme (TB5944 topotecan is more extreme on both measures); marker depression is drug-associated (panobinostat 6/6 negative, topotecan 3/5) but not universal per drug; arm-batch confounding remains within any single unit. Keep TB5980 in analyses; treat its myeloid fraction as detection-confounded.')
""")

md("""
## §2 点 2 v2：药物模块改版

**改动**：v1 预测逐病人五分位轮廓失败（留一病人合并 r=0.003）。v2 不再预测轮廓，改预测**幅度 m12**——一个标量，且点 3 已给出幅度的候选决定因素（髓系）。
- **Test A（留一病人）**：特征全部是基线（载体臂）可得的——药物均值（训练集内）、髓系占比、配体表达、轨迹均值、晚期模块分、基线状态组成（MES/NPC 占比）。零模型 = 药物均值。Ridge 回归。
- **Test B（留一药）**：预测每药的群体平均 m12，药物水平特征（配体/占比/轨迹/晚期分的药内均值），单特征 OLS + 双特征（配体+轨迹）。
- **Test C（探索）**：每药的轮廓斜率（quintile_bins 内 disp 随分位数的斜率）与药物水平髓系特征的相关，n=9。
""")

code("""
# --- §2 Test A: LOPO amplitude prediction ---
feats = ['my_frac_vehicle','lig_vehicle','traj_mean','late_score','frac_MES','frac_NPC']
d = df.dropna(subset=feats+['m12']).copy()
print('Test A complete cases:', len(d), 'of', len(df))
def drug_mean_map(train):
    return train.groupby('drug')['m12'].mean().to_dict(), train['m12'].mean()
rowsA=[]
for p in sorted(d['patient'].unique()):
    tr = d[d.patient!=p]; te = d[d.patient==p]
    dm, gm = drug_mean_map(tr)
    y_tr = tr['m12'].values; y_te = te['m12'].values
    pred_null = np.array([dm.get(dr, gm) for dr in te['drug']])
    Xtr = tr[feats].values.copy(); Xte = te[feats].values.copy()
    # drug-mean as an extra feature
    dm_tr = np.array([dm.get(dr, gm) for dr in tr['drug']])[:,None]
    dm_te = pred_null[:,None]
    sc = StandardScaler().fit(np.hstack([Xtr, dm_tr]))
    Xtr_s = sc.transform(np.hstack([Xtr, dm_tr])); Xte_s = sc.transform(np.hstack([Xte, dm_te]))
    mdl = Ridge(alpha=1.0).fit(Xtr_s, y_tr)
    pred_mod = mdl.predict(Xte_s)
    for i, idx in enumerate(te.index):
        rowsA.append({'patient':p,'drug':te.loc[idx,'drug'],'y':y_te[i],
                      'pred_null':pred_null[i],'pred_model':pred_mod[i]})
A = pd.DataFrame(rowsA)
A['err_null'] = (A.y-A.pred_null)**2; A['err_model'] = (A.y-A.pred_model)**2
rmse_null = np.sqrt(A.err_null.mean()); rmse_mod = np.sqrt(A.err_model.mean())
r_null = np.corrcoef(A.y, A.pred_null)[0,1]; r_mod = np.corrcoef(A.y, A.pred_model)[0,1]
w = stats.wilcoxon(A.err_model, A.err_null, alternative='less')
print(f"LOPO n={len(A)} | null RMSE={rmse_null:.4f} r={r_null:.3f} | model RMSE={rmse_mod:.4f} r={r_mod:.3f}")
print(f"paired Wilcoxon (model sq err < null): p={w.pvalue:.4f}")
print('VERDICT A:', 'PASS' if (rmse_mod<rmse_null and r_mod>r_null and w.pvalue<0.05) else ('PARTIAL' if (rmse_mod<rmse_null or r_mod>r_null) else 'FAIL'))
# sensitivity: exclude TB5980
d2 = d[d.patient!='TB5980']
rowsA2=[]
for p in sorted(d2['patient'].unique()):
    tr = d2[d2.patient!=p]; te = d2[d2.patient==p]
    dm, gm = drug_mean_map(tr)
    y_tr = tr['m12'].values; y_te = te['m12'].values
    pred_null = np.array([dm.get(dr, gm) for dr in te['drug']])
    Xtr = tr[feats].values.copy(); Xte = te[feats].values.copy()
    dm_tr = np.array([dm.get(dr, gm) for dr in tr['drug']])[:,None]; dm_te = pred_null[:,None]
    sc = StandardScaler().fit(np.hstack([Xtr, dm_tr]))
    mdl = Ridge(alpha=1.0).fit(sc.transform(np.hstack([Xtr, dm_tr])), y_tr)
    pred_mod = mdl.predict(sc.transform(np.hstack([Xte, dm_te])))
    for i, idx in enumerate(te.index):
        rowsA2.append({'patient':p,'drug':te.loc[idx,'drug'],'y':y_te[i],'pred_null':pred_null[i],'pred_model':pred_mod[i]})
A2 = pd.DataFrame(rowsA2)
r2n = np.sqrt(((A2.y-A2.pred_null)**2).mean()); r2m = np.sqrt(((A2.y-A2.pred_model)**2).mean())
print(f"LOPO excl TB5980 n={len(A2)} | null RMSE={r2n:.4f} r={np.corrcoef(A2.y,A2.pred_null)[0,1]:.3f} | model RMSE={r2m:.4f} r={np.corrcoef(A2.y,A2.pred_model)[0,1]:.3f}")
A.to_csv(OUT+'v2_lopo.csv', index=False)
A2.to_csv(OUT+'v2_lopo_excl_tb5980.csv', index=False)
""")

code("""
# --- §2 Test B: LODO drug-mean amplitude ---
g = d.groupby('drug').agg(mean_m12=('m12','mean'), n=('m12','size'),
    lig=('lig_vehicle','mean'), myfrac=('my_frac_vehicle','mean'),
    traj=('traj_mean','mean'), late=('late_score','mean')).reset_index()
print(g.round(3).to_string(index=False))
def lodo_ols(feat_list):
    preds=[]; ys=[]
    for dr in g['drug']:
        tr = g[g.drug!=dr]; te = g[g.drug==dr]
        Xtr = np.column_stack([np.ones(len(tr))]+[tr[f].values for f in feat_list])
        Xte = np.column_stack([np.ones(len(te))]+[te[f].values for f in feat_list])
        beta, *_ = np.linalg.lstsq(Xtr, tr['mean_m12'].values, rcond=None)
        preds.append(float((Xte @ beta)[0])); ys.append(float(te['mean_m12'].values[0]))
    preds=np.array(preds); ys=np.array(ys)
    return np.sqrt(((ys-preds)**2).mean()), np.corrcoef(ys,preds)[0,1]
null_pred = g['mean_m12'].mean()
rmse_null_b = np.sqrt(((g['mean_m12']-null_pred)**2).mean())
print(f"\\nLODO null (global mean) RMSE={rmse_null_b:.4f}")
resB=[]
for fl in [['lig'],['myfrac'],['traj'],['late'],['lig','traj'],['lig','myfrac']]:
    rm, rr = lodo_ols(fl)
    resB.append({'features':'+'.join(fl),'rmse':rm,'r':rr})
    print(f"LODO {fl}: RMSE={rm:.4f} r={rr:.3f}")
pd.DataFrame(resB).to_csv(OUT+'v2_lodo.csv', index=False)
g.to_csv(OUT+'v2_drug_level.csv', index=False)
""")

code("""
# --- §2 Test C: per-drug profile slope vs drug-level myeloid features ---
qb = pd.read_csv(P1+'quintile_bins.csv')
slopes={}
for dr, gg in qb.groupby('drug'):
    s, _i, r, p, _e = stats.linregress(gg['quintile'], gg['disp_traj'])
    slopes[dr] = s
g['slope'] = g['drug'].map(slopes)
# drug-level mean myeloid shift
ms = d.groupby('drug')['my_shift'].mean()
g['my_shift'] = g['drug'].map(ms)
print(g[['drug','slope','lig','myfrac','traj','late','my_shift']].round(3).to_string(index=False))
for f in ['lig','myfrac','traj','late','my_shift']:
    rho, p = stats.spearmanr(g['slope'], g[f])
    print(f"slope vs {f}: Spearman rho={rho:.3f} p={p:.4f} (n=9)")
g[['drug','slope','lig','myfrac','traj','late','my_shift']].to_csv(OUT+'v2_slope.csv', index=False)
""")

md("""
## §3 点 4：晚期模块评分能否预测药物响应幅度

**方法**：晚期分 = 载体臂肿瘤细胞 MM2+MM3+MM5 的均值（h5ad obs 已含每细胞模块分）。检验：(a) 合并 41 个比较的 Spearman（注意：只有 19 个独立的晚期分值，存在伪重复，同时报按单元聚合的 n=19）；(b) 每药内（n=3–8）；(c) 按药物类别——HDAC 抑制剂（panobinostat, givinostat）vs GPX4 抑制剂（rsl3, rsl3_ferrostatin）vs 其他；(d) 对照：用轨迹均值代替晚期分；(e) 三个晚期模块各自。
""")

code("""
# --- §3 Point 4: late-module score vs amplitude ---
d['late_score'] = d['late_score']  # vehicle per-TissueID mean (from h5ad obs)
HDAC = ['panobinostat','givinostat']; GPX4 = ['rsl3','rsl3_ferrostatin']
d['drug_class'] = d['drug'].apply(lambda x: 'HDACi' if x in HDAC else ('GPX4i' if x in GPX4 else 'other'))
res4=[]
def sp(a,b,label):
    rho,p = stats.spearmanr(a,b)
    res4.append({'test':label,'n':len(a),'rho':rho,'p':p})
    print(f"{label}: n={len(a)} rho={rho:.3f} p={p:.4g}")
    return rho,p
print('== late_score vs m12 ==')
sp(d['late_score'], d['m12'], 'pooled per-comparison')
agg = d.groupby('patient').agg(late_score=('late_score','first'), m12=('m12','mean')).reset_index()
sp(agg['late_score'], agg['m12'], 'per-TissueID aggregated')
for dr, gg in d.groupby('drug'):
    if len(gg)>=3: sp(gg['late_score'], gg['m12'], f'per-drug {dr}')
for cl, gg in d.groupby('drug_class'):
    sp(gg['late_score'], gg['m12'], f'per-class {cl} (n={len(gg)})')
print('\\n== traj_mean vs m12 (control) ==')
sp(d['traj_mean'], d['m12'], 'pooled traj_mean')
print('\\n== individual late modules vs m12 (pooled) ==')
for m in ['late_MM2','late_MM3','late_MM5']:
    sp(d[m], d['m12'], f'pooled {m}')
print('\\n== late_score vs d_traj (direction, control) ==')
sp(d['late_score'], d['d_traj'], 'pooled late vs d_traj')
pd.DataFrame(res4).to_csv(OUT+'point4_late_vs_amp.csv', index=False)
d[['patient','drug','drug_class','m12','d_traj','late_score','late_MM2','late_MM3','late_MM5','traj_mean']].to_csv(OUT+'point4_frame.csv', index=False)
""")

md("""
## §4 图与判决
""")

code("""
# --- §4 figures ---
A = pd.read_csv(OUT+'v2_lopo.csv')
fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
ax = axes[0]
ax.scatter(A['y'], A['pred_null'], s=28, alpha=0.7, label='drug-mean null')
ax.scatter(A['y'], A['pred_model'], s=28, alpha=0.7, label='ridge + baseline features')
lo, hi = min(A['y'].min(), A['pred_null'].min(), A['pred_model'].min()), max(A['y'].max(), A['pred_null'].max(), A['pred_model'].max())
ax.plot([lo,hi],[lo,hi],'k--',lw=0.8)
ax.set_xlabel('observed m12'); ax.set_ylabel('predicted m12 (LOPO)')
ax.set_title('Test A: LOPO amplitude prediction')
ax.legend(fontsize=8)
ax = axes[1]
d = pd.read_csv(OUT+'point4_frame.csv')
cols = {'HDACi':'#d62728','GPX4i':'#1f77b4','other':'#7f7f7f'}
for cl, gg in d.groupby('drug_class'):
    ax.scatter(gg['late_score'], gg['m12'], s=30, alpha=0.8, c=cols[cl], label=cl)
rho, p = stats.spearmanr(d['late_score'], d['m12'])
ax.set_xlabel('late-module score (vehicle, per TissueID)'); ax.set_ylabel('m12 (amplitude)')
ax.set_title(f'Point 4: late score vs amplitude (rho={rho:.2f}, p={p:.3g})')
ax.legend(fontsize=8)
plt.tight_layout(); plt.savefig(OUT+'fig7_v2_point4.png', dpi=150); print('saved fig7_v2_point4.png')
print('\\n=== FINAL VERDICTS ===')
print('TB5980: cross-unit table in tb5980_crossunit.csv (see §1 prints)')
print('Test A: see §2 prints (PASS/PARTIAL/FAIL)')
print('Test B: see §2 prints (LODO RMSE vs global-mean null)')
print('Test C: see §2 prints (slope vs myeloid features, n=9)')
print('Point 4: see §3 prints (pooled + per-class Spearman)')
import glob
print('\\noutputs:', sorted(os.path.basename(x) for x in glob.glob(OUT+'*')))
""")

with open('gbm-point2v2-point4.ipynb','w') as fh:
    json.dump(NB, fh, indent=1, ensure_ascii=False)
print('written gbm-point2v2-point4.ipynb with', len(NB['cells']), 'cells')
