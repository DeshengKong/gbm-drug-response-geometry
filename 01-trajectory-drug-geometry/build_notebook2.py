"""Build gbm-point2-drug-module.ipynb: Point 2 prototype of the virtual-drug module.

Reads Point 1 outputs already on the Colab VM (/content/point1_outputs/):
  quintile_bins.csv  (214 rows: per-comparison 5-bin displacement/magnitude)
  traj_projection.csv (43 rows: per-comparison d_traj, m12, C, M, per-state disp)

Tests (all on held-out data, no refitting on the test fold):
  LOPO: leave-one-patient-out -- predict a patient's per-drug 5-bin
        displacement profile from the same drug in OTHER patients.
        Nulls: flat profile (drug mean) and global profile (all drugs).
  LODO: leave-one-drug-out -- predict a drug's CENTERED quintile profile
        from the global disp-vs-quintile slope learned on other drugs.
        Null: flat (zero slope). Permutation check (quintile shuffle).
  Magnitude LOPO: predict per-comparison mean magnitude from same-drug
        other-patient mean vs global mean (verdict says magnitude is flat
        along trajectory, so the model is per-drug mean, no quintile term).

Verdict rule (pre-registered in the notebook):
  trajectory-aware model beats flat null on held-out data
  (median per-comparison correlation higher AND paired Wilcoxon on RMSE p<0.05)
  in LOPO direction test -> module v1 is predictive for new patients.

Fixes vs v1 (2026-11-08 08:09): one comparison has only 4 bins (a bin with
<10 cells was dropped in Point 1), so pred/obs must align on the quintiles
actually present instead of assuming 5. Permutation now shuffles quintile
labels WITHIN each comparison (global shuffle emptied some drug x quintile
cells -> NaN). Pearson guarded against constant inputs.
"""
import nbformat as nbf

nb = nbf.v4.new_notebook()
nb.metadata = {
    "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
    "language_info": {"name": "python", "version": "3.12"},
}

md0 = """# 点 2：虚拟药物模块 v1 原型

**问题**：点 1 的 verdict（方向是轨迹位置的连续函数、幅度与位置无关）能不能变成一个可验证预测的模型？

**模型 v1（按 verdict 建）**
- 方向：同一药物的 5 箱位移轮廓（disp_traj 随轨迹五分位的变化），新病人用其他病人同药物的平均轮廓预测。
- 幅度：同一药物的平均幅度（位置无关，按 verdict 不设五分位项）。
- 跨药泛化：中心化轮廓的斜率（disp 随五分位变化的全局斜率）是否在没见过的药上也成立。

**检验（全部留出法，测试折不参与训练）**
1. **LOPO**（留一病人）：预测某病人某药物的 5 箱轮廓；对照：平坦轮廓（药物均值）、全局轮廓（所有药物平均）。
2. **LODO**（留一药）：用其他药物的全局斜率预测该药物的中心化轮廓；对照：平坦（斜率 0）。另做五分位标签置换检验。
3. **幅度 LOPO**：同药物其他病人平均幅度 vs 全局平均幅度。

**判决标准（事先定好）**：轨迹模型在留出数据上同时满足——中位相关高于平坦对照、且配对 Wilcoxon（RMSE）p<0.05——才算对新病人有预测力。
"""

c1 = """import numpy as np, pandas as pd
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

OUT = Path('/content/point2_outputs'); OUT.mkdir(exist_ok=True)
qb = pd.read_csv('/content/point1_outputs/quintile_bins.csv')
tj = pd.read_csv('/content/point1_outputs/traj_projection.csv')

def _pr(a, b):
    a = np.asarray(a, dtype=float); b = np.asarray(b, dtype=float)
    if np.std(a) == 0 or np.std(b) == 0: return float('nan')
    return stats.pearsonr(a, b).statistic

def _sp(a, b):
    a = np.asarray(a, dtype=float); b = np.asarray(b, dtype=float)
    if np.std(a) == 0 or np.std(b) == 0: return float('nan')
    return stats.spearmanr(a, b).statistic

print('quintile_bins:', qb.shape, '| comparisons:', qb[['patient','drug']].drop_duplicates().shape[0])
print('patients:', qb.patient.nunique(), '| drugs:', sorted(qb.drug.unique()))
print('per-drug comparison counts:'); print(qb[['patient','drug']].drop_duplicates().groupby('drug').size())
print('bins per comparison:'); print(qb.groupby(['patient','drug']).size().value_counts())
assert qb.shape[0] >= 200 and tj.shape[0] == 43
qb.head(3)
"""

md1 = """## §1 LOPO：留一病人——新病人的轮廓能预测吗？
对每个（病人，药物）比较：预测 = 该药物在**其他病人**身上的平均轮廓，按该比较实际有的五分位对齐（个别比较只有 4 箱）。对照：(a) 平坦轮廓 = 该药物平均位移（常数）；(b) 全局轮廓 = 所有药物、所有病人的平均轮廓。指标：每比较的 Pearson r、Spearman rho、RMSE；汇总中位数 + 全部配对点的合并相关；配对 Wilcoxon（模型 RMSE vs 平坦 RMSE）。
"""

c2 = """rows = []
for (p, d), g in qb.groupby(['patient', 'drug']):
    g = g.sort_values('quintile')
    qs = g['quintile'].to_numpy()
    obs = g['disp_traj'].to_numpy()
    train = qb[(qb.drug == d) & (qb.patient != p)]
    if train.patient.nunique() == 0:
        rows.append(dict(patient=p, drug=d, skipped='no other patient on drug'))
        continue
    pred = train.groupby('quintile')['disp_traj'].mean().reindex(qs).to_numpy()
    flat = np.full(len(qs), train['disp_traj'].mean())
    glob = qb.groupby('quintile')['disp_traj'].mean().reindex(qs).to_numpy()
    rmse = lambda a, b: float(np.sqrt(np.mean((a-b)**2)))
    rows.append(dict(patient=p, drug=d, skipped='', n_bins=len(qs),
                     pearson=_pr(pred, obs), spearman=_sp(pred, obs),
                     rmse_model=rmse(pred, obs), rmse_flat=rmse(flat, obs),
                     rmse_global=rmse(glob, obs),
                     pred=pred.tolist(), obs=obs.tolist()))
lopo = pd.DataFrame(rows)
sk = lopo[lopo.skipped != '']; lopo = lopo[lopo.skipped == ''].copy()
print(f'LOPO comparisons: {len(lopo)} evaluated, {len(sk)} skipped')
print('median Pearson r (model):', round(float(np.nanmedian(lopo.pearson)), 3),
      '| median Spearman rho:', round(float(np.nanmedian(lopo.spearman)), 3))
print('median RMSE  model / flat / global:',
      round(lopo.rmse_model.median(), 4), round(lopo.rmse_flat.median(), 4), round(lopo.rmse_global.median(), 4))
print('model beats flat on RMSE in', (lopo.rmse_model < lopo.rmse_flat).mean().round(3), 'of comparisons')
w = stats.wilcoxon(lopo.rmse_model, lopo.rmse_flat)
print('paired Wilcoxon model vs flat RMSE: p =', f'{w.pvalue:.3e}')
pred_all = np.concatenate(lopo.pred.to_numpy()); obs_all = np.concatenate(lopo.obs.to_numpy())
print('pooled Pearson r (all bins):', round(stats.pearsonr(pred_all, obs_all).statistic, 3))
lopo.to_csv(OUT/'lopo_direction.csv', index=False)
"""

md2 = """## §2 LODO：留一药——斜率能泛化到没见过的药吗？
对每个药物：用**其他药物**中心化轮廓的合并斜率（disp 随五分位线性回归）预测该药物的中心化轮廓。对照：平坦（斜率 0）。置换检验：同一流程下把五分位标签在**每个比较内部**随机打乱 500 次，看中位相关的零分布。
"""

c3 = """# per-drug mean profile (centered) + global slope from other drugs
prof = qb.groupby(['drug','quintile'])['disp_traj'].mean().unstack()
prof_c = prof.sub(prof.mean(axis=1), axis=0)  # center each drug
res = []
for d in prof.index:
    others = prof_c.drop(index=d)
    x = np.tile([1,2,3,4,5], len(others)); y = others.to_numpy().ravel()
    slope = stats.linregress(x, y).slope
    pred_c = slope * (np.array([1,2,3,4,5]) - 3)
    obs_c = prof_c.loc[d].to_numpy()
    rmse = lambda a, b: float(np.sqrt(np.mean((a-b)**2)))
    res.append(dict(drug=d, slope_train=slope,
                    pearson=_pr(pred_c, obs_c),
                    spearman=_sp(pred_c, obs_c),
                    rmse_model=rmse(pred_c, obs_c), rmse_flat=rmse(np.zeros(5), obs_c)))
lodo = pd.DataFrame(res)
print('median slope (train):', round(lodo.slope_train.median(), 5))
print('median Pearson r (centered profile):', round(float(np.nanmedian(lodo.pearson)), 3))
print('median RMSE model / flat:', round(lodo.rmse_model.median(), 4), round(lodo.rmse_flat.median(), 4))
print('model beats flat in', (lodo.rmse_model < lodo.rmse_flat).mean().round(3), 'of drugs')
print('paired Wilcoxon model vs flat RMSE: p =', f"{stats.wilcoxon(lodo.rmse_model, lodo.rmse_flat).pvalue:.3e}")
lodo.to_csv(OUT/'lodo_shape.csv', index=False)

rng = np.random.default_rng(0)
null_r = []
for _ in range(500):
    sh = qb.copy()
    sh['quintile'] = sh.groupby(['patient','drug'])['quintile'].transform(lambda s: rng.permuted(s.to_numpy()))
    pr = sh.groupby(['drug','quintile'])['disp_traj'].mean().unstack()
    pr_c = pr.sub(pr.mean(axis=1), axis=0)
    rs = []
    for d in pr.index:
        o = pr_c.drop(index=d); x = np.tile([1,2,3,4,5], len(o)); y = o.to_numpy().ravel()
        s = stats.linregress(x, y).slope
        rs.append(_pr(s*(np.array([1,2,3,4,5])-3), pr_c.loc[d].to_numpy()))
    null_r.append(np.nanmedian(rs))
null_r = np.array(null_r)
print('permutation null: median r =', round(float(np.nanmean(null_r)), 3), '| 95% CI [', round(float(np.nanquantile(null_r, .025)), 3), ',', round(float(np.nanquantile(null_r, .975)), 3), ']')
print('empirical p (observed median r >= null):', float((null_r >= np.nanmedian(lodo.pearson)).mean()))
"""

md3 = """## §3 幅度 LOPO：幅度能按药物预测吗？
按 verdict 幅度与轨迹位置无关，所以幅度模型 = 该药物在其他病人身上的平均幅度（不设五分位项）。对照：全局平均幅度。
"""

c4 = """qb['mean_mag'] = qb.groupby(['patient','drug'])['mag12'].transform('mean')
mrows = []
for (p, d), g in qb.groupby(['patient','drug']):
    obs = g['mean_mag'].iloc[0]
    train = qb[(qb.drug == d) & (qb.patient != p)]
    if train.patient.nunique() == 0:
        continue
    pred = train.groupby(['patient','drug'])['mean_mag'].first().mean()
    mrows.append(dict(patient=p, drug=d, obs=obs, pred_drug=pred, pred_global=qb['mean_mag'].mean()))
mag = pd.DataFrame(mrows)
for c in ['pred_drug', 'pred_global']:
    r = _pr(mag[c], mag.obs)
    print(c, 'RMSE:', round(float(np.sqrt(np.mean((mag[c]-mag.obs)**2))), 4),
          '| Pearson r:', 'n/a (constant)' if np.isnan(r) else round(r, 3))
print('paired Wilcoxon drug-mean vs global RMSE: p =',
      f"{stats.wilcoxon(np.abs(mag.pred_drug-mag.obs), np.abs(mag.pred_global-mag.obs)).pvalue:.3e}")
mag.to_csv(OUT/'lopo_magnitude.csv', index=False)
"""

md4 = """## §4 图与 verdict
"""

c5 = """fig, ax = plt.subplots(1, 2, figsize=(11, 4.5))
ax[0].scatter(np.concatenate(lopo.pred.to_numpy()), np.concatenate(lopo.obs.to_numpy()), s=14, alpha=.6)
lo, hi = -0.06, 0.04; ax[0].plot([lo, hi], [lo, hi], 'r--', lw=1)
ax[0].set_xlabel('predicted disp_traj (LOPO, same drug other patients)')
ax[0].set_ylabel('observed disp_traj'); ax[0].set_title('LOPO direction: predicted vs observed')
ax[0].set_xlim(lo, hi); ax[0].set_ylim(lo, hi)
x = np.arange(1, 6)
for d in prof.index:
    ax[1].plot(x, prof_c.loc[d], marker='o', ms=3, lw=1, alpha=.7, label=d)
ax[1].axhline(0, color='k', lw=.8); ax[1].set_xlabel('trajectory quintile'); ax[1].set_ylabel('centered disp_traj')
ax[1].set_title('Centered quintile profiles per drug (LODO input)'); ax[1].legend(fontsize=7, ncol=2)
plt.tight_layout(); plt.savefig(OUT/'fig5_lopo_lodo.png', dpi=150); print('saved fig5_lopo_lodo.png')

print('=== 点 2 v1 verdict ===')
print('LOPO direction: median Pearson r =', round(float(np.nanmedian(lopo.pearson)), 3),
      '| model beats flat RMSE in', round((lopo.rmse_model < lopo.rmse_flat).mean(), 2),
      'of comparisons | Wilcoxon p =', f"{w.pvalue:.2e}")
print('LODO shape: median Pearson r =', round(float(np.nanmedian(lodo.pearson)), 3),
      '| model beats flat in', round((lodo.rmse_model < lodo.rmse_flat).mean(), 2), 'of drugs')
print('Magnitude LOPO: drug-mean RMSE =', round(float(np.sqrt(np.mean((mag.pred_drug-mag.obs)**2))), 4),
      'vs global', round(float(np.sqrt(np.mean((mag.pred_global-mag.obs)**2))), 4))
print('outputs in', OUT)
"""

nb.cells = [
    nbf.v4.new_markdown_cell(md0),
    nbf.v4.new_code_cell(c1),
    nbf.v4.new_markdown_cell(md1),
    nbf.v4.new_code_cell(c2),
    nbf.v4.new_markdown_cell(md2),
    nbf.v4.new_code_cell(c3),
    nbf.v4.new_markdown_cell(md3),
    nbf.v4.new_code_cell(c4),
    nbf.v4.new_markdown_cell(md4),
    nbf.v4.new_code_cell(c5),
]
nbf.write(nb, 'gbm-point2-drug-module.ipynb')
print('written')
