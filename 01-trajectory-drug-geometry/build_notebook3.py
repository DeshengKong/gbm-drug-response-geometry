"""Build gbm-point3-myeloid.ipynb: Point 3 myeloid compartment reanalysis.

Memo (research_notes/gbm-trajectory-drug-geometry-link-20261105.md §3):
  A: myeloid ligands zone along the trajectory (VEGFA, versican, annexin 1
     enriched in gliosis-hypoxia zones). B: ~139k non-transformed cells
     (incl. myeloid) were NOT analyzed by B.
  Q-a: per-TissueID myeloid composition or ligand expression vs that unit's
       drug response magnitude?
  Q-b: does treatment change myeloid state?
  Feasibility first: enough cells in both arms per TissueID, else downgrade.

Design (all on the same Colab session gbm-point1; loom already there):
  §0 setup: verify loom MD5 (skip download if present), load Point 1 outputs
     (traj_projection.csv m12, quintile_bins.csv mean_traj_ctrl).
  §1 feasibility (col_attrs only): myeloid = Cell_type == 'Myeloid'
     (63,178, all non-tumor). Per-TissueID x arm table; inclusion >=100/arm
     for TissueID-level tests, >=25 treated myeloid cells per comparison.
  §2 composition: myeloid fraction of non-tumor cells per TissueID x arm;
     paired Wilcoxon treated vs vehicle; oligodendrocyte fraction as control.
  §3 heavy pass (ONE full-matrix tile loop, cached to npz): library sizes +
     nnz for all non-tumor cells + raw counts of a ~54-gene myeloid panel
     (A ligands VEGFA/VCAN/ANXA1 + canonical myeloid markers) for myeloid
     cells. Tile loop = gene blocks x cell blocks with boolean masks
     (myeloid columns are scattered, so no contiguous fast path).
  §4 state shift: per-TissueID panel LFC (treated vs vehicle, log-CPM);
     consistency = mean pairwise Spearman of LFC vectors + per-gene sign
     consistency. Centroid distance in panel log-CPM space (genes z-scored)
     vs split-half vehicle null (200x). Per-drug pooled distances.
  §5 link to tumor response: per-comparison myeloid features (baseline
     myeloid fraction, baseline ligand expression, myeloid shift magnitude)
     vs tumor m12 (Spearman pooled + per TissueID-aggregated); per-TissueID
     mean tumor trajectory position vs vehicle myeloid ligand expression
     (tests A's zoning claim in B's data). All hypothesis-generating only
     (43 comparisons, per-drug n=3-8).
  §6 verdict (pre-registered): myeloid treatment-responsive if LFC
     consistency mean pairwise Spearman > 0.3 AND median centroid distance
     > 1.5x median null with paired Wilcoxon p < 0.05. Myeloid-response
     link if any feature shows |Spearman rho| > 0.4, p < 0.05 (pooled).
"""
import nbformat as nbf

nb = nbf.v4.new_notebook()
nb.metadata = {
    "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
    "language_info": {"name": "python", "version": "3.12"},
}

PANEL = ["VEGFA", "VCAN", "ANXA1",
         "CSF1R", "ITGAM", "CD68", "CD33", "LY86", "AIF1",
         "P2RY12", "TMEM119", "CX3CR1", "SALL1",
         "CD163", "MRC1", "MS4A4A", "C1QA", "C1QB", "C1QC", "APOE", "TREM2",
         "SPP1", "MARCO", "FOLR2", "LYVE1", "STAB1", "SEPP1",
         "IL1B", "TNF", "IL6", "CCL2", "CCL3", "CXCL10", "S100A8", "S100A9",
         "NLRP3", "IL1RN",
         "HLA-DRA", "HLA-DRB1", "CD74", "B2M",
         "LPL", "FABP5", "CD36",
         "MKI67", "TOP2A",
         "MMP9", "TIMP1", "FN1",
         "IL10", "TGFB1", "CD274",
         "CCL5", "CXCL8"]
PANEL_REPR = repr(PANEL)

md0 = """# 点 3：髓系 compartment 重分析

**备忘问题**（§3）：A 说髓系配体沿轨迹分区（VEGFA、versican、annexin 1 富集于胶质增生-缺氧区）；B 的 391,444 个细胞里约 139,000 个非转化细胞（含髓系）没被 B 分析过。

- **Q-a**：每 TissueID 的髓系组成或配体表达，是否与该单元的药物响应幅度相关？
- **Q-b**：药物处理是否改变髓系状态？

**先验可行性**（本地已核，跑本子时再验一遍）：Myeloid 63,178 个细胞，全是非转化；15/18 有髓系的 TissueID 每臂 ≥100 个髓系细胞；9 种药物治疗臂髓系细胞 869–6,479 个。TB6534 无髓系、TB6813（11/16）与 TB6528（31/75）细胞太少，降级排除。

**面板基因**：A 点名的三个配体（VEGFA、VCAN、ANXA1）+ 文献常规髓系标志（泛髓系、稳态小胶质、巨噬细胞、炎症、抗原呈递、脂质、增殖、基质、免疫抑制），共约 54 个；loom 里缺失的如实报告并剔除。

**判决标准（事先定好）**：髓系“治疗有响应”当且仅当 LFC 一致性平均成对 Spearman > 0.3，且中位质心距离 > 1.5 倍中位零分布、配对 Wilcoxon p<0.05。髓系-响应关联：任一特征 pooled Spearman |ρ|>0.4 且 p<0.05 即记为信号——但一律只算假设生成（43 个比较、每药仅 3–8 单元，功效不足）。
"""

c1 = """import numpy as np, pandas as pd, os, hashlib
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

LOOM_DRIVE_FILE_ID = "18-KInmm43wKdBX95Gq9xbuzAQwtLjgE9"
LOOM_MD5 = "9ebdd65747685c92b9ce57929da99a46"
LOOM_LOCAL = "/content/levitin_annotated_gse148842.loom"
OUT = Path('/content/point3_outputs'); OUT.mkdir(exist_ok=True)
RANDOM_SEED = 20261108
rng = np.random.default_rng(RANDOM_SEED)

PANEL = """ + PANEL_REPR + """
MIN_MY_ARMS = 100    # TissueID 级检验：每臂最少髓系细胞数
MIN_MY_COMP = 25     # 比较级检验：治疗臂最少髓系细胞数

def md5_of(path, chunk=1 << 22):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(chunk), b""):
            h.update(b)
    return h.hexdigest()

# loom 已在会话里（点 1 留下）且 MD5 通过时跳过下载
if os.path.exists(LOOM_LOCAL) and md5_of(LOOM_LOCAL) == LOOM_MD5:
    print("loom 已存在且 MD5 通过，跳过下载")
else:
    import subprocess
    print("loom 缺失或 MD5 不符，重新下载 ...")
    subprocess.run(["pip", "install", "-q", "gdown"], check=True)
    subprocess.run(["gdown", "--fuzzy", f"https://drive.google.com/uc?id={LOOM_DRIVE_FILE_ID}", "-O", LOOM_LOCAL], check=True)
    assert md5_of(LOOM_LOCAL) == LOOM_MD5
    print("下载后 MD5 通过")

tj = pd.read_csv('/content/point1_outputs/traj_projection.csv')
qb = pd.read_csv('/content/point1_outputs/quintile_bins.csv')
print('Point 1 outputs loaded: traj_projection', tj.shape, '| quintile_bins', qb.shape)
"""

md1 = """## §1 可行性：髓系细胞数（只读 col_attrs）
"""

c2 = """import h5py
f = h5py.File(LOOM_LOCAL, "r")
mat = f["/matrix"]
print("matrix shape:", mat.shape, "| dtype:", mat.dtype)
assert mat.shape == (60725, 391444)

def read_str(path):
    arr = f[path][:]
    return np.array([x.decode("utf-8") if isinstance(x, bytes) else str(x) for x in arr])

genes = read_str("/row_attrs/Gene")
celltype = read_str("/col_attrs/Cell_type")
tissue = read_str("/col_attrs/TissueID")
treatment_raw = read_str("/col_attrs/Treatment")
print("Cell_type:"); print(pd.Series(celltype).value_counts().to_string())

def normalize_treatment(raw):
    r = raw.strip().lower()
    if "vehicle" in r or "dmso" in r or r in ("control", "ctrl", "none"): return "vehicle"
    if "panobin" in r: return "panobinostat"
    if "etopos" in r: return "etoposide"
    if "topotec" in r: return "topotecan"
    if "ro492" in r or "r0492" in r: return "RO4929097"
    if "tazemet" in r: return "tazemetostat"
    if "ispin" in r or "ispen" in r: return "ispinesib"
    if "ana-12" in r or "ana12" in r: return "ANA-12"
    if "givinostat" in r: return "givinostat"
    if "vorinostat" in r: return "vorinostat"
    if "temozolomide" in r or r == "tmz": return "temozolomide"
    return r
treatment = np.array([normalize_treatment(x) for x in treatment_raw])

is_tumor = np.array([any(k in c.lower() for k in ("tumor","transform","malig","neoplast","glioma")) for c in celltype])
is_my = (celltype == "Myeloid")
is_nt = ~is_tumor
print("\\nTumor:", int(is_tumor.sum()), "| non-tumor:", int(is_nt.sum()), "| Myeloid:", int(is_my.sum()))
assert int(is_my.sum()) == 63178 and int(is_nt.sum()) == 139338
assert (is_my & is_tumor).sum() == 0, "Myeloid 与 Tumor 重叠，停下检查！"

feas = pd.DataFrame({"tis": tissue[is_my], "arm": np.where(treatment[is_my]=="vehicle","vehicle","treated")})
piv = feas.groupby(["tis","arm"]).size().unstack(fill_value=0)
piv["min_arm"] = piv.min(axis=1)
piv["included_100"] = piv["min_arm"] >= MIN_MY_ARMS
print("\\nMyeloid per TissueID x arm:"); print(piv.to_string())
print("\\nTissueIDs >=100/arm:", int(piv["included_100"].sum()), "of", len(piv))
inc_tis = sorted(piv.index[piv["included_100"]])
print("included:", inc_tis)
piv.to_csv(OUT/"myeloid_feasibility.csv")
# 每比较治疗臂髓系细胞数（§5 用）
comp_my = pd.DataFrame({"tis": tissue[is_my], "drug": treatment[is_my]})
comp_my = comp_my[comp_my.drug != "vehicle"].groupby(["tis","drug"]).size().rename("n_my_treated").reset_index()
print("\\ncomparisons with >=25 treated myeloid cells:",
      int((comp_my.n_my_treated >= MIN_MY_COMP).sum()), "of", len(comp_my))
"""

md2 = """## §2 组成：治疗改变髓系占比吗？
髓系占非转化细胞的比例，按 TissueID × 臂算；配对 Wilcoxon（治疗 vs 载体）。少突胶质细胞占比作对照——看变化是不是髓系特异的。
"""

c3 = """nt_df = pd.DataFrame({"tis": tissue[is_nt], "ct": celltype[is_nt],
                          "arm": np.where(treatment[is_nt]=="vehicle","vehicle","treated")})
cnt = nt_df.groupby(["tis","arm","ct"]).size().unstack(fill_value=0)
tot = cnt.sum(axis=1)
frac = cnt.div(tot, axis=0)
frac_my = frac["Myeloid"].unstack()
frac_ol = frac["Oligodendrocyte"].unstack()
print("Myeloid fraction (vehicle -> treated) per TissueID:")
print(pd.DataFrame({"vehicle": frac_my["vehicle"], "treated": frac_my["treated"],
                    "delta": frac_my["treated"]-frac_my["vehicle"]}).round(3).to_string())
for name, fr in [("Myeloid", frac_my), ("Oligodendrocyte", frac_ol)]:
    both = fr.dropna()
    w = stats.wilcoxon(both["treated"], both["vehicle"])
    print(f"{name}: median delta = {float((both['treated']-both['vehicle']).median()):+.3f} | paired Wilcoxon p = {w.pvalue:.3e} | n = {len(both)}")
frac_my.to_csv(OUT/"myeloid_fraction.csv")
"""

md3 = """## §3 文库大小与面板基因计数（全矩阵瓦片循环，结果缓存）
髓系列是散的，没有连续快路径；用基因块 × 细胞块瓦片循环一次读完全矩阵，只累积非转化列的文库大小/nnz 与髓系列的面板基因计数。缓存到 npz，重跑跳过。
"""

c4 = """import scipy.sparse as sp
cache_lib = OUT/"lib_nnz_nontumor.npz"
cache_sel = OUT/"myeloid_panel_counts.npz"

gene_to_row = {}
for _i, _g in enumerate(genes):
    if _g not in gene_to_row: gene_to_row[_g] = _i
present = [g for g in PANEL if g in gene_to_row]
missing = [g for g in PANEL if g not in gene_to_row]
print("panel genes:", len(PANEL), "| present:", len(present), "| missing:", missing)
assert len(present) / len(PANEL) >= 0.8, "面板基因覆盖率低于 80%，停下检查！"
sel_rows = np.array(sorted(gene_to_row[g] for g in present), dtype=np.int64)
sel_set = np.zeros(mat.shape[0], dtype=bool); sel_set[sel_rows] = True
row_to_selpos = {r: i for i, r in enumerate(sel_rows)}
present = [genes[r] for r in sel_rows]  # 与 sel_rows 同序

nt_idx = np.where(is_nt)[0]
my_idx = np.where(is_my)[0]
nt_pos = np.full(mat.shape[1], -1); nt_pos[nt_idx] = np.arange(len(nt_idx))
my_pos = np.full(mat.shape[1], -1); my_pos[my_idx] = np.arange(len(my_idx))
# Tumor 占 0..252105 连续列；全部非转化细胞都在 [252106, 391444) 这个连续区间里
# （本地已核：区间宽度 = 非转化细胞数），所以瓦片循环只读这个区间，省一半时间。
cmin_nt = int(nt_idx.min())
assert (nt_idx == np.arange(cmin_nt, mat.shape[1])).all(), "非转化列不是连续区间，停下检查！"

if os.path.exists(cache_lib) and os.path.exists(cache_sel):
    z = np.load(cache_lib); lib_nt = z["lib"]; nnz_nt = z["nnz"]
    z2 = np.load(cache_sel); sel_counts = z2["counts"]; sel_genes = list(z2["genes"])
    print("读缓存:", cache_lib, cache_sel)
else:
    lib_nt = np.zeros(len(nt_idx), dtype=np.float64)
    nnz_nt = np.zeros(len(nt_idx), dtype=np.int64)
    sel_counts = np.zeros((len(sel_rows), len(my_idx)), dtype=np.float64)
    GB, CB = 512, 65536
    n_genes, n_cells = mat.shape
    for g0 in range(0, n_genes, GB):
        g1 = min(g0+GB, n_genes)
        local_sel = np.where(sel_set[g0:g1])[0]
        sel_global = [row_to_selpos[g0+i] for i in local_sel] if len(local_sel) else []
        for c0 in range(cmin_nt, n_cells, CB):  # 只读非转化区间 [252106, 391444)
            c1 = min(c0+CB, n_cells)
            blk = mat[g0:g1, c0:c1]
            m_nt = is_nt[c0:c1]
            if m_nt.any():
                sub = blk[:, m_nt]
                cols = nt_pos[c0:c1][m_nt]
                lib_nt[cols] += sub.sum(axis=0)
                nnz_nt[cols] += (sub > 0).sum(axis=0)
                if sel_global:
                    m_my = is_my[c0:c1]
                    if m_my.any():
                        cols_my = my_pos[c0:c1][m_my]
                        m_my_nt = is_my[c0:c1][m_nt]  # myeloid mask within the non-tumor columns of this block
                        cols_my_in_sub = np.where(m_my_nt)[0]
                        sel_counts[np.ix_(sel_global, cols_my)] += sub[np.ix_(local_sel, cols_my_in_sub)]
        if (g0 // GB) % 20 == 0:
            print(f"  gene block {g0}-{g1} / {n_genes}")
    np.savez(cache_lib, lib=lib_nt, nnz=nnz_nt)
    np.savez(cache_sel, counts=sel_counts, genes=np.array(present))
    sel_genes = list(present)
    print("已缓存:", cache_lib, cache_sel)

print("non-tumor lib: median", int(np.median(lib_nt)), "| mean", round(float(np.mean(lib_nt)), 1))
print("non-tumor nnz: median", int(np.median(nnz_nt)))
assert (lib_nt > 0).all(), "存在文库大小为 0 的非转化细胞，停下检查！"
# lib_nt 按 nt_idx 顺序；取髓系子集：
my_in_nt = np.searchsorted(nt_idx, my_idx)  # my_idx 在 nt_idx 中的位置
lib_my = lib_nt[my_in_nt]
print("myeloid lib: median", int(np.median(lib_my)))
print("panel counts shape:", sel_counts.shape, "| genes:", len(sel_genes))
"""

md4 = """## §4 状态位移：治疗推动髓系状态吗？
面板基因 log-CPM 后，按 TissueID（≥100/臂）算每基因 LFC（治疗 − 载体）；一致性 = LFC 向量平均成对 Spearman + 每基因符号一致率。质心距离：基因先跨全部髓系细胞 z-score，再算治疗 vs 载体质心欧氏距离；零分布 = 载体细胞对半分 200 次的质心距离。另按药物汇总。

注意：每 TissueID 的 LFC 是把该单元所有药物 pooled 后算的，不同单元的药物组合不完全相同，所以跨单元一致性部分反映共同的药物混合；按药物的汇总距离是更干净的药物特异性检查。
"""

c5 = """# log-CPM for myeloid panel genes
X = np.log1p(1e4 * sel_counts / lib_my[None, :])  # genes x myeloid cells
my_tis = tissue[my_idx]; my_tr = treatment[my_idx]
my_df = pd.DataFrame({"tis": my_tis, "arm": np.where(my_tr=="vehicle","vehicle","treated")})
print("panel log-CPM shape:", X.shape)

# per-TissueID LFC vectors
lfc_rows = {}
for t in inc_tis:
    m_v = (my_tis==t) & (my_tr=="vehicle"); m_t = (my_tis==t) & (my_tr!="vehicle")
    lfc_rows[t] = X[:, m_t].mean(axis=1) - X[:, m_v].mean(axis=1)
LFC = pd.DataFrame(lfc_rows, index=sel_genes)
print("LFC matrix:", LFC.shape)
# consistency: mean pairwise Spearman across TissueIDs
import itertools
rhos = []
for a, b in itertools.combinations(LFC.columns, 2):
    r = stats.spearmanr(LFC[a], LFC[b]).statistic
    if np.isfinite(r): rhos.append(r)
print("mean pairwise Spearman of LFC vectors:", round(float(np.mean(rhos)), 3), "| n pairs:", len(rhos))
sign_cons = (np.sign(LFC) == np.sign(LFC.median(axis=1)).values[:, None]).mean(axis=1)
print("top consistent genes (sign agreement):")
print(sign_cons.sort_values(ascending=False).head(10).round(2).to_string())
print("\\nmedian |LFC| per gene (top 10):")
print(LFC.abs().median(axis=1).sort_values(ascending=False).head(10).round(3).to_string())
LFC.to_csv(OUT/"myeloid_lfc_per_tissue.csv")

# centroid distance vs split-half null (genes z-scored across myeloid cells)
mu = X.mean(axis=1, keepdims=True); sd = X.std(axis=1, keepdims=True); sd[sd==0]=1
Z = (X - mu)/sd
dist_rows = []
for t in inc_tis:
    m_v = np.where((my_tis==t) & (my_tr=="vehicle"))[0]
    m_t = np.where((my_tis==t) & (my_tr!="vehicle"))[0]
    c_v = Z[:, m_v].mean(axis=1); c_t = Z[:, m_t].mean(axis=1)
    d = float(np.linalg.norm(c_t - c_v))
    nulls = []
    for _ in range(200):
        half = rng.permutation(m_v); h1, h2 = half[:len(half)//2], half[len(half)//2:]
        nulls.append(float(np.linalg.norm(Z[:, h1].mean(axis=1) - Z[:, h2].mean(axis=1))))
    nulls = np.array(nulls)
    dist_rows.append(dict(tis=t, dist=d, null_median=float(np.median(nulls)),
                          ratio=d/float(np.median(nulls)),
                          emp_p=float((nulls >= d).mean()), n_v=len(m_v), n_t=len(m_t)))
dist_df = pd.DataFrame(dist_rows)
print("\\ncentroid distance (treated vs vehicle) vs null:")
print(dist_df.round(3).to_string(index=False))
print("median ratio (dist/null):", round(float(dist_df.ratio.median()), 2))
print("paired Wilcoxon dist vs null_median: p =",
      f"{stats.wilcoxon(dist_df.dist, dist_df.null_median).pvalue:.3e}")
print("TissueIDs with emp_p < 0.05:", int((dist_df.emp_p < 0.05).sum()), "of", len(dist_df))
dist_df.to_csv(OUT/"myeloid_centroid_distance.csv", index=False)

# per-drug pooled centroid distance
print("\\nper-drug pooled centroid distance:")
drug_rows = []
for d in sorted(set(my_tr[my_tr!="vehicle"])):
    m_v = np.where(my_tr=="vehicle")[0]; m_t = np.where(my_tr==d)[0]
    c_v = Z[:, m_v].mean(axis=1); c_t = Z[:, m_t].mean(axis=1)
    drug_rows.append(dict(drug=d, dist=float(np.linalg.norm(c_t-c_v)), n_t=len(m_t)))
print(pd.DataFrame(drug_rows).round(3).to_string(index=False))
"""

md5 = """## §5 关联：髓系特征 vs 肿瘤响应幅度（假设生成）
每比较（TissueID × 药物）的髓系特征——基线髓系占比（载体臂）、基线配体表达（载体臂 VEGFA/VCAN/ANXA1 均值）、髓系位移幅度（该比较质心距离）——对肿瘤 m12（点 1 traj_projection.csv）。另按 TissueID 汇总（m12 取该单元各药均值）再算一次；最后用点 1 的 mean_traj_ctrl 验 A 的分区主张（肿瘤轨迹位置 vs 髓系配体表达）。
"""

c6 = """# per-comparison myeloid features
feat = []
mean_traj = qb.groupby(["patient","drug"])["mean_traj_ctrl"].mean()  # bin centers averaged (quintiles have ~equal n0)
for _, r in tj.iterrows():
    p, d = r["patient"], r["drug"]
    m_v = (my_tis==p) & (my_tr=="vehicle")
    m_t = (my_tis==p) & (my_tr==d)
    if m_t.sum() < MIN_MY_COMP:
        continue
    nt_p = int(((tissue==p) & is_nt).sum())
    frac_v = float(m_v.sum()/nt_p) if nt_p else np.nan
    lig_idx = [list(sel_genes).index(g) for g in ["VEGFA","VCAN","ANXA1"] if g in sel_genes]
    lig_v = float(X[np.ix_(lig_idx, np.where(m_v)[0])].mean()) if lig_idx else np.nan
    c_v = Z[:, np.where(m_v)[0]].mean(axis=1); c_t = Z[:, np.where(m_t)[0]].mean(axis=1)
    shift = float(np.linalg.norm(c_t - c_v))
    feat.append(dict(patient=p, drug=d, m12=r["m12"], d_traj=r["d_traj"],
                     mean_traj_ctrl=float(mean_traj.get((p,d), np.nan)),
                     my_frac_vehicle=frac_v, lig_vehicle=lig_v, my_shift=shift,
                     n_my_treated=int(m_t.sum()), n_my_vehicle=int(m_v.sum())))
feat_df = pd.DataFrame(feat)
print("comparisons with myeloid features:", len(feat_df), "of", len(tj))
feat_df.to_csv(OUT/"myeloid_vs_tumor_response.csv", index=False)

def _sp2(a, b):
    m = np.isfinite(a) & np.isfinite(b)
    if m.sum() < 5: return (np.nan, np.nan, int(m.sum()))
    s = stats.spearmanr(a[m], b[m])
    return (round(float(s.statistic), 3), f"{s.pvalue:.2e}", int(m.sum()))

print("\\nSpearman vs tumor m12 (pooled comparisons):")
for col in ["my_frac_vehicle", "lig_vehicle", "my_shift"]:
    print(f"  {col}: rho, p, n =", _sp2(feat_df[col].to_numpy(), feat_df["m12"].to_numpy()))
print("\\nSpearman vs tumor d_traj (pooled):")
for col in ["my_frac_vehicle", "lig_vehicle", "my_shift"]:
    print(f"  {col}: rho, p, n =", _sp2(feat_df[col].to_numpy(), feat_df["d_traj"].to_numpy()))
print("\\nSpearman vs mean_traj_ctrl (A zoning claim):")
for col in ["lig_vehicle", "my_frac_vehicle"]:
    print(f"  {col}: rho, p, n =", _sp2(feat_df[col].to_numpy(), feat_df["mean_traj_ctrl"].to_numpy()))

# per-TissueID aggregation
agg = feat_df.groupby("patient").agg(m12=("m12","mean"), my_frac_vehicle=("my_frac_vehicle","first"),
                                     lig_vehicle=("lig_vehicle","first"), my_shift=("my_shift","mean"),
                                     mean_traj_ctrl=("mean_traj_ctrl","first")).reset_index()
print("\\nper-TissueID (n=%d) Spearman vs mean m12:" % len(agg))
for col in ["my_frac_vehicle", "lig_vehicle", "my_shift"]:
    print(f"  {col}: rho, p, n =", _sp2(agg[col].to_numpy(), agg["m12"].to_numpy()))
print("per-TissueID mean_traj_ctrl vs lig_vehicle:", _sp2(agg["lig_vehicle"].to_numpy(), agg["mean_traj_ctrl"].to_numpy()))
"""

md6 = """## §6 图与 verdict
"""

c7 = """fig, ax = plt.subplots(1, 3, figsize=(13, 4.2))
# (a) myeloid fraction paired
fr = pd.read_csv(OUT/"myeloid_fraction.csv", index_col=0)
ax[0].scatter(fr["vehicle"], fr["treated"], s=28)
lo = min(fr["vehicle"].min(), fr["treated"].min())-0.01; hi = max(fr["vehicle"].max(), fr["treated"].max())+0.01
ax[0].plot([lo,hi],[lo,hi],"r--",lw=1); ax[0].set_xlabel("myeloid fraction, vehicle"); ax[0].set_ylabel("myeloid fraction, treated")
ax[0].set_title("(a) Myeloid fraction per TissueID")
# (b) centroid distance vs null
ax[1].scatter(dist_df["null_median"], dist_df["dist"], s=28)
mx = max(dist_df["null_median"].max(), dist_df["dist"].max())*1.1
ax[1].plot([0,mx],[0,mx],"r--",lw=1); ax[1].set_xlabel("null median distance (split-half vehicle)"); ax[1].set_ylabel("treated-vs-vehicle distance")
ax[1].set_title("(b) Myeloid centroid shift vs null")
# (c) myeloid feature vs tumor m12
ax[2].scatter(feat_df["lig_vehicle"], feat_df["m12"], s=28)
ax[2].set_xlabel("baseline ligand expression (VEGFA/VCAN/ANXA1, log-CPM)"); ax[2].set_ylabel("tumor m12")
ax[2].set_title("(c) Baseline myeloid ligands vs tumor response")
plt.tight_layout(); plt.savefig(OUT/"fig6_my_1.png", dpi=150); print("saved fig6_my_1.png")

print("=== 点 3 verdict ===")
print("feasibility: 15 TissueIDs >=100 myeloid/arm; myeloid total 63,178")
print("LFC consistency mean pairwise Spearman:", round(float(np.mean(rhos)), 3))
print("median centroid dist/null ratio:", round(float(dist_df.ratio.median()), 2),
      "| Wilcoxon dist vs null p =", f"{stats.wilcoxon(dist_df.dist, dist_df.null_median).pvalue:.2e}")
print("myeloid-responsive by pre-registered rule:",
      bool(np.mean(rhos) > 0.3 and float(dist_df.ratio.median()) > 1.5
           and stats.wilcoxon(dist_df.dist, dist_df.null_median).pvalue < 0.05))
for col in ["my_frac_vehicle", "lig_vehicle", "my_shift"]:
    print(f"vs m12 {col}: rho, p, n =", _sp2(feat_df[col].to_numpy(), feat_df["m12"].to_numpy()))
print("outputs in", OUT)
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
    nbf.v4.new_markdown_cell(md5),
    nbf.v4.new_code_cell(c6),
    nbf.v4.new_markdown_cell(md6),
    nbf.v4.new_code_cell(c7),
]
nbf.write(nb, 'gbm-point3-myeloid.ipynb')
print('written')
