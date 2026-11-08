# Meylan MM2-缺氧 vs T 细胞空间排斥分析结果（2026-11-08）

**问题**：缺氧（MM2-high）肿瘤细胞周围的 T 细胞是不是更少（免疫排斥）？
**数据**：Meylan GSE296577 Xenium，8 位病人，治疗前/后配对。肿瘤细胞
1,591,862 个，T 细胞 172,386 个。空间距离在每个 orig.ident 切片内
用 cKDTree 计算（肿瘤细胞到最近 T 细胞的距离，µm）。

## 方法学关键点：4 基因 MM2 是弱代理，作者的 mes2 才是正解
- MM2 的 6 个非 VEGFA 基因中，**CAV1 与 LOX 不在 Meylan panel 上**，
  只能用 4 个基因（MET、ICAM1、CXCL10、CD44）算 MM2 分。
- 验证：4 基因 MM2 与作者的 mes2 分的相关只有 **rho=0.158**；
  按 tier2 分，4 基因 MM2 在 MES-like 1（0.331）反而高于
  MES-like 2（hypoxia，0.239），而作者 mes2 在 MES-like 2 是 0.715、
  远高于其他态（0.19–0.31）。**4 基因 MM2 不能干净地识别缺氧态。**
- 可能的原因：CXCL10 是 T 细胞趋化因子（招募 T 细胞），CD44 表达广泛；
  用含 CXCL10 的分数去检验"T 细胞排斥"，本身就被趋化信号污染，
  会把排斥效应抵消掉。
- 因此主检验改用两个更可靠的缺氧定义：**作者 mes2 分**（连续）与
  **tier2 的 MES-like 2（hypoxia）离散态**。4 基因 MM2 只作为阴性对照报告。

## 主结果（POST 治疗后样本，n=8 病人配对单位）

### A. 作者 mes2（中位二分，per-section）
- Pooled：mes2-high 中位距离 **60.7 µm** vs mes2-low 49.5 µm（均值 92.7 vs 65.5）；
  50 µm 内有 T 细胞的占比 0.408 vs 0.506。
- **每病人配对 Wilcoxon（n=8）：距离中位差 +10.4 µm，7/8 为正，p=0.0234**；
  占比中位差 −0.105，p=0.0234。
- 每细胞 Spearman（mes2 vs 距离）：rho=0.098，p≈0（效应小但方向一致）。
- 三分位敏感性：距离 Wilcoxon p=0.0234，中位差 **+16.9 µm**；
  占比 p=0.0234，中位差 −0.160。效应随缺氧程度加深而加大。

### B. tier2 MES-like 2（hypoxia）离散态 vs 其他肿瘤态（最干净的定义）
- Pooled：缺氧态中位距离 **90.7 µm** vs 其他态 51.2 µm（均值 129.8 vs 70.4）；
  50 µm 内有 T 细胞的占比 **0.281 vs 0.487**。
- **每病人配对 Wilcoxon（n=8）：距离中位差 +31.5 µm，8/8 全为正，
  p=0.0078（n=8 的最小值）**；占比中位差 −0.216，p=0.0156。
- 各病人距离差（µm）：P05 +13.9、P21 +113.5、P28 +41.9、P32 +5.1、
  P34 +21.2、P45 +55.0、P46 +0.04、P51 +76.2。P21、P51、P45 效应最大。

### 结论：T 细胞确实被排除在缺氧肿瘤区之外（治疗后）
- 三种缺氧定义（mes2 二分、三分位、tier2 离散态）方向完全一致，
  且离散态的效应最强（8/8 病人、中位 +31.5 µm）。
- 与 Q1 的发现互补：T 细胞治疗后大量涌入组织，但主要分布在间质/
  血管周围；进入肿瘤巢后，又进一步避开缺氧（MES-like 2）区域。
  这是"缺氧 niche 免疫冷"的空间证据。

## 阴性对照：4 基因 MM2（POST）无发现
- Pooled MW 距离：MM2-high 中位 54.2 µm vs MM2-low 54.7 µm，p=0.50。
- 每病人 Wilcoxon（n=8）：距离中位差 −0.9 µm，4/8 为正，p=0.46；
  占比中位差 −0.001，p=0.74。
- 每细胞 Spearman：rho=0.013。
- **解读**：不是"没有排斥"，而是 4 基因分数测不准缺氧（见方法学节）。
  以 mes2/tier2 为准，4 基因结果不作为阴性证据。

## 次要：PRE 治疗前样本
- 作者 mes2：每病人 Wilcoxon 距离中位差 +0.5 µm，p=0.64——**无一致效应**。
  治疗前 T 细胞极少（占比 0.06%–2.7%），距离主要由整体 T 密度决定，
  空间模式测不出来；排斥效应只在治疗后（T 细胞大量涌入后）才可测。
- 4 基因 MM2 在 PRE 反而显示 MM2-high 更近（Wilcoxon p=0.0234，
  中位差 −7.7 µm），方向与 POST 相反；鉴于 4 基因分数不可靠、
  且 PRE 的 T 密度极低，不作解释，只记录。

## 产物
- `meylan_analysis/mm2_genes_long.csv.gz`：4 基因每细胞表达（长格式，
  1,532,556 行非零记录；CAV1、LOX 不在 panel 上）。
- `meylan_analysis/mm2_tcell_tumor_cells.csv.gz`：每肿瘤细胞表
  （1,591,862 行：orig.ident、timepoint、tier2、x、y、mes2、MM2_4gene、
  dist_T、has_T_50）。
- `meylan_analysis/mm2_tcell_per_section_post.csv` /
  `mm2_tcell_per_section_pre.csv`：4 基因 MM2 的 per-section 表。
- `meylan_analysis/mm2_tcell_per_section_mes2_post.csv` /
  `mm2_tcell_per_section_tier2_post.csv`：mes2 与 tier2 的 per-section 表。
- 脚本：`extract_mm2_genes.R`（Colab 会话 gbm-meylan 内跑）、
  `meylan_mm2_tcell.py`、`meylan_mm2_tcell_followup.py`。

## 边界与注意事项
- 本检验是描述性空间分析，无生存/疗效标签，不建立"排斥导致耐药"的因果。
- 效应只在 POST 样本内可测；PRE 因 T 细胞太少而无法检验，不等于
  "治疗前无排斥"。
- 距离是二维切片内的欧氏距离，Xenium 切片厚度与三维结构未计入。
- 4 基因 MM2 的教训：panel 缺失关键基因时，基因集分数可能失真；
  有作者提供的模块分或离散态时，应优先用后者。
- 计算注：首版脚本在跑完主检验后因内存不足被系统杀死（exit 137），
  但每肿瘤细胞表已落盘；mes2/tier2 的跟进分析用落盘表重跑，内存占用小，
  正常完成。数字均来自落盘表，无缺失。
