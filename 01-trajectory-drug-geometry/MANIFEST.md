# GBM 轨迹 × 药物响应重分析（点 1）— 数据清单 MANIFEST

> 生成时间：2026-11-07（UTC）；下载执行者：Lotus 子代理
> 用户 2026-11-07 22:37 CST 拍板“先拉点 1 的数据”。
> 原则：只下公开的处理后数据；不注册、不联系作者、不下 EGA 原始 reads。
> 本 VM：磁盘 150 GB（起始可用 93 GB），内存 7 GB，4 核，无 GPU。
> 下载一律 curl（走代理环境变量）；并行分片（range）+ 断点续传。

---

## 1. A 侧 — de Jong 等 Nature 接受稿（DOI 10.1038/s41586-026-11234-8）

12 例 IDH 野生型原发 GBM，单细胞 + 空间多组学；1,025,329 个核。

### 1.1 核心文件：处理后 snRNA 对象（E-MTAB-17183）

| 文件 | 来源 URL | 大小 | 状态 |
|---|---|---|---|
| `GBM_space_snRNA.h5ad.gz` | https://www.ebi.ac.uk/biostudies/files/E-MTAB-17183/GBM_space_snRNA.h5ad.gz（302 跳转至 https://ftp.ebi.ac.uk/pub/databases/biostudies/E-MTAB-/183/E-MTAB-17183/Files/GBM_space_snRNA.h5ad.gz） | 11,267,029,458 B（10.5 GB） | ✅ 已完成并核验（见下文） |
| `GBM_space_snRNA.h5ad`（解压后） | 由上解压（`gunzip -k`，保留 .gz） | 62,410,688,915 B（58.1 GiB） | ✅ HDF5；X 为 CSR 稀疏矩阵 (1,025,329 细胞 × 36,601 基因)，nnz 2,571,309,456，float32 整数值计数；obs 20 列（见 §3），**无轨迹注释列**；obsm 仅 X_umap；uns 仅 hvg |
| `GBM_space_snRNA_README.md` | https://www.ebi.ac.uk/biostudies/files/E-MTAB-17183/GBM_space_snRNA_README.md | 2,001 B | ✅ |
| `GBM_space_ATAC_filtered_peaks_README.md` | 同上目录 | 1,646 B | ✅（仅说明文件，ATAC 对象本身未下，点 1 不需要） |
| `README_cellranger_arc.md` | 同上目录 | 1,196 B | ✅ |
| `README_metacell.md` | 同上目录 | 2,691 B | ✅ |
| `E-MTAB-17183.idf.txt` | 同上目录 | 9,491 B | ✅ |
| `E-MTAB-17183.sdrf.txt` | 同上目录 | 312,662 B | ✅ |

SDRF 核实：12 个 donor（AT3, AT4, AT5, AT6, AT7, AT9, AT10–AT15），155 个 multiome 反应（RNA+ATAC 各 155），310 行。

**⚠️ snRNA 下载事故与处理（2026-11-07 17:00–18:10 UTC）：**
- 首轮 4 分片并行下载（14:46–16:56，中间为让带宽给 loom 暂停一次），合并后字节数正好 11,267,029,458，但 `gzip -t` 报错 `invalid compressed data--format violated`。
- 根因定位：逐块解压（Python zlib，`gz_locate.py`）在压缩偏移 ~4,155,506,688 处失败，即 **part1 内偏移 1,338,749,323——恰为第一轮被 kill 时 part1 的断点续传接缝**。`curl -C -`（断点续传）与 `-r`（range）联用时续传偏移错位，导致接缝处数据损坏；part0 的接缝经解压验证是干净的（解压一路通过 part0 全程 0–2.8 GB）。
- 附带发现：EBI BioStudies URL 返回 302 到 ftp.ebi.ac.uk；经实测 **ftp.ebi.ac.uk 直连 URL 支持 range（206 + Content-Range）**，Sanger 镜像亦支持；分片头字节检查证实各 part 数据本身是对齐的（part0 以 gzip 魔数 `1f8b` 开头，part1–3 以中段压缩字节开头，非文件头重复）。
- 处理：保留已验证干净的 part0；删除 part1–3，2026-11-07 17:12 起用**无 `-C -` 的全新单次 range 请求**从直连 ftp URL 重下 part1–3（脚本 `redownload_parts123.sh`）。重下中 part1、part3 在 curl 退出码 0、无报错的情况下静默短缺 ~5.27 MB（疑似连接提前关闭而 curl 未察觉），已用精确 range 请求补齐两条尾巴（part1 尾 5,269,690 B、part3 尾 5,269,689 B）并核对四片字节数精确到位；重新合并后 `gzip -t` **全量通过（GZIP INTEGRITY OK）**。**教训：range 分片下载中途不要 kill 后用 `-C -` 续传；要续传应删掉该片整片重下；且每次下载后必须逐片核对字节数，curl 退出 0 不等于文件完整。**

**核验结果（2026-11-07 18:10–18:45 UTC，`inspect_h5ad.py` + `inspect_obs.py` + `inspect_x.py`，只读元数据）：**
- `file` = Hierarchical Data Format (version 5)；解压后 62,410,688,915 B。
- X：CSR 稀疏矩阵，shape **(1,025,329 细胞 × 36,601 基因)**，nnz **2,571,309,456**，float32，值为整数计数（首 10 个值 2,2,1,1,3,1,2,1,2,4；min 1.0、max 3434.0）。`raw/X` 与 X 同形同 nnz（同为 2,571,309,456），raw/var 另有 `feature_types`、`mt` 两列。
- layers：空。obsm：仅 `X_umap` (1,025,329, 2)。obsp：`connectivities`、`distances`。uns：仅 `hvg`（高变基因）。
- obs 20 列（见 §3）：cell_id、donor_id（12 例 AT3–AT15）、site_id、sample（155 反应）、n_genes_by_counts、total_counts、mt_frac、doublet_scores、phase、CNV_signal_mean、cnv_corr、neftel、celltypist、scPoli、ontology keywords、top_markers、cell_status（Malignant/TME）、annotation_coarse、annotation_granular、index。**其中没有任何轨迹位置 / pseudotime / latent time / 元模块评分列（§3 最终确认）。**
- var 3 列：ENSEMBL、SYMBOL、_index。
- ⚠️ 待分析时注意：obs 的 `n_genes_by_counts`（均值 6,693，min 4,884）与 X 的实际每细胞非零数（均值 2,507；文件首部细胞 ~10k、尾部细胞 ~800）对不上——说明 obs 的 QC 指标是在与 X 不同的（更完整的）计数表上算的，或 X 经过滤后未重算 QC。分析时以 X 的实际 nnz 为准，不要直接用 obs 的 n_genes_by_counts 做细胞过滤阈值。

### 1.2 补充表格（gbmspace.org，全部已下，字节数与服务端 Content-Length 一致）

| 文件 | 大小 | 内容 |
|---|---|---|
| `supplementary/EDFig1c1.csv` | 41,985,554 B | ATAC QC：log10 nFrags vs TSSEnrichment（图源数据） |
| `supplementary/EDFig1c2.csv` | 24,003,694 B | 同上 |
| `supplementary/EDFig1c3.csv` | 35,710,626 B | ATAC QC：log10 nFrags vs PromoterRatio |
| `supplementary/EDFig5b1_cellClass_L0_L1.csv` | 383 B | 克隆 × 细胞类别计数 |
| `supplementary/EDFig5b2_clone_bySite.csv` | 379 B | 克隆 × 取样位点计数 |
| `supplementary/EDFig5c_sankey_AI_vs_multiome.csv` | 144 B | AI vs multiome 克隆对照 |
| `supplementary/EDFig5h2_new.csv` | 3,863,262 B | 空间聚类表 |
| `supplementary/Fig4d_BAF_L0_L1.csv` | 770,070 B | BAF（等位基因分数）表 |

注：前三个大 CSV 首次下载曾被截断（字节数不足），已用 `curl -C -` 续传至与 Content-Length 完全一致。

### 1.3 Metacell 对象（gbmspace.org）

| 文件 | 大小 | 核验 |
|---|---|---|
| `metacell/rna_metacell_scdori_11_03_24.h5ad` | 216,042,856 B（与 Content-Length 一致） | ✅ `file` = HDF5；h5py 读到 X shape (16761, 3192) float32；obs index 为 `<cluster>\|<annotation>\|<batch>\|<donor>` 格式，例 `9\|Hypoxic\|cellranger-arc201_count_6404016915495b342ccbe43e4ad068f8\|AT12`，与 README 描述一致 |

`genescore_metacell_scdori_11_03_24.h5ad`（ATAC gene-activity，1,673,587,956 B）未下：点 1 不需要。

### 1.4 代码（已克隆）

| 仓库 | 本地路径 | 大小 | 说明 |
|---|---|---|---|
| github.com/PMBio/spaceTree | `code/a-side-spaceTree/` | 538 MB | 轨迹推断方法（spaceTree） |
| github.com/BayraktarLab/GBM_space_paper | `code/a-side-GBM_space_paper/` | 461 MB | 论文全部分析代码 + 图源数据（`data/` 255 MB，`metadata/` 904 KB，含 `metadata/2.TM/factor_id_map.tsv`——时间元模块→细胞簇→环境映射表） |

### 1.6 Nature 论文附表（Supplementary Tables 1–16，全部已下）

来源：https://www.nature.com/articles/s41586-026-11234-8 的 Supplementary information（static-content.springer.com），15 个 xlsx，字节数与服务端 Content-Length 全部一致。

| 文件 | 大小 | 内容 |
|---|---|---|
| `41586_2026_11234_MOESM5_ESM.xlsx` | 335,069 B | **Supplementary Tables S1–S16**；其中 **Table S10 = "Cell2fate recurrent module annotation"（12 个时间元模块的基因列表）** |
| `41586_2026_11234_MOESM7_ESM.xlsx` | 62,000 B | Fig1 图源数据（含 Neftel 各状态评分） |
| `41586_2026_11234_MOESM8_ESM.xlsx` | 70,000 B | Fig2 图源数据（空间 niche 因子） |
| `41586_2026_11234_MOESM9_ESM.xlsx` | 25,000 B | Fig3 图源数据 |
| `41586_2026_11234_MOESM10_ESM.xlsx` | 6,774,921 B | Fig4 图源数据（克隆组成） |
| `41586_2026_11234_MOESM11_ESM.xlsx` | 3,537,352 B | Fig5 图源数据；**Fig5a2 sheet 含 6,520 个细胞的每细胞 latent time（`obs_names, Patient_Clone, time, cell_class_no_prolif`）** |
| `41586_2026_11234_MOESM12_ESM.xlsx` | 5,212 B | Fig6 图源数据 |
| `41586_2026_11234_MOESM13_ESM.xlsx` | 34,750 B | EDFig3 图源数据 |
| `41586_2026_11234_MOESM14_ESM.xlsx` | 21,422 B | EDFig4 图源数据 |
| `41586_2026_11234_MOESM15_ESM.xlsx` | 6,778,588 B | EDFig5 图源数据（克隆系统发育） |
| `41586_2026_11234_MOESM16_ESM.xlsx` | 3,226,967 B | EDFig6 图源数据；EDFig6c 为每 donor × 12 元模块的细胞数 |
| `41586_2026_11234_MOESM17_ESM.xlsx` | 6,282,635 B | EDFig7 图源数据（每克隆 latent time 分布） |
| `41586_2026_11234_MOESM18_ESM.xlsx` | 5,900,820 B | EDFig8 图源数据（元模块 marker 表达分数） |
| `41586_2026_11234_MOESM19_ESM.xlsx` | 35,555 B | EDFig9 图源数据 |
| `41586_2026_11234_MOESM20_ESM.xlsx` | 46,723 B | EDFig10 图源数据 |

**Table S10 提取结果（已另存为 `a-side/trajectory_metamodules_s10.csv` 与 `.json`）：**

| MetaModule | annotation | 基因数 |
|---|---|---|
| 0 | Cycling | 193 |
| 1 | NPC-neuronal-like 1 | 66 |
| 2 | Gliosis-hypoxia | 83 |
| 3 | AC-gliosis-like 2 | 86 |
| 4 | NPC-OPC-like 2 | 77 |
| 5 | AC-gliosis-like 1 | 68 |
| 6 | NPC-neuronal-like 2 | 77 |
| 7 | NPC-OPC-like 1 | 138 |
| 8 | AC-NPC-OPC-like 1 | 88 |
| 9 | AC-progenitor-like | 137 |
| 10 | NPC-neuronal-like 3 | 148 |
| 11 | OPC-like | 207 |

共 12 个模块、1,368 个基因条目（去重后 1,200 个基因），每模块另有富集通路（GO）一并存于 CSV/JSON。

### 1.7 未下载的 A 侧文件及原因

| 文件 | 大小 | 原因 |
|---|---|---|
| `GBM_space_ATAC_filtered_peaks.h5ad.gz`（E-MTAB-17183） | 22,413,672,480 B（20.9 GB） | ATAC 对象，点 1（RNA 轨迹 × 药物响应）不需要 |
| 约 100 个 `cellranger-arc201_count_*.gz`（E-MTAB-17183） | 每个 4–16 GB，总计约 1 TB | 每反应的原始 Cell Ranger ARC 输出；处理后 h5ad 已含注释，不需要 |
| E-MTAB-17164 全部文件（空间 Visium + Space Ranger 输出） | 总计约 840 GB | 空间数据；B 侧无空间坐标，点 1 不需要 |
| `GBM_LEAP_annotations_v2.18_12_24.h5ad`（gbmspace.org） | 62,974,906,629 B（58.7 GB） | 完整注释大对象；先看 snRNA h5ad 是否够用，不够再议 |
| `EDFig9b_atlas_GBO_integrated.h5ad`（补充） | 40,815,560,940 B（38 GB） | 类器官（GBO）整合图谱，点 1 不需要 |
| `Fig4a_AT10_clone_prediction.h5ad` / 两份 ST_map h5ad | 2.8 GB / 236 MB / 580 MB | 克隆预测与空间图，点 1 不需要 |
| EGA 原始 reads（EGAD00001015526/15527/16155/16126） | — | **点 1 不需要，绝对未碰** |

---

## 2. B 侧 — Lavrenteva & Kim 预印本（DOI 10.64898/2026.11.03.763852）

重分析 Zhao 2021 切片培养药物扰动数据（GSE148842；391,444 细胞，252,106 转化胶质瘤细胞进入分析，43 个治疗-载体比较，19 个 TissueID，9 种药物）。

### 2.1 首选：作者整合 loom（Levitin 等 2023，bioRxiv 10.1101/2023.12.05.570193）

| 文件 | 来源 | 大小 | 状态 |
|---|---|---|---|
| `levitin_annotated_gse148842.loom` | Google Drive 公开分享链接（file id `18-KInmm43wKdBX95Gq9xbuzAQwtLjgE9`，scArchon 论文引用的“Levitin 注释版”同一链接；Last-Modified 2023-12-05，与预印本日期一致；无登录，仅 confirm-token 流程） | 7,416,810,607 B（6.9 GB，服务端 Content-Length；8 分片并行下载后合并，字节数精确一致） | ✅ 已完成并核验 |

**核验结果（`inspect_loom.py` + `loom_detail.py`，只读元数据）：**
- 文件头为 HDF5 魔数；MD5：`9ebdd65747685c92b9ce57929da99a46`。
- `/matrix`：shape (60725 基因 × 391444 细胞)，int64，gzip 压缩，chunks (64,64)，存储占用 2.1 GB。**细胞数 391,444 与预印本完全一致。**
- `/row_attrs`：`Accession`（ENSG…，例 ENSG00000000003.14）、`Gene`（例 TSPAN6、TNMD、DPM1）——B 代码的基因字段检测可命中。
- `/col_attrs`：
  - `Cell_type`（8 类）：Tumor **252,106**（与预印本“252,106 转化胶质瘤细胞”完全一致）、Myeloid 63,178、Oligodendrocyte 39,957、Unspecific 31,577、TCell 2,058、Pericyte 1,312、Endothelial 1,142、Neuron 114。
  - `TissueID`（19 个）：TB5884（43,528）、TB5886、TB5944、TB5966、TB5974、TB5980、TB6140、TB6181、TB6224、TB6328、TB6393、TB6458、TB6488、TB6505、TB6528、TB6534、TB6545、TB6813、TB6814——19 个生物学单元，与预印本一致。
  - `Treatment`（10 个）：DMSO（载体，191,230）、Etoposide 36,068、Ispenisib 32,741、Topotecan 31,120、Panobinostat 21,646、Ana12 19,108、RO4929097 18,606、RSL3 18,554、Givinostat 14,160、RSL3_Ferrostatin 8,211——9 种药物 + DMSO 载体。
  - `Sequencing_Name`（108 个，形如 PW034-N703）：B 代码的患者字段检测（`pw\d+`）可命中。
- `/layers` 为空；`col_graphs`/`row_graphs` 为空。
- 与 B 代码 `inspectAndDetectLoom` 的检测逻辑逐项对照：基因字段 ✅、患者字段 ✅（Sequencing_Name）、治疗字段 ✅（Treatment，含 panobin/etopos/vehicle/dmso/topotec/ispin/ana-12/ro492/givinostat；仅 tazemet 未出现，不影响）、细胞类型字段 ✅（Cell_type，含 Tumor/Myeloid/Oligodendrocyte）、**必需的 `/col_attrs/TissueID` ✅**、library-size 字段未见（代码会退回用矩阵自行计算总 counts，不阻塞）。
- 各 TissueID 的药物臂：每个 TissueID 都有 DMSO + 1–6 种药物（例 TB5974 有 7 个臂，TB5980/TB6140/TB6224/TB6488/TB6505/TB6528/TB6534/TB6813 只有 DMSO + 1 种药）。

B 代码需要的 loom 字段（`Glioma_drug_geometry_analysis.m` 第 362–480 行）：
- `/matrix`（基因 × 细胞或细胞 × 基因）
- `/row_attrs/` 含基因名字段（须命中 marker 基因，如 CST3/GFAP/OLIG2/SOX4/VIM 等）
- `/col_attrs/` 含：患者字段（值形如 `PW030`）、治疗字段（值含 panobin/etopos/vehicle/dmso/topotec/ispin/ana-12/ro492/givinostat/tazemet）、细胞类型字段（值含 transform/myeloid/oligodend/malig），以及 **必需的 `/col_attrs/TissueID`**；可选 library-size 字段。

### 2.2 备选：GEO 原始文件（已下，作后备）

| 文件 | 来源 URL | 大小 | 状态 |
|---|---|---|---|
| `GSE148842_RAW.tar` | https://ftp.ncbi.nlm.nih.gov/geo/series/GSE148nnn/GSE148842/suppl/GSE148842_RAW.tar | 358,000,640 B（342 MB，与 filelist.txt 一致） | ✅ 已完成并核验 |

内含 43 个样本的 `.cts.txt.gz` 计数表（GSM4483741–GSM4483773，2020 年 33 个；GSM5073418–GSM5073427，2021 年 10 个），与 filelist.txt 逐项核对无缺失；`tar -tf` 列出 43 个条目；抽查 GSM4483741_PW029-701：60,726 行（基因）× 1,698 列（gid + gene + 1,696 个细胞 barcode），首行 `gid gene CAGAGGAAAATG …`，首基因 ENSG00000000003.14 TSPAN6。MD5：`fcdd65f5f8e99b6c52a1c2675a73109a`。

注意：RAW.tar 只有计数表，无细胞类型注释、无 TissueID；若 loom 不可用，需从 GEO 样本元数据 + 原文方法重建整合对象，工作量大。下载中曾在 98.5% 处断流，用 `curl -C -` 三次续传补齐。

### 2.3 代码与作者结果（已克隆）

| 仓库/文件 | 本地路径 | 说明 |
|---|---|---|
| github.com/bigbio-el/Glioma-drug-response-geometry | `code/b-side-drug-geometry/` | `Glioma_drug_geometry_analysis.m`（主分析）、`Inspect_glioma_loom.m`（loom 检查）、`results/cross_drug_unit_summary.csv`（11 个 TissueID 单元汇总）、`results/cross_drug_global_inference.csv`（全局推断） |

---

## 3. 关键问题：E-MTAB-17183 是否直接含有每细胞的轨迹位置注释？

**结论：h5ad 文件本身不含；但轨迹轴可从论文附表 S10 的 12 个元模块基因列表重建——点 1 的关键缺口已补上。**（2026-11-07 18:45 UTC 以实际文件最终确认。）

证据链：

1. **实测 obs 列名（`inspect_obs.py`，`GBM_space_snRNA.h5ad`，1,025,329 行）：共 20 列**——`cell_id, donor_id, site_id, sample`、QC 指标（`n_genes_by_counts, total_counts, mt_frac`）、`doublet_scores, phase, CNV_signal_mean, cnv_corr, neftel`（Neftel 2019 评分）、`celltypist, scPoli`、ontology keywords、`top_markers, cell_status, annotation_coarse, annotation_granular`，另有 `index`（与 cell_id 相同）。**其中没有任何轨迹位置 / pseudotime / latent time / 时间元模块评分列。**obsm 仅 `X_umap`（1,025,329 × 2）；uns 仅 `hvg`；obsp 为 `connectivities`、`distances`。`raw` 组与 X 同形同 nnz，无额外注释。
   - 各列示例值：`annotation_coarse` 12 类（AC-gliosis-like、AC-progenitor-like、Astrocytes、Gliosis-like、Hypoxic、Lymphocytes、Myeloid、NPC-neuronal-like、Neurons (Exc/Inh)、OPC-NPC-like、OPC-like）；`cell_status` = Malignant / TME；`donor_id` = AT3, AT4, AT5, AT6, AT7, AT9, AT10–AT15；`neftel` = AC、MES、NPC1、NPC2、OPC、Cycling 等；`phase` = G1/G2M/S。
2. 论文仓库的轨迹分析代码（`code/5.CELLULAR_TRAJECTORIES/Analysis/TemporalAnalysis/5.TemporalAnalysis.2.meta_module_clustering.ipynb`）显示：12 个时间元模块的基因列表（`factor_gene_dict`）是由每克隆的 cell2fate 输出（`../data/trajectories/` 下的 `*_res_tab*.csv`，**未公开**）聚类得到；每细胞的模块评分是分析时用 `sc.tl.score_genes(adata, factor_gene_dict[m], score_name=str(m)+'_score')` **现算的**，不是文件里存好的。
3. **但是，Nature 论文附表 Table S10（MOESM5）公开了 12 个元模块的完整基因列表**（见 §1.6；已提取为 `a-side/trajectory_metamodules_s10.csv/.json`，共 1,200 个去重基因）。点 1 的做法与论文一致：用这 12 个基因列表给 B 侧每个细胞打分（scanpy `score_genes`），得到 12 维模块评分，再投影到轨迹轴上。**不需要重跑 cell2fate。**
4. 另有验证参照：MOESM11 Fig5a2 提供 6,520 个图谱细胞的每细胞 latent time（cell2fate 潜时间），可用于校验“模块评分→轨迹位置”的映射是否与原文一致；MOESM16 EDFig6c 提供每 donor × 12 元模块的细胞数。
5. 公开仓库里只有 `metadata/2.TM/factor_id_map.tsv`（元模块→细胞簇→环境映射，如 `Dev-like (AC)`），没有元模块基因列表；`factor_genes_latest.pkl` / `all_top_genes_latest.csv` 等中间文件不在公开仓库中——但附表 S10 已等价覆盖。

---

## 4. 点 3 可行性核查（B 侧非转化细胞，含髓系）

**结论：可行。**（`inspect_loom.py` 实测，细胞类型字段以 B 代码同款逻辑判定：值含 transform/malig/neoplast/tumor/glioma 者为转化细胞。）

- 非转化细胞总数 **139,338**（与备忘的“约 139,000”一致）：载体臂（DMSO）70,261，治疗臂 69,077。
- 髓系（Myeloid）总数 **63,178**：载体臂 32,861，治疗臂 30,317。
- 每 TissueID × 每臂的髓系细胞数：

| TissueID | 载体臂髓系 | 治疗臂髓系 | 判定 |
|---|---|---|---|
| TB5884 | 4,028 | 7,185 | ✅ |
| TB5886 | 149 | 1,752 | ✅ |
| TB5944 | 2,247 | 2,674 | ✅ |
| TB5966 | 577 | 4,404 | ✅ |
| TB5974 | 2,037 | 3,707 | ✅ |
| TB5980 | 4,233 | 390 | ✅ |
| TB6140 | 621 | 122 | ✅ |
| TB6181 | 76 | 62 | ⚠️ 偏少（每臂 <100） |
| TB6224 | 927 | 107 | ✅ |
| TB6328 | 805 | 1,277 | ✅ |
| TB6393 | 2,027 | 5,337 | ✅ |
| TB6458 | 179 | 178 | ✅ |
| TB6488 | 7,782 | 1,024 | ✅ |
| TB6505 | 1,191 | 286 | ✅ |
| TB6528 | 31 | 75 | ⚠️ 偏少 |
| TB6534 | 0 | 0 | ❌ 无髓系细胞 |
| TB6545 | 5,454 | 1,244 | ✅ |
| TB6813 | 16 | 11 | ❌ 太少 |
| TB6814 | 481 | 482 | ✅ |

- **15/19 个 TissueID 每臂髓系 ≥100 个**，可支撑“每单元髓系组成/配体表达 vs 药物响应幅度”的关联分析；TB6181、TB6528、TB6813 三个单元偏少（建议降权或剔除），**TB6534 完全没有髓系细胞，必须剔除**。
- 其他非转化组分（供参考）：Oligodendrocyte 39,957（载体 23,744 / 治疗 16,213）、Unspecific 31,577（未分类；按 B 代码逻辑计入非转化）、TCell 2,058、Pericyte 1,312、Endothelial 1,142、Neuron 114。
- 注意：各 TissueID 的药物臂不同（每单元只有 DMSO + 1–6 种药），“某药的髓系响应”只能在有该药臂的单元内比较，不能跨单元硬拼。

---

## 5. 磁盘占用

| 项目 | 占用 |
|---|---|
| 起始可用 | 93 GB |
| A 侧已落盘（README/idf/sdrf + 补充 CSV + metacell + 代码仓库 + 附表） | 约 1.5 GB |
| A 侧 snRNA（.gz 压缩 + .h5ad 解压后，两份并存） | .gz 11,267,029,458 B（10.5 GB）+ .h5ad 62,410,688,915 B（58.1 GiB）= 约 69 GB |
| B 侧 loom（已完成并核验） | 7,416,810,607 B（6.9 GB） |
| B 侧 GEO RAW.tar（已完成并核验） | 358,000,640 B（0.34 GB） |
| gbm-point1 目录合计 | 约 78 GB |
| 当前磁盘 | 150 GB 总盘，2026-11-07 18:45 时可用约 69 GB |
| 说明 | .gz 与 .h5ad 并存是故意的：.gz 为下载存档（防 .h5ad 损坏可重解），.h5ad 为分析工作文件。磁盘吃紧时可删 .gz（11 GB），.h5ad 不可删。分析本身（scanpy 读 59 GB 稀疏矩阵）建议上 Colab（用户 AI Ultra 权益）或加内存，本 VM 仅 7 GB 内存直接全量加载会吃力。 |

---

## 6. 未下到的文件及原因

A 侧未下文件见 §1.7（ATAC 对象、约 100 个 cellranger-arc 原始输出、E-MTAB-17164 空间数据、LEAP 大对象、GBO 图谱等，均为点 1 不需要；EGA 原始 reads 绝对未碰）。

B 侧：作者整合 loom 与 GEO RAW.tar 均已下完并核验，无缺失。

**全部下载已完成（2026-11-07 18:45 UTC）。**snRNA 经事故重下后 `gzip -t` 全量通过并已解压核验；B 侧两文件均已核验；代码三仓库已克隆。无遗留缺口。
