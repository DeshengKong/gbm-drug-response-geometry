# Meylan GSE296577 空间分析设计（预注册）

## 数据
- Seurat 对象 `GSE296577_GBM_oHSV_seurat_Full_20260107.rds`（3.8 GB，2,401,946 细胞，480 基因，16 个 FOV，8 位病人，治疗前/后配对）。
- meta.data 实际 **23 列**（strings 扫描只看到 16 列，R 提取补全）：orig.ident、nCount/nFeature_RNA、ControlCodeword/ControlProbe、timepoint、**QC（全部 passed QC）**、seurat_clusters、**npc1、npc2、opc、ac、mes1、mes2（6 个模块分）**、HSV_total、batch（batch1-2 / batch3-4）、tier1、tier2、**x、y（坐标直接在 meta.data 里）**、tier1_T。
- 坐标在 meta.data 的 x、y 列（µm），centroids 提取文件只是交叉验证。
- tier1 共 12 类：Tumor 1,591,862、Macrophages/Monocytes 270,055、Vasculature 123,038、Neutrophils 82,424、Plasma cells 67,164、Stromal 36,676、Low count/Necrotic 35,340、T cells 172,386、B 10,709、NK 6,663、Dendritic 5,629。tier2 肿瘤 6 态：AC-like 489,779、NPC-like 1 378,173、OPC-like 265,100、MES-like 2（hypoxia）221,417、MES-like 1 130,127、NPC-like 2 107,263。tier1_T = tier1 + T 细胞 10 功能态。
- 目标基因 counts（22 个，稀疏矩阵 7.2M 非零）：VEGFA、VCAN、ANXA1、CD3D、CD3E、CD8A、CD4、GZMB、PDCD1、CD68、CD163、MS4A1、NCAM1、KLRD1、FOXP3、IFNG、CXCL10、CXCL9、CCL5、CD69、NR4A1、MKI67。缺失：HLA-DRA、TOP2A、HSV 基因（病毒载量改用 HSV_total 列）。
- HSV_total > 0：post 19,303 细胞，pre 1,075 细胞（pre 有少量背景），最大值 12.1。
- **Q3（TCR）已判死刑**：RDS 内无可用的克隆型映射（递归搜索全部 slot、commands、assay meta.data/misc、DimReduc misc、neighbors、FOV misc，只找到 RunPCA 命令里 `genes_wo_TCR_HSV` 一句）；GEO filelist 无 TCR/clonotype 文件；panel 无 TRBV/TRBJ 基因、无法从表达推断克隆型。README 的“TCR clonotype mapping”与可下载数据不符。RDS 序列化里 43 个 TCR-P28/P34-* 字符串（仅 P28、P34 两病人）无细胞映射，不可用。

## 红线（沿用 Points 1–3）
- 不做轨迹投影：panel 对 S10 的覆盖只有 4.8%，score_genes 不可行。
- 不把转录响应读成疗效；无生存标签，不做“预测疗效”。
- 病人是配对检验的单位（n=8，有前后配对者更少），不把细胞当独立样本。
- 空间结论须有置换/旋转零模型对照，不只报观测值。

## 问题与预注册检验

### Q1 溶瘤病毒后 T 细胞浸润是否增加（空间版）
- 指标：每病人每时间点，T 细胞占全部细胞比例；及肿瘤邻域内 T 细胞比例（距任一肿瘤细胞 ≤50 µm 的细胞中 T 细胞占比）。
- 主检验：8 位病人配对 Wilcoxon（post − pre），双侧。
- 零模型：病人内随机置换 pre/post 标签，复算差值分布。

### Q2 配体分区（Point 3 的空间延伸）
- 指标：肿瘤细胞 VEGFA/VCAN/ANXA1 表达（每基因二值化：表达 >0 vs 0，因 Xenium 计数稀疏）。
- 检验 A：配体阳性率在 tier2 各态间的分布，重点是 hypoxic MES-like 2 是否最高（Kruskal-Wallis + 事后两两）。
- 检验 B：配体阳性肿瘤细胞周围 50 µm 内髓系细胞密度 vs 配体阴性肿瘤细胞（病人内配对比较）。
- 零模型：随机打乱配体标签（保持阳性率），复算 B 的差值。

### Q3 TCR 克隆扩增与空间聚集
- **已判死刑，不做**：RDS 内无可用克隆型映射（见“数据”节）；GEO 无 TCR 文件；panel 无 TRBV/TRBJ。README 的说法与可下载数据不符。如实报告，不强行做。

### Q4 髓系细胞的空间重排
- 指标：Macrophages/Monocytes 在肿瘤邻域内的占比，前后配对比较；及髓系细胞与配体阳性肿瘤细胞的共定位（Q2B 的镜像）。
- 主检验：同 Q1 的配对 Wilcoxon。

## 计算路线
- Colab 高内存 + Python（pyreadr 读 RDS，scipy.sparse + cKDTree 做空间邻域查询）。不装 Seurat（省 10–20 分钟安装时间；只需表达矩阵子集 + meta.data + centroids）。
- 先提取轻量工作表：每细胞（patient、timepoint、FOV、tier1、tier2、tier1_T、x、y）+ 目标基因表达（VEGFA、VCAN、ANXA1、CD3D、CD3E、CD8A、CD4、GZMB、PDCD1、HSV 相关 + 内参），存为 parquet；后续检验只读此表，不再碰 RDS。
- 细胞数 ~240 万，cKDTree 邻域查询在 Colab 高内存可行；按 FOV 分块做，避免一次性全量距离矩阵。

## 产出
- 每问题一张主图 + 一张表，数字与源文件逐项核对；失败的检验如实报告（Point 2 v1/v2 的先例）。
- 结论只写数据支持的部分；“疗效预测”缺标签，明确写为不可做。
