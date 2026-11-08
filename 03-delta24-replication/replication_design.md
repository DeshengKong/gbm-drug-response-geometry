# GSE330794 复制分析设计（预注册）

## 数据
- MD Anderson Delta-24-RGD 队列（Fueyo），GEO SuperSeries
  "Clinical response to glioma virotherapy is linked to immune organization
  and T cell diversity"。本地已下齐 3.4 GB（61 文件）。
- TRB：16 个 MiXCR 文件（cloneId、cloneCount、cloneFraction、V/D/J、
  nSeqCDR3、aaSeqCDR3），每文件 30–353 个克隆。配对 pre/post **5 位病人**
  （P15、P21、P22、P27、P42）——P17、P20 只有 post、无 pre TRB，
  之前“7 位配对”的说法是错的；另有单时间点
  （P13-post、P17-post、P20-post、P31-pre、P39-post、P41-post）。
- Xenium：15 个人源 section（h5 counts + cells parquet + transcripts），
  panel 541 特征。配对 pre/post 仅 2 位病人（P17、P20）——Xenium 的
  pre/post 配对检验 n=2，不做 p 值，只做描述 + 横断面。
- panel 关键基因：VEGFA、VCAN、ANXA1、CD3D、CD3E、CD8A、CD4、CD68、
  CD163、MS4A1、CD79A、OLIG2、SOX2、EGFR、HIF1A、GZMB、PDCD1、FOXP3、
  NCAM1、KLRD1、CXCL10、CXCL9、CCL5、MKI67 均在；缺 PECAM1、VWF、
  GFAP、CD69、NR4A1、TOP2A、HLA-DRA、IFNG、CA9。
- 无生存/疗效标签（GEO 里没有，须从论文取）；S10 覆盖 3.8%，轨迹投影仍不可行。

## 复制目标（对应 Meylan 的结论）
- R-TCR（Meylan Q3 做不了的，在这里能做）：TCR 多样性与克隆扩增。
- R-Q1：T 细胞浸润（Meylan：8/8 上升）。
- R-Q2：VEGFA 追踪低氧（Meylan：VEGFA 专属，92.2%）。
- R-Q4：髓系无一致变化（Meylan：5 上 3 下，p=0.25）。

## TRB 检验（配对 n=7，预注册）
- T1 多样性：Shannon 熵（归一化）与克隆性（1 − Shannon/log N）、
  top 克隆占比，pre vs post 配对 Wilcoxon（**n=5**，非 7）。
- T2 重叠：Morisita-Horn 指数 + 共享克隆比例（post 克隆中有多少 pre 已存在，
  反之亦然）；扩增克隆是 pre 已有的还是新出现的。
- T3 描述：P13-post 的极端单克隆（top 克隆 0.9999）是否孤立现象；
  各样本克隆数与多样性的分布。
- 零模型：病人内随机置换 pre/post 标签，复算中位差（10k 次）。

## Xenium 检验
- 细胞分型：marker 评分法。T = CD3D/CD3E/CD8A/CD4；髓系 = CD68/CD163；
  B = MS4A1/CD79A；NK = NCAM1/KLRD1；肿瘤 = OLIG2/SOX2/EGFR 高分
  （阈值法，检查分离度与占比是否合理）。
- R-Q1：P17、P20 两对 pre/post 的 T 细胞占比（描述，不做 p 值）；
  另报告 15 个 section 的横断面 T 占比分布。
- R-Q2：肿瘤细胞内 VEGFA 阳性率 vs 低氧分（HIF1A + panel 内 HIF 靶基因）；
  VEGFA+ vs VEGFA− 肿瘤细胞的低氧分，section 内配对比较。
- R-Q4：P17、P20 两对的髓系占比（描述）。
- 边界：Xenium 配对 n=2，不做统计推断；不读疗效；不做轨迹投影。

## 顺序
1. TRB 先行（快、干净、直接补 Meylan 的 TCR 缺口）。
2. Xenium 其后（分型 + R-Q1/R-Q2/R-Q4）。
