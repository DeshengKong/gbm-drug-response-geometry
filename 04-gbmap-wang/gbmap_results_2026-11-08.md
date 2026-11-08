# GBmap 分析结果（2026-11-08）

**数据**：GBmap core H5AD（8.13 GB，校验通过），338,564 细胞 × 27,632 基因，
16 数据集、110 供体。X 为对数归一化表达。
**独立性**：排除 author='Zhao2020'（= GSE148842，GBmap 内仅 PW032 一人，
5,182 细胞），剩余分析与 B 侧（GSE148842）独立。
**方法**：见 analysis_design.md。模块分数 = 模块基因对数表达的简单均值
（未做 scanpy 参考减法，方法差异已记录）；TRAJ = mean(LATE) − mean(EARLY)。
分块读 CSR（20k 细胞/块），仅保留 1,094 个 S10 基因 + 3 配体。
S10 覆盖 1,094/1,200（91%）。

## Check 1：轨迹轴有效性
- 肿瘤细胞（n=125,136）12 维模块分数的 PC1 与 TRAJ：
  Spearman **rho = 0.979**，p ≈ 0；PC1 解释模块分数方差的 45.3%。
- **轨迹轴在 33 万细胞大图谱中成立**，且为模块分数空间的主导方向。
  简单均值版 TRAJ 与 canonical 定义一致性得到验证。

## Check 2：Point 3 复现（供体级，n=62 供体）
每供体肿瘤 TRAJ 均值 vs 该供体髓系细胞配体表达均值：
| 配体 | rho | p |
|---|---|---|
| VEGFA | **0.334** | **0.0079** |
| VCAN | 0.175 | 0.17 |
| ANXA1 | 0.219 | 0.088 |

- **VEGFA 复现**：这是第五个独立数据集（B 体外、Meylan 空间、
  GSE330794 空间、Wang 解离 snRNA、GBmap 图谱）。
- VCAN 在 GBmap 不复现（Wang 中曾显著 rho=0.470）；ANXA1 边缘不显著。

## Check 3：VEGFA–MM2（肿瘤细胞，n=125,136）
按 MM2 中位分高低两组：
- VEGFA 均值：MM2-high **0.937** vs MM2-low 0.218，MW p ≈ 0
- VEGFA 阳性率：MM2-high **66.2%** vs MM2-low 25.0%
- VCAN：MM2-high 0.350 vs MM2-low 0.496——**VCAN 反而在 MM2-low 更高**
- ANXA1：MM2-high 0.772 vs MM2-low 0.199——ANXA1 也富集于 MM2-high

## 结论
- 轨迹轴在 GBmap 大图谱中强验证（rho=0.979）。
- **配体故事收窄到 VEGFA**：VEGFA 在五个数据集中一致地
  （a）与肿瘤轨迹位置相关、（b）在缺氧（MM2-high）肿瘤细胞中富集。
  VCAN 与 ANXA1 跨数据集不一致（VCAN 在 Wang 显著、GBmap/Meylan 不显著；
  ANXA1 在 GBmap 富集于 MM2-high、Meylan 中不区分），
  不宜作为与 VEGFA 并列的"配体分区"证据。
- 与 Meylan 的结论一致：zoning 是 **VEGFA-specific**。
