# -*- coding: utf-8 -*-
"""Build gbm-point1-colab.ipynb for GBM point-1 analysis (trajectory position vs drug-response magnitude).

Reads the local S10 JSON and embeds it into the notebook so the notebook is
self-contained on Colab. Neftel 2019 modules + the 34-gene compact panel are
transcribed verbatim from the author's MATLAB code
(code/b-side-drug-geometry/Glioma_drug_geometry_analysis.m):
  - compact panel: lines 32-35 (marker.AC/OPC/NPC/MES)
  - full Neftel modules: getNeftelModules2019() (line 1265)
No gene list is hand-written from memory.
"""
import json
import nbformat as nbf

BASE = "/home/hatch/workspace/bioinfo_downloads/gbm-point1"
S10_PATH = BASE + "/a-side/trajectory_metamodules_s10.json"
OUT_PATH = BASE + "/gbm-point1-colab.ipynb"

with open(S10_PATH) as f:
    S10 = json.load(f)
S10_JSON_STR = json.dumps(S10)  # embedded verbatim into the notebook

# ---- B-side gene programs, transcribed from the author's MATLAB code ----
PANEL_34 = {
    "AC": ["CST3", "GFAP", "S100B", "HOPX", "SLC1A3", "MLC1", "AQP4"],
    "OPC": ["PLP1", "ALCAM", "OLIG1", "OLIG2", "OMG", "PLLP", "PDGFRA"],
    "NPC": ["SOX4", "SOX11", "DCX", "CD24", "DLL3", "STMN2", "STMN4", "RND3", "DLX5", "DLX6-AS1"],
    "MES": ["VIM", "ANXA1", "ANXA2", "CHI3L1", "CD44", "DDIT3", "ADM", "HILPDA", "ENO2", "LDHA"],
}

NEFTEL_2019 = {
    "AC": ["CST3", "S100B", "SLC1A3", "HEPN1", "HOPX", "MT3", "SPARCL1", "MLC1", "GFAP", "FABP7", "BCAN", "PON2", "METTL7B", "SPARC", "GATM", "RAMP1", "PMP2", "AQP4", "DBI", "EDNRB", "PTPRZ1", "CLU", "PMP22", "ATP1A2", "S100A16", "HEY1", "PCDHGC3", "TTYH1", "NDRG2", "PRCP", "ATP1B2", "AGT", "PLTP", "GPM6B", "F3", "RAB31", "PPAP2B", "ANXA5", "TSPAN7"],
    "OPC": ["BCAN", "PLP1", "GPR17", "FIBIN", "LHFPL3", "OLIG1", "PSAT1", "SCRG1", "OMG", "APOD", "SIRT2", "TNR", "THY1", "PHYHIPL", "SOX2-OT", "NKAIN4", "LPPR1", "PTPRZ1", "VCAN", "DBI", "PMP2", "CNP", "TNS3", "LIMA1", "CA10", "PCDHGC3", "CNTN1", "SCD5", "P2RX7", "CADM2", "TTYH1", "FGF12", "TMEM206", "NEU4", "FXYD6", "RNF13", "RTKN", "GPM6B", "LMF1", "ALCAM", "PGRMC1", "HRASLS", "BCAS1", "RAB31", "PLLP", "FABP5", "NLGN3", "SERINC5", "EPB41L2", "GPR37L1"],
    "NPC1": ["DLL3", "DLL1", "SOX4", "TUBB3", "HES6", "TAGLN3", "NEU4", "MARCKSL1", "CD24", "STMN1", "TCF12", "BEX1", "OLIG1", "MAP2", "FXYD6", "PTPRS", "MLLT11", "NPPA", "BCAN", "MEST", "ASCL1", "BTG2", "DCX", "NXPH1", "HN1", "PFN2", "SCG3", "MYT1", "CHD7", "GPR56", "TUBA1A", "PCBP4", "ETV1", "SHD", "TNR", "AMOTL2", "DBN1", "HIP1", "ABAT", "ELAVL4", "LMF1", "GRIK2", "SERINC5", "TSPAN13", "ELMO1", "GLCCI1", "SEZ6L", "LRRN1", "SEZ6", "SOX11"],
    "NPC2": ["STMN2", "CD24", "RND3", "HMP19", "TUBB3", "MIAT", "DCX", "NSG1", "ELAVL4", "MLLT11", "DLX6-AS1", "SOX11", "NREP", "FNBP1L", "TAGLN3", "STMN4", "DLX5", "SOX4", "MAP1B", "RBFOX2", "IGFBPL1", "STMN1", "HN1", "TMEM161B-AS1", "DPYSL3", "SEPT3", "PKIA", "ATP1B1", "DYNC1I1", "CD200", "SNAP25", "PAK3", "NDRG4", "KIF5A", "UCHL1", "ENO2", "KIF5C", "DDAH2", "TUBB2A", "LBH", "LOC150568", "TCF4", "GNG3", "NFIB", "DPYSL5", "CRABP1", "DBN1", "NFIX", "CEP170", "BLCAP"],
    "MES1": ["CHI3L1", "ANXA2", "ANXA1", "CD44", "VIM", "MT2A", "C1S", "NAMPT", "EFEMP1", "C1R", "SOD2", "IFITM3", "TIMP1", "SPP1", "A2M", "S100A11", "MT1X", "S100A10", "FN1", "LGALS1", "S100A16", "CLIC1", "MGST1", "RCAN1", "TAGLN2", "NPC2", "SERPING1", "C8orf4", "EMP1", "APOE", "CTSB", "C3", "LGALS3", "MT1E", "EMP3", "SERPINA3", "ACTN1", "PRDX6", "IGFBP7", "SERPINE1", "PLP2", "MGP", "CLIC4", "GFPT2", "GSN", "NNMT", "TUBA1C", "GJA1", "TNFRSF1A", "WWTR1"],
    "MES2": ["HILPDA", "ADM", "DDIT3", "NDRG1", "HERPUD1", "DNAJB9", "TRIB3", "ENO2", "AKAP12", "SQSTM1", "MT1X", "ATF3", "NAMPT", "NRN1", "SLC2A1", "BNIP3", "LGALS3", "INSIG2", "IGFBP3", "PPP1R15A", "VIM", "PLOD2", "GBE1", "SLC2A3", "FTL", "WARS", "ERO1L", "XPOT", "HSPA5", "GDF15", "ANXA2", "EPAS1", "LDHA", "P4HA1", "SERTAD1", "PFKP", "PGK1", "EGLN3", "SLC6A6", "CA9", "BNIP3L", "RPL21", "TRAM1", "UFM1", "ASNS", "GOLT1B", "ANGPTL4", "SLC39A14", "CDKN1A", "HSPA9"],
}

nb = nbf.v4.new_notebook()
nb.metadata = {
    "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
    "language_info": {"name": "python", "version": "3.10.0"},
}
cells = []


def md(src):
    cells.append(nbf.v4.new_markdown_cell(src))


def code(src):
    cells.append(nbf.v4.new_code_cell(src))


# ============================ Section 0 ============================
md("""# GBM 点 1：药物响应幅度是不是轨迹位置的函数？

**研究问题**（备忘 `gbm-trajectory-drug-geometry-link-20261105.md` 点 1）：
B（Lavrenteva & Kim 预印本）把胶质瘤状态当离散分箱，发现各状态对同一药物的转录响应**方向高度一致、只在幅度上分歧**；
A（de Jong 等 Nature）说状态其实是一条连续时空轨迹（发育样 → 胶质增生 → 缺氧）。
本 notebook 检验：**响应幅度是否随轨迹位置连续变化，而不是按离散状态分段？**

## 本 notebook 做什么
1. 复现 B 的几何指标 C（方向一致性）与 M（幅度），与作者发表的结果对照 —— 先证明管线算得对，再谈新结论。
2. 把 B 的 252,106 个转化胶质瘤细胞投影到 A 的轨迹轴上（A 的 12 个时间元模块，附表 S10）。
3. 在 43 个治疗–载体比较（TissueID × 药物）中检验：沿轨迹的有向位移、按轨迹位置分箱的位移趋势、幅度分歧在轨迹坐标下是否收敛。
4. 拆分混杂：质心位移里多少是“状态内移动”、多少是“细胞跨状态重分布”。

## 需要什么
- 一个 Google Colab 运行时（CPU 即可，不需要 GPU；AI Ultra 的高内存运行时更稳）。
- 本 notebook 自动下载：B 侧整合 loom（Google Drive 公开链接，约 7 GB）+ 作者的结果 CSV（GitHub raw）。
- A 侧 12 个元模块的基因列表已嵌进本 notebook（附表 S10 提取，1,200 个去重基因），不需要再下载 Nature 附表。

## 不需要什么
- **不需要 62 GB 的 A 侧 h5ad**（E-MTAB-17183 snRNA 对象）：点 1 只需要它的元模块基因列表，不需要它的细胞。
- **不碰 EGA 原始 reads**，不注册、不联系作者。全部是公开处理后数据。
- **不要在本 VM（7 GB 内存）上跑全量**：本 notebook 的目标环境是 Colab；本 VM 只用于 1 万细胞冒烟测试。

## 预计耗时（Colab CPU）
下载 loom 约 5–10 分钟 → MD5 校验约 1 分钟 → 全矩阵流式读（算文库大小）约 20–45 分钟（gzip 解压速度决定；结果缓存到 npz，重跑跳过这步）
→ 读选定基因建 AnnData 约 20–45 分钟（同一矩阵第二遍；也缓存，重跑跳过）→ 打分与几何约 5–10 分钟。
**总计约 60–100 分钟**，按顺序 Run all 即可，中间有断言，失败会停下报错，不会静默出空结果。

## 输出什么
- `point1_outputs/` 下的 4 张图（PNG）+ 4 张结果表（CSV）。
- 最后一段 **verdict**：幅度是不是轨迹位置的连续函数（是 / 否 / 不确定）+ 证据数字 + 对虚拟肿瘤药物模块的含义（连续轨迹输入 vs 离散状态分型）。

## 两条红线（备忘强制要求，全程遵守）
1. **只建模方向 C 与幅度 M，不建模状态分离 K** —— K 已被作者自己证明不稳定（依赖特征空间）。本 notebook 只在复现对照里报告 K，不做任何基于 K 的结论。
2. **绝不把转录响应读成疗效** —— C/M 是转录响应表型，不是活力、克隆形成或临床响应。任何“阻断轨迹能治病”的说法都需要功能与体内验证，本文不提供。
""")

# ============================ Config ============================
code("""# ===== 配置：改这里即可 =====
LOOM_DRIVE_FILE_ID = "18-KInmm43wKdBX95Gq9xbuzAQwtLjgE9"   # B 侧整合 loom（Levitin 注释版，Google Drive 公开链接）
LOOM_MD5 = "9ebdd65747685c92b9ce57929da99a46"             # MANIFEST 已核验
LOOM_EXPECTED_SHAPE = (60725, 391444)                    # (基因, 细胞)
LOOM_LOCAL = "/content/levitin_annotated_gse148842.loom"
GITHUB_RAW = "https://raw.githubusercontent.com/bigbio-el/Glioma-drug-response-geometry/main/results"
OUT_DIR = "/content/point1_outputs"

SMOKE = False        # 设为 True：只取 10,000 个 Tumor 细胞做冒烟测试（快速验证管线，不出全量结论）
SMOKE_N = 10000
RANDOM_SEED = 20261108

# 作者几何的阈值（照抄 MATLAB cfg，不改）
MIN_CELLS_PATIENT_CONTROL = 100   # 每个 TissueID 载体臂最少 Tumor 细胞数
MIN_CELLS_PER_STATE = 25          # 每状态、每臂最少细胞数
MIN_STATES_FOR_METRIC = 3         # 一个比较至少要有 3 个状态可用
EPS = 1e-8

import os
os.makedirs(OUT_DIR, exist_ok=True)
print("OUT_DIR:", OUT_DIR)
print("SMOKE:", SMOKE)
""")

# ============================ Section 1 ============================
md("""## 1. 环境与数据准备

**目的**：装好依赖，拿到 B 侧 loom 与作者的结果表，并逐项核验（MD5、shape、关键字段），对不上就停下报错。
**判断标准**：loom MD5 与 `9ebdd657…` 一致、shape 为 (60725, 391444)；作者两个 CSV 都能读到且列名符合预期。
""")

code("""# Colab：安装依赖（约 2 分钟）
# 本 notebook 用 h5py 直接读 loom（可流式分块读，不把 7 GB 全载入内存），所以不装 loompy。
# scanpy 只用来打分（score_genes）和存 AnnData。
!pip install -q scanpy gdown
import sys
print(sys.version)
""")

code("""# --- 下载 B 侧整合 loom（Google Drive 公开链接，无需登录）---
# 已存在且 MD5 通过时跳过下载（Colab 会话重启后不用重下）。
import os, hashlib

def md5_of(path, chunk=1 << 22):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(chunk), b""):
            h.update(b)
    return h.hexdigest()

need_download = True
if os.path.exists(LOOM_LOCAL):
    print("文件已存在，先核 MD5 ...")
    if md5_of(LOOM_LOCAL) == LOOM_MD5:
        print("MD5 通过，跳过下载")
        need_download = False
    else:
        print("MD5 不符，重新下载")
        os.remove(LOOM_LOCAL)

if need_download:
    # gdown 处理 Google Drive 的 confirm-token 流程；--fuzzy 接受 uc?id= 链接
    !gdown --fuzzy "https://drive.google.com/uc?id={LOOM_DRIVE_FILE_ID}" -O {LOOM_LOCAL}
    md5 = md5_of(LOOM_LOCAL)
    print("下载后 MD5:", md5)
    assert md5 == LOOM_MD5, "MD5 不符：文件可能损坏或下错，停止！请检查网络后重跑本格"
    print("MD5 通过")

print("loom 大小 (GB):", round(os.path.getsize(LOOM_LOCAL) / 1e9, 2))
assert os.path.getsize(LOOM_LOCAL) == 7416810607, "文件大小与 MANIFEST 不符，停止！"
print("大小核验通过")
""")

code("""# --- 下载作者的结果 CSV 作对照（GitHub raw）---
!curl -sL -o {OUT_DIR}/cross_drug_unit_summary.csv {GITHUB_RAW}/cross_drug_unit_summary.csv
!curl -sL -o {OUT_DIR}/cross_drug_global_inference.csv {GITHUB_RAW}/cross_drug_global_inference.csv

import pandas as pd
auth_unit = pd.read_csv(f"{OUT_DIR}/cross_drug_unit_summary.csv")
auth_glob = pd.read_csv(f"{OUT_DIR}/cross_drug_global_inference.csv")
print("auth_unit:", auth_unit.shape)
print(auth_unit.head(12).to_string())
print()
print("auth_glob:", auth_glob.shape)
print(auth_glob.to_string())
assert list(auth_unit.columns)[:3] == ["patient", "n_drug_pairs", "median_within_drug_C"], "作者 unit 表列名变了，停止！"
assert list(auth_glob.columns)[:3] == ["comparison", "n_units", "median_delta"], "作者 global 表列名变了，停止！"
print("作者结果表核验通过")
""")

# ---- embed gene programs ----
s10_cell = (
    "# --- A 侧 12 个时间元模块（附表 S10，Table S10 = Cell2fate recurrent module annotation）---\n"
    "# 来源：Nature 论文附表 MOESM5（41586_2026_11234_MOESM5_ESM.xlsx）Table S10，\n"
    "# 已提取为本地 a-side/trajectory_metamodules_s10.json，本格把 JSON 原文嵌进来，\n"
    "# 所以本 notebook 不需要再下载 Nature 附表。共 12 个模块、1,368 条目、1,200 个去重基因。\n"
    "import json\n"
    "S10_JSON = " + repr(S10_JSON_STR) + "\n"
    "S10 = json.loads(S10_JSON)\n"
    "assert len(S10) == 12, 'S10 应为 12 个模块'\n"
    "MM_GENES = {m['metamodule']: m['genes'] for m in S10}\n"
    "MM_ANNOT = {m['metamodule']: m['annotation'] for m in S10}\n"
    "n_dedup = len(set(g for v in MM_GENES.values() for g in v))\n"
    "print('S10 模块数:', len(S10), '| 去重基因数:', n_dedup)\n"
    "assert n_dedup == 1200, 'S10 去重基因数应为 1200'\n"
    "for m in S10:\n"
    "    print(f\"  MM{m['metamodule']:>2} {m['annotation']:<24} n={m['n_genes']}\")\n"
)
code(s10_cell)

neftel_cell = (
    "# --- B 侧基因程序（逐字转录自作者 MATLAB 代码，不凭记忆手写） ---\n"
    "# 来源：code/b-side-drug-geometry/Glioma_drug_geometry_analysis.m\n"
    "#   * 34 基因紧凑 panel：marker.AC/OPC/NPC/MES（第 32-35 行）\n"
    "#   * 完整 Neftel 2019 六模块：getNeftelModules2019()（第 1265 行起）\n"
    "# B 的状态定义：紧凑 panel 打分取 argmax；完整模块先算 6 个模块分，再合成 4 维\n"
    "# S4 = [AC, OPC, max(NPC1,NPC2), max(MES1,MES2)] 取 argmax（assignNeftelFourStates）。\n"
    "PANEL_34 = " + repr(PANEL_34) + "\n"
    "NEFTEL_2019 = " + repr(NEFTEL_2019) + "\n"
    "STATES = ['AC', 'OPC', 'NPC', 'MES']\n"
    "n_panel = sum(len(v) for v in PANEL_34.values())\n"
    "n_neftel_unique = len(set(g for v in NEFTEL_2019.values() for g in v))\n"
    "print('紧凑 panel 基因数:', n_panel)\n"
    "print('Neftel 六模块去重基因数:', n_neftel_unique)\n"
    "assert n_panel == 34, '紧凑 panel 应为 34 基因'\n"
    "for k, v in PANEL_34.items():\n"
    "    print(f'  panel {k}: {len(v)} 基因')\n"
    "for k, v in NEFTEL_2019.items():\n"
    "    print(f'  Neftel {k}: {len(v)} 基因')\n"
)
code(neftel_cell)

# ============================ Section 2 ============================
md("""## 2. 读 loom 建 AnnData

**目的**：从 loom 只取出 Tumor 细胞（252,106），建成 log 归一化的 AnnData，后面全部分析都在这个小对象上做。
**判断标准**：loom shape 为 (60725, 391444)；`Cell_type` / `TissueID` / `Treatment` / `Sequencing_Name` 四个 col_attrs 都存在；
Tumor 细胞数 = 252,106（冒烟模式下 = 10,000）；QC 用矩阵实际 nnz（不要用 obs 里现成的指标 —— MANIFEST 已证实 obs 的 n_genes 与实际 nnz 对不上）。
基因名有重复符号（60,725 行、58,828 个唯一符号，主要是 Y_RNA、snoRNA 等多拷贝非编码 RNA 家族；必需基因里只有 SERPINA3 一个重复），
一律按首次出现取用（与 MATLAB ismember 的行为一致）。

**关键技术点**：loom 里没有文库大小字段（作者代码也说明了这一点），所以文库大小必须从全矩阵流式算出来。
本格用 h5py 分块读 `/matrix`（基因块 × 细胞块），累积每细胞的总 counts 与非零基因数，**任何时刻都不把 7 GB 全载入内存**。
这一步在 Colab 上约 20–45 分钟（gzip 解压速度决定；结果会缓存到 `point1_outputs/lib_nnz.npz`，重跑直接读缓存，不用重算）。
""")

code("""# --- 打开 loom，核验 shape 与字段（MANIFEST §2.1 字段词表）---
import h5py, numpy as np, pandas as pd

f = h5py.File(LOOM_LOCAL, "r")
mat = f["/matrix"]
print("matrix shape (基因 x 细胞):", mat.shape)
print("matrix dtype:", mat.dtype, "| 压缩:", mat.compression, "| chunk:", mat.chunks)
assert mat.shape == LOOM_EXPECTED_SHAPE, "loom shape 与 MANIFEST 不符，停止！"
assert mat.shape[0] == 60725 and mat.shape[1] == 391444

def read_str(path):
    arr = f[path][:]
    return np.array([x.decode("utf-8") if isinstance(x, bytes) else str(x) for x in arr])

genes = read_str("/row_attrs/Gene")
celltype_raw = read_str("/col_attrs/Cell_type")
tissue = read_str("/col_attrs/TissueID")
treatment_raw = read_str("/col_attrs/Treatment")
seq_name = read_str("/col_attrs/Sequencing_Name")
print("基因数:", len(genes), "| 细胞数:", len(celltype_raw))
# 注意：loom 的 Gene 有重复符号（60,725 行、58,828 个唯一符号；313 个符号重复、共 2,210 行，
# 主要是 Y_RNA、snoRNA、U1/U2/U6 等多拷贝非编码 RNA 家族）。这不是数据损坏，是上游注释的固有属性。
# 处理规则：基因名 -> 行号只保留首次出现（与 MATLAB ismember 的行为一致）。
# 经核查，必需基因里只有一个重复：SERPINA3（Neftel MES1，2 行），按首次出现取用；其余重复符号都不是必需基因。
import collections
_gene_counts = collections.Counter(genes)
dup_syms = sorted(g for g, c in _gene_counts.items() if c > 1)
print(f"重复基因符号：{len(dup_syms)} 个，共 {sum(_gene_counts[g] for g in dup_syms)} 行；唯一符号数：{len(_gene_counts)}")
assert len(_gene_counts) == 58828, "唯一基因符号数与预期 58,828 不符，停止！"
_needed = set()
for _v in MM_GENES.values():
    _needed.update(_v)
for _v in NEFTEL_2019.values():
    _needed.update(_v)
for _v in PANEL_34.values():
    _needed.update(_v)
_needed_dups = sorted(_needed & set(dup_syms))
print("必需基因里的重复符号：", _needed_dups)
assert _needed_dups == ["SERPINA3"], f"必需基因里出现未预期的重复符号 {_needed_dups}，停下检查！"
print("基因名重复情况已确认（仅 SERPINA3 一个必需基因重复，取首次出现）")
print()
print("Cell_type 值分布：")
print(pd.Series(celltype_raw).value_counts().to_string())
print()
print("TissueID 数:", len(np.unique(tissue)), "| Treatment 值:", sorted(set(treatment_raw)))
print("Treatment 值分布：")
print(pd.Series(treatment_raw).value_counts().to_string())
""")

code("""# --- 药物名归一化（照抄作者 normalizeTreatment 的映射规则，不用正则）---
# DMSO -> vehicle；panobin -> panobinostat；etopos -> etoposide；topotec -> topotecan；
# ro492 -> RO4929097；tazemet -> tazemetostat；ispin/ispen -> ispinesib；ana-12/ana12 -> ANA-12；
# givinostat / vorinostat / temozolomide 保持原名；其余去剂量后缀后小写保留。
def normalize_treatment(raw):
    r = raw.strip().lower()
    if "vehicle" in r or "dmso" in r or r in ("control", "ctrl", "none"):
        return "vehicle"
    if "panobin" in r:
        return "panobinostat"
    if "etopos" in r:
        return "etoposide"
    if "topotec" in r:
        return "topotecan"
    if "ro492" in r or "r0492" in r:
        return "RO4929097"
    if "tazemet" in r:
        return "tazemetostat"
    if "ispin" in r or "ispen" in r:
        return "ispinesib"
    if "ana-12" in r or "ana12" in r:
        return "ANA-12"
    if "givinostat" in r:
        return "givinostat"
    if "vorinostat" in r:
        return "vorinostat"
    if "temozolomide" in r or r == "tmz":
        return "temozolomide"
    return r

treatment = np.array([normalize_treatment(x) for x in treatment_raw])
print("归一化后 Treatment:", sorted(set(treatment)))
print(pd.Series(treatment).value_counts().to_string())
assert "vehicle" in set(treatment), "没有 vehicle 对照，停止！"

# Tumor = 转化/恶性细胞（作者规则：celltype 含 transform/malig/neoplast/tumor/glioma）
is_tumor = np.array([("tumor" in c.lower()) or ("transform" in c.lower())
                     or ("malig" in c.lower()) or ("neoplast" in c.lower())
                     or ("glioma" in c.lower()) for c in celltype_raw])
print()
print("Tumor 细胞数:", int(is_tumor.sum()))
assert int(is_tumor.sum()) == 252106, "Tumor 细胞数与 MANIFEST (252,106) 不符，停止！"
tumor_idx_all = np.where(is_tumor)[0]
print("TissueID (Tumor) 数:", len(np.unique(tissue[is_tumor])))
assert len(np.unique(tissue[is_tumor])) == 19, "Tumor 的 TissueID 应为 19 个"
print("每个 TissueID 的 Tumor 细胞数：")
print(pd.Series(tissue[is_tumor]).value_counts().sort_index().to_string())

# 冒烟模式：取按列序排在最前面的 10,000 个 Tumor 细胞（loom 里 Tumor 连续占第 0..252105 列，
# 所以这 10,000 个也在一个连续列区间里，流式读走连续 hyperslab、很快；只验证管线，不出结论）
if SMOKE:
    tumor_idx = tumor_idx_all[:SMOKE_N]
    print("SMOKE 模式：取按列序最前面的 Tumor 细胞", len(tumor_idx), "个")
else:
    tumor_idx = tumor_idx_all
print("本次分析用细胞数:", len(tumor_idx))
""")

code("""# --- 流式算文库大小与 nnz（全矩阵分块读；结果缓存到 npz）---
# 注意：loom 没有文库大小字段，必须从全矩阵算（作者代码也是这样做的）。
# 正常模式：基因块 256 x 细胞块 65536 的瓦片循环，只累积 Tumor 列。
# 冒烟模式：只对抽样的 10,000 列做花式索引读（快得多）。
import os
cache_tag = "smoke" if SMOKE else "full"
cache_path = f"{OUT_DIR}/lib_nnz_{cache_tag}.npz"
n_genes, n_cells = mat.shape
cmin, cmax = int(tumor_idx.min()), int(tumor_idx.max()) + 1
contiguous = (cmax - cmin) == len(tumor_idx)

if os.path.exists(cache_path):
    z = np.load(cache_path)
    lib = z["lib"]; nnz = z["nnz"]
    print("读缓存:", cache_path)
else:
    lib = np.zeros(len(tumor_idx), dtype=np.float64)
    nnz = np.zeros(len(tumor_idx), dtype=np.int64)
    print(f"列区间 [{cmin}, {cmax})，细胞数 {len(tumor_idx)}，连续：{contiguous}")
    if contiguous:
        # 连续列区间：基因块 x 连续列（快；实测 Tumor 连续占第 0..252105 列）
        # 注意：不要用 h5py 花式索引 mat[genes, tumor_idx] —— 实测 512x10k 花式读要 45 秒/次，全量要 90 分钟，是病态的。
        GB = 512 if SMOKE else 128
        for g0 in range(0, n_genes, GB):
            g1 = min(g0 + GB, n_genes)
            blk = mat[g0:g1, cmin:cmax]
            lib += blk.sum(axis=0)
            nnz += (blk > 0).sum(axis=0)
            if (g0 // GB) % 20 == 0:
                print(f"  基因块 {g0}-{g1} / {n_genes}")
    else:
        # 退路：基因块 x 细胞块瓦片循环 + 布尔掩码（慢但正确）
        tumor_mask = np.zeros(n_cells, dtype=bool)
        tumor_mask[tumor_idx] = True
        GB, CB = 256, 65536
        for g0 in range(0, n_genes, GB):
            g1 = min(g0 + GB, n_genes)
            for c0 in range(0, n_cells, CB):
                c1 = min(c0 + CB, n_cells)
                blk = mat[g0:g1, c0:c1]
                m = tumor_mask[c0:c1]
                if not m.any():
                    continue
                sub = blk[:, m]
                lib += sub.sum(axis=0)
                nnz += (sub > 0).sum(axis=0)
            if (g0 // GB) % 40 == 0:
                print(f"  基因块 {g0}-{g1} / {n_genes}")
    np.savez(cache_path, lib=lib, nnz=nnz)
    print("已缓存:", cache_path)

print("文库大小：中位数", int(np.median(lib)), "| 均值", round(float(np.mean(lib)), 1))
print("nnz（每细胞非零基因数）：中位数", int(np.median(nnz)), "| 均值", round(float(np.mean(nnz)), 1))
assert np.median(lib) > 1000, "文库大小中位数异常低，停止！"
assert np.median(nnz) > 500, "nnz 中位数异常低，停止！"
assert (lib > 0).all(), "存在文库大小为 0 的细胞，停止！"
print("文库大小 / nnz 流式计算通过")
""")

code("""# --- 选基因：S10 + Neftel + 紧凑 panel 的并集，只在 loom 里存在的基因 ---
# 基因名 -> 行号：只保留首次出现（与 MATLAB ismember 一致；见 §2 的重复符号说明）
gene_to_row = {}
for _i, _g in enumerate(genes):
    if _g not in gene_to_row:
        gene_to_row[_g] = _i
assert len(gene_to_row) == 58828, "唯一基因符号数与预期 58,828 不符，停止！"
all_needed = set()
for v in MM_GENES.values():
    all_needed.update(v)
for v in NEFTEL_2019.values():
    all_needed.update(v)
for v in PANEL_34.values():
    all_needed.update(v)
present = sorted(g for g in all_needed if g in gene_to_row)
missing = sorted(all_needed - set(present))
print("需要基因:", len(all_needed), "| loom 中存在:", len(present), "| 缺失:", len(missing))
print("缺失基因示例:", missing[:20])

# 每模块覆盖率（低于 50% 的模块要警惕）
print()
print("S10 每模块覆盖率：")
for mm in sorted(MM_GENES):
    tot = len(MM_GENES[mm])
    cov = sum(1 for g in MM_GENES[mm] if g in gene_to_row)
    print(f"  MM{mm:>2} {MM_ANNOT[mm]:<24} {cov}/{tot} = {cov/tot:.0%}")
    assert cov / tot >= 0.5, f"MM{mm} 覆盖率低于 50%，停止！"
print()
print("Neftel / panel 覆盖率：")
for k, v in NEFTEL_2019.items():
    cov = sum(1 for g in v if g in gene_to_row)
    print(f"  Neftel {k}: {cov}/{len(v)} = {cov/len(v):.0%}")
    assert cov / len(v) >= 0.7, f"Neftel {k} 覆盖率低于 70%，停止！"
for k, v in PANEL_34.items():
    cov = sum(1 for g in v if g in gene_to_row)
    print(f"  panel {k}: {cov}/{len(v)}")
    assert cov == len(v), f"紧凑 panel {k} 有基因缺失，停止！"
print("基因覆盖率核验通过")
sel_rows = np.array([gene_to_row[g] for g in present], dtype=np.int64)
# h5py 花式索引要求行号递增：按行号重排 present/sel_rows（AnnData 的 var 顺序跟着走即可）
_order = np.argsort(sel_rows)
sel_rows = sel_rows[_order]
present = [present[i] for i in _order]
assert (np.diff(sel_rows) > 0).all(), "选定基因行号不唯一或不递增，停止！"
""")

code("""# --- 读选定基因的 Tumor 细胞，log 归一化，建 AnnData（缓存到 h5ad）---
# 归一化与作者一致：Xlog = log1p(1e4 * X / lib)，lib 用上一步从全矩阵算的文库大小。
import scipy.sparse as sp
import anndata as ad

adata_path = f"{OUT_DIR}/adata_tumor_selected_{'smoke' if SMOKE else 'full'}.h5ad"
if os.path.exists(adata_path):
    adata = ad.read_h5ad(adata_path)
    print("读缓存 AnnData:", adata_path)
else:
    # 瓦片循环读（与上一步相同的走法）：基因块 x 连续 Tumor 列区间，块内再取选定基因行。
    # 不要用 mat[sel_rows] 花式行读 —— 选定基因散落在全基因范围，花式行读要解压数倍于结果的数据，实测极慢。
    try:
        _ = mat.shape
    except (NameError, ValueError):
        f = h5py.File(LOOM_LOCAL, "r")
        mat = f["/matrix"]
    sel_set = np.zeros(mat.shape[0], dtype=bool)
    sel_set[sel_rows] = True
    GB = 512 if SMOKE else 128
    chunks = []
    for g0 in range(0, mat.shape[0], GB):
        g1 = min(g0 + GB, mat.shape[0])
        local = np.where(sel_set[g0:g1])[0]
        if len(local) == 0:
            continue
        blk = mat[g0:g1, cmin:cmax][local, :]
        blk = blk.astype(np.float32)
        blk = np.log1p(1e4 * blk / lib[None, :])       # 按列（细胞）除文库大小
        chunks.append(sp.csr_matrix(blk.T))            # 转置成 细胞 x 基因，AnnData 要求观测在行
        if (g0 // GB) % 20 == 0:
            print(f"  基因块 {g0}-{g1} / {mat.shape[0]}，已取 {len(chunks)} 块")
    X = sp.hstack(chunks, format="csr")   # 每块是 细胞 x 块内基因，横向拼成 细胞 x 全部选定基因
    del chunks
    obs = pd.DataFrame({
        "TissueID": tissue[tumor_idx],
        "Treatment": treatment[tumor_idx],
        "Sequencing_Name": seq_name[tumor_idx],
        "lib": lib,
        "nnz": nnz,
    }, index=[f"cell_{i}" for i in tumor_idx])
    var = pd.DataFrame(index=pd.Index(present, name="gene"))
    adata = ad.AnnData(X=X, obs=obs, var=var)
    adata.write_h5ad(adata_path)
    print("已缓存 AnnData:", adata_path)

print(adata)
print("X nnz:", int(adata.X.nnz), "| 密度:", round(float(adata.X.nnz) / (adata.n_obs * adata.n_vars), 4))
expected_n = SMOKE_N if SMOKE else 252106
assert adata.n_obs == expected_n, f"AnnData 细胞数应为 {expected_n}，实际 {adata.n_obs}"
assert adata.n_vars == len(present)
assert (adata.obs["lib"].values > 0).all()
print("AnnData 构建通过")
print()
print("每 TissueID x Treatment 的 Tumor 细胞数（前 20 行）：")
print(pd.crosstab(adata.obs["TissueID"], adata.obs["Treatment"]).to_string())
""")

# ============================ Section 3 ============================
md("""## 3. 状态重映射：B 的状态 ↔ A 的 12 个元模块

**目的**：建一座桥 —— B 的离散状态（紧凑 panel 的 AC/OPC/NPC/MES，以及完整 Neftel 六模块合成的四状态）
与 A 的 12 个元模块分数之间，到底对应得上吗？对应不上，后面的轨迹分析就站不住。
**判断标准**：
- B 自己的两种状态定义（紧凑 panel vs 完整 Neftel）一致率应与作者报告的数值接近（作者说两套方法中等一致，约 50% 左右 —— 注意：作者原文说 Neftel 完整模块与紧凑 panel 的一致率就是中等水平，这是已知事实，不是本 notebook 的失败）。
- B 的 MES 状态应富集 A 的晚期模块（Gliosis-hypoxia、AC-gliosis-like），NPC/OPC 状态应富集早期模块（NPC/OPC 样模块）；
  若完全反过来或毫无结构，**停下**，先检查基因映射，不要往下做。
- 基因重叠：B 的状态基因集与 A 的模块基因集之间应有可解释的重叠（Jaccard），而不是全零。
""")

code("""# --- 患者内 z-score（照抄作者逻辑）：每基因、每 TissueID，以本患者 vehicle 细胞为参考 ---
# Z = (Xlog - mu_vehicle) / sd_vehicle；sd < 1e-6 时按 1 处理（作者同款 floor）。
# 只对 panel + Neftel 的并集基因算（约 280 个基因），不需要整矩阵。
import numpy as np

def patient_zscores(adata, gene_list):
    # 返回 (n_cells x len(gene_list)) 的 float32 稠密矩阵，按 TissueID 分别 z-score。
    idx = [adata.var_names.get_loc(g) for g in gene_list]
    X = adata[:, idx].X
    X = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
    X = X.astype(np.float32)
    Z = np.zeros_like(X)
    tid = adata.obs["TissueID"].values
    trt = adata.obs["Treatment"].values
    for p in np.unique(tid):
        sel = tid == p
        ref = sel & (trt == "vehicle")
        assert ref.sum() >= MIN_CELLS_PATIENT_CONTROL, f"{p} 的 vehicle Tumor 细胞不足 {MIN_CELLS_PATIENT_CONTROL}"
        mu = X[ref].mean(axis=0)
        sd = X[ref].std(axis=0)
        sd = np.where(sd < 1e-6, 1.0, sd)
        Z[sel] = (X[sel] - mu) / sd
    return Z

panel_genes = sorted(set(g for v in PANEL_34.values() for g in v))
neftel_genes = sorted(set(g for v in NEFTEL_2019.values() for g in v))
# 只对 loom 里存在的基因算 z-score（Neftel 表里有个别旧符号 loom 里没有，如 MES2 的 ERO1L；
# 覆盖率已在 §2 核过：Neftel 各模块 97–100%，缺 1 个基因不影响模块均值）
panel_genes = [g for g in panel_genes if g in adata.var_names]
neftel_genes = [g for g in neftel_genes if g in adata.var_names]
print("panel 基因（loom 中存在）:", len(panel_genes), "| Neftel 去重基因（loom 中存在）:", len(neftel_genes))

Z_panel = patient_zscores(adata, panel_genes)
Z_neftel = patient_zscores(adata, neftel_genes)
print("Z_panel:", Z_panel.shape, "| Z_neftel:", Z_neftel.shape)
assert np.isfinite(Z_panel).all() and np.isfinite(Z_neftel).all()
print("患者内 z-score 完成")
""")

code("""# --- B 状态（两种定义）---
# 紧凑 panel：每状态 = 该状态基因的平均 Z，取 argmax；margin = 最高 - 次高。
# 完整 Neftel：6 个模块分 -> S4 = [AC, OPC, max(NPC1,NPC2), max(MES1,MES2)] -> argmax。
def module_means(Z, gene_list, modules):
    out = {}
    pos = {g: i for i, g in enumerate(gene_list)}
    for name, gs in modules.items():
        cols = [pos[g] for g in gs if g in pos]
        out[name] = Z[:, cols].mean(axis=1)
    return out

panel_scores = module_means(Z_panel, panel_genes, PANEL_34)
neftel_scores = module_means(Z_neftel, neftel_genes, NEFTEL_2019)

panel_mat = np.vstack([panel_scores[s] for s in STATES]).T
state_compact_idx = panel_mat.argmax(axis=1)
sorted_panel = np.sort(panel_mat, axis=1)
margin = sorted_panel[:, -1] - sorted_panel[:, -2]
adata.obs["state_compact"] = pd.Categorical([STATES[i] for i in state_compact_idx], categories=STATES)
adata.obs["state_margin"] = margin

s4 = np.vstack([
    neftel_scores["AC"],
    neftel_scores["OPC"],
    np.maximum(neftel_scores["NPC1"], neftel_scores["NPC2"]),
    np.maximum(neftel_scores["MES1"], neftel_scores["MES2"]),
]).T
state_neftel_idx = s4.argmax(axis=1)
adata.obs["state_neftel"] = pd.Categorical([STATES[i] for i in state_neftel_idx], categories=STATES)
for s in STATES:
    adata.obs[f"panel_{s}"] = panel_mat[:, STATES.index(s)]
    adata.obs[f"neftel_{s}"] = s4[:, STATES.index(s)]

print("紧凑 panel 状态分布：")
print(adata.obs["state_compact"].value_counts().to_string())
print()
print("完整 Neftel 状态分布：")
print(adata.obs["state_neftel"].value_counts().to_string())
print()
print("margin 分布：中位数", round(float(np.median(margin)), 3),
      "| margin>0.25 的比例（作者 confident 阈值）:", round(float((margin > 0.25).mean()), 3))
adata.write_h5ad(adata_path)
print("状态已写入 AnnData")
""")

code("""# --- A 的 12 个元模块打分（scanpy score_genes，在 log 归一化 X 上）---
import scanpy as sc
for mm in sorted(MM_GENES):
    gs = [g for g in MM_GENES[mm] if g in adata.var_names]
    sc.tl.score_genes(adata, gs, score_name=f"MM{mm}_{MM_ANNOT[mm]}", random_state=0, use_raw=False)
    print(f"MM{mm} {MM_ANNOT[mm]}: 用 {len(gs)}/{len(MM_GENES[mm])} 个基因打分完成")
mm_cols = [c for c in adata.obs.columns if c.startswith("MM")]
print("模块分数列:", len(mm_cols))
assert len(mm_cols) == 12
print(adata.obs[mm_cols].describe().T[["mean", "std", "min", "max"]].round(3).to_string())
adata.write_h5ad(adata_path)
print("12 个元模块分数已写入 AnnData")
""")

code("""# --- 列联表 / 重叠系数：B 状态 vs A 模块 ---
from sklearn.metrics import cohen_kappa_score
import numpy as np

# 1) B 的两种状态定义是否一致（作者报告约中等一致）
agree = (adata.obs["state_compact"] == adata.obs["state_neftel"]).mean()
kappa = cohen_kappa_score(adata.obs["state_compact"], adata.obs["state_neftel"])
print(f"紧凑 panel vs 完整 Neftel：一致率 {agree:.1%}，Cohen kappa {kappa:.3f}")
print("（作者原文：两套方法中等一致 —— 若本格一致率远低于 30% 或高于 90%，先检查基因映射）")
print()
print("列联表（行=紧凑 panel，列=完整 Neftel）：")
ct = pd.crosstab(adata.obs["state_compact"], adata.obs["state_neftel"], normalize="index").round(2)
print(ct.to_string())
assert agree >= 0.30, "两种 B 状态定义一致率过低，基因映射可能有问题，停下检查！"

# 2) 每 B 状态的平均 S10 模块分（按模块标准化后看结构）
mm_cols = [c for c in adata.obs.columns if c.startswith("MM")]
means = adata.obs.groupby("state_compact", observed=True)[mm_cols].mean()
zmeans = (means - means.mean(axis=0)) / means.std(axis=0)
print()
print("每 B 状态的平均 S10 模块分（按列 z-score，便于看结构）：")
print(zmeans.round(2).to_string())

# 3) 基因重叠 Jaccard：B 状态基因集 vs A 模块
def jaccard(a, b):
    a, b = set(a), set(b)
    return len(a & b) / max(1, len(a | b))
print()
print("Jaccard（行=B 状态基因集，列=A 模块；只显示 >= 0.05 的格子）：")
b_sets = {}
for s in STATES:
    b_sets[s] = set(PANEL_34[s]) | set(NEFTEL_2019.get(s, [])) | set(NEFTEL_2019.get(s + "1", [])) | set(NEFTEL_2019.get(s + "2", []))
for s in STATES:
    row = []
    for mm in sorted(MM_GENES):
        j = jaccard(b_sets[s], MM_GENES[mm])
        row.append(f"{j:.2f}" if j >= 0.05 else " .  ")
    print(f"  {s:4} " + " ".join(row))
print("     列序：MM0..MM11")
print()
print("判断：若 MES 行在 MM2/MM3/MM5 列有明显重叠、NPC/OPC 行在 MM1/MM4/MM6/MM7/MM11 列有重叠，")
print("      则 B 状态与 A 轨迹在基因层面有对应关系，可以往下做；否则先停下检查。")
""")

# ============================ Section 4 ============================
md("""## 4. 轨迹打分：一维轨迹位置的定义

**目的**：把 12 维模块分数压缩成一维“轨迹位置”—— 每个细胞在“发育样 → 胶质增生 → 缺氧”这条轴上走到哪了。
**判断标准**：TRAJ 与数据驱动的 PC1 轴高度相关（Spearman > 0.7 视为定义合理）；MM9 归入早期 vs 排除的敏感性分析结果方向一致。

> **显式假设（需要用户拍板）**：S10 的 12 个模块本身没有“早/晚”顺序标签。
> 本 notebook 按模块注释 + 论文 `factor_id_map.tsv`（code/a-side-GBM_space_paper/metadata/2.TM/factor_id_map.tsv；
> 该表映射的是*空间* TM 因子到环境，Dev-like (AC/OPC/Neuron) → Gliosis → Hypoxic，与论文正文的时序一致）
> 做如下分组 —— **这是人为分组，不是论文直接给出的标签**：
> - **早期（发育样）**：MM1 NPC-neuronal-like 1、MM4 NPC-OPC-like 2、MM6 NPC-neuronal-like 2、MM7 NPC-OPC-like 1、MM8 AC-NPC-OPC-like 1、MM10 NPC-neuronal-like 3、MM11 OPC-like
> - **晚期（胶质增生/缺氧）**：MM2 Gliosis-hypoxia、MM3 AC-gliosis-like 2、MM5 AC-gliosis-like 1
> - **主定义排除**：MM0 Cycling（周期与轨迹正交，论文中细胞周期是独立维度）、MM9 AC-progenitor-like（祖细胞样，位置不确定）
> - **敏感性分析**：(a) MM9 归入早期；(b) 12 维分数的 PC1 作为数据驱动轴（符号定向为与 TRAJ 正相关）。
> 若用户认为 MM9 应属晚期、或 MM0 不应排除，改下面 `EARLY` / `LATE` 两个列表重跑即可。
""")

code("""# --- 轨迹位置 TRAJ = mean(晚期模块分) - mean(早期模块分) ---
import numpy as np
from scipy.stats import spearmanr
from sklearn.decomposition import PCA

MM_ORDER = list(range(12))  # MM0..MM11
EARLY = [1, 4, 6, 7, 8, 10, 11]
LATE = [2, 3, 5]
EXCLUDED = [0, 9]  # MM0 Cycling, MM9 AC-progenitor-like（主定义排除）
print("早期模块:", [(m, MM_ANNOT[m]) for m in EARLY])
print("晚期模块:", [(m, MM_ANNOT[m]) for m in LATE])
print("主定义排除:", [(m, MM_ANNOT[m]) for m in EXCLUDED])

S = adata.obs[[f"MM{m}_{MM_ANNOT[m]}" for m in MM_ORDER]].values.astype(np.float32)
assert S.shape == (adata.n_obs, 12)

def traj_score(S, early, late):
    return S[:, late].mean(axis=1) - S[:, early].mean(axis=1)

adata.obs["TRAJ"] = traj_score(S, EARLY, LATE)
adata.obs["TRAJ_MM9EARLY"] = traj_score(S, EARLY + [9], LATE)  # 敏感性 (a)

# 敏感性 (b)：12 维分数的 PC1（标准化后），符号定向为与 TRAJ 正相关
S_z = (S - S.mean(axis=0)) / (S.std(axis=0) + 1e-8)
pc1 = PCA(n_components=1, random_state=0).fit_transform(S_z).ravel()
if np.corrcoef(pc1, adata.obs["TRAJ"].values)[0, 1] < 0:
    pc1 = -pc1
adata.obs["TRAJ_PC1"] = pc1

r1, p1 = spearmanr(adata.obs["TRAJ"], adata.obs["TRAJ_MM9EARLY"])
r2, p2 = spearmanr(adata.obs["TRAJ"], adata.obs["TRAJ_PC1"])
print(f"TRAJ vs TRAJ_MM9EARLY：Spearman rho = {r1:.3f} (p={p1:.2e})")
print(f"TRAJ vs TRAJ_PC1：      Spearman rho = {r2:.3f} (p={p2:.2e})")
assert r2 > 0.7, "TRAJ 与 PC1 相关性低于 0.7 —— 晚期/早期分组可能有问题，停下检查！"
print("轨迹位置定义通过合理性检查（TRAJ 与数据驱动 PC1 高度一致）")
print()
print("TRAJ 分布：均值", round(float(adata.obs["TRAJ"].mean()), 3),
      "| 标准差", round(float(adata.obs["TRAJ"].std()), 3))
print("按 Treatment 的 TRAJ 均值（仅供参考，不作结论）：")
print(adata.obs.groupby("Treatment")["TRAJ"].agg(["mean", "count"]).round(3).to_string())
adata.write_h5ad(adata_path)

# 轨迹轴的单位向量（12 维模块空间）：晚期 +1，早期 -1，其余 0，再归一化
u_traj = np.zeros(12)
u_traj[LATE] = 1.0
u_traj[EARLY] = -1.0
u_traj = u_traj / np.linalg.norm(u_traj)
print("u_traj =", np.round(u_traj, 3))
""")

# ============================ Section 5 ============================
md("""## 5. 核心检验：响应幅度是不是轨迹位置的函数？

**目的**：在 43 个治疗–载体比较（TissueID × 药物）中逐个检验点 1。
**方法**：
- **(0) 先复现 B**：在 34 基因 panel 的 Z 空间里算 C（方向一致性）与 M（幅度），与作者结果对照。
  预期：43 个比较的 C 中位数 ≈ 0.608、43/43 为正；每患者的中位数与作者 `cross_drug_unit_summary.csv` 的
  `median_within_drug_C` 接近（Pearson r > 0.8）。对不上说明管线有问题，**先停下修管线，不要信后面的轨迹结果**。
- **(a) 沿轨迹的有向位移**：在 12 维模块分数空间算载体 vs 治疗的质心位移向量，投影到轨迹轴 u_traj 上，
  得有向位移 d_traj（正 = 向晚期移动，负 = 向早期移动）与位移幅度 m12。
- **(b) 轨迹分箱趋势**：按载体细胞的 TRAJ 五分位分箱，每箱内算位移，看位移幅度/方向是否随轨迹位置平滑变化（Spearman 趋势检验）。
- **(c) 幅度分歧是否收敛**：同一空间（12 维模块空间）里比较离散状态分箱 vs 轨迹分箱的幅度离散度 ——
  若轨迹分箱的幅度离散度显著小于状态分箱，说明 B 的“幅度分歧”部分是离散分箱的人为产物，轨迹坐标更能解释它。

**红线**：只建模 C 与 M，不建模 K（K 只在复现对照里报告数值，不做结论）；绝不把转录响应读成疗效。
""")

code("""# --- (0) 复现 B：C / M / K / TV（照抄作者 computeGeometryLocal 的定义）---
# 34 基因 panel 的 Z 空间；每 (TissueID, 药物) 比较：每状态在两臂各 >=25 细胞才可用，至少 3 个状态。
import numpy as np

def pairwise_mean_cos(D):
    # D: (n_states, n_genes)，返回 Delta 向量间的平均两两余弦
    n = D.shape[0]
    vals = []
    for i in range(n):
        for j in range(i + 1, n):
            ni, nj = np.linalg.norm(D[i]), np.linalg.norm(D[j])
            if ni < 1e-12 or nj < 1e-12:
                continue
            vals.append(float(D[i] @ D[j] / (ni * nj)))
    return float(np.mean(vals)) if vals else np.nan

def pairwise_mean_dist(C):
    n = C.shape[0]
    vals = [float(np.linalg.norm(C[i] - C[j])) for i in range(n) for j in range(i + 1, n)]
    return float(np.mean(vals)) if vals else np.nan

tid = adata.obs["TissueID"].values
trt = adata.obs["Treatment"].values
stt = adata.obs["state_compact"].values
drugs = sorted(set(trt) - {"vehicle"})
patients = sorted(set(tid))
print("药物（除 vehicle）:", drugs)
print("TissueID 数:", len(patients))

rows = []
for p in patients:
    for d in drugs:
        ctrl = (tid == p) & (trt == "vehicle")
        trtd = (tid == p) & (trt == d)
        if ctrl.sum() < MIN_CELLS_PATIENT_CONTROL:
            continue
        if trtd.sum() < MIN_CELLS_PATIENT_CONTROL // 2:
            continue
        kept = []
        for s in STATES:
            if ((ctrl) & (stt == s)).sum() >= MIN_CELLS_PER_STATE and ((trtd) & (stt == s)).sum() >= MIN_CELLS_PER_STATE:
                kept.append(s)
        if len(kept) < MIN_STATES_FOR_METRIC:
            continue
        C0 = np.vstack([Z_panel[ctrl & (stt == s)].mean(axis=0) for s in kept])
        C1 = np.vstack([Z_panel[trtd & (stt == s)].mean(axis=0) for s in kept])
        Delta = C1 - C0
        C = pairwise_mean_cos(Delta)
        D0, D1 = pairwise_mean_dist(C0), pairwise_mean_dist(C1)
        K = float(np.log2((D1 + EPS) / (D0 + EPS)))
        M = float(np.mean(np.linalg.norm(Delta, axis=1)))
        p0 = np.array([((ctrl) & (stt == s)).sum() / ctrl.sum() for s in STATES])
        p1 = np.array([((trtd) & (stt == s)).sum() / trtd.sum() for s in STATES])
        TV = float(0.5 * np.abs(p1 - p0).sum())
        rows.append(dict(patient=p, drug=d, n_states=len(kept), n_control=int(ctrl.sum()),
                         n_treated=int(trtd.sum()), C=C, K=K, M=M, TV=TV))
geom_b = pd.DataFrame(rows)
print("比较数:", len(geom_b))
print(geom_b.to_string())
if not SMOKE:
    assert len(geom_b) >= 40, f"比较数只有 {len(geom_b)} 个，预期 43 个左右 —— 检查 TissueID/药物/状态阈值！"
med_C = geom_b["C"].median()
n_pos = int((geom_b["C"] > 0).sum())
print()
print(f"复现结果：C 中位数 = {med_C:.3f}（作者：0.608）| C>0 的比较 = {n_pos}/{len(geom_b)}（作者：43/43）")
print(f"          M 中位数 = {geom_b['M'].median():.3f} | K 中位数 = {geom_b['K'].median():.3f}（仅报告，不做结论）")
print(f"          TV 中位数 = {geom_b['TV'].median():.3f}")
if not SMOKE:
    assert 0.40 < med_C < 0.80, "复现的 C 中位数与作者 0.608 相差太大，管线可能有问题，停下检查！"
    assert n_pos >= 40, "C>0 的比较数远低于作者的 43/43，停下检查！"
    print("B 的 C/M 复现通过 —— 管线可信，继续做轨迹检验")
else:
    print("SMOKE 模式：比较数少是预期的（1 万细胞通不过 43 比较的阈值），跳过复现断言；")
    print("          冒烟只验证管线能跑通，不验证复现数字。")
""")

code("""# --- (0) 续：每患者的中位数 C 与作者表对照 ---
merged = geom_b.groupby("patient")["C"].median().rename("C_reproduced").reset_index()
merged = merged.merge(auth_unit[["patient", "median_within_drug_C", "n_drug_pairs"]],
                      left_on="patient", right_on="patient", how="left")
merged["abs_diff"] = (merged["C_reproduced"] - merged["median_within_drug_C"]).abs()
print(merged.round(3).to_string(index=False))
if len(merged) >= 3 and merged["median_within_drug_C"].notna().all():
    r = np.corrcoef(merged["C_reproduced"], merged["median_within_drug_C"])[0, 1]
    print()
    print(f"每患者中位数 C：与作者的相关系数 r = {r:.3f}，中位绝对差 = {merged['abs_diff'].median():.3f}")
    if not SMOKE:
        assert r > 0.8, "与作者每患者中位数 C 的相关性低于 0.8，停下检查！"
        print("作者对照通过")
else:
    r = float("nan")
    print("可对照的患者太少（SMOKE 模式预期），跳过作者相关性断言")
""")

code("""# --- (a) 沿轨迹的有向位移：12 维模块空间的质心位移投影到 u_traj ---
# S：(n_cells, 12) 模块分数，列序 MM0..MM11；u_traj：12 维单位向量
traj = adata.obs["TRAJ"].values
rows = []
for _, g in geom_b.iterrows():
    p, d = g["patient"], g["drug"]
    ctrl = (tid == p) & (trt == "vehicle")
    trtd = (tid == p) & (trt == d)
    c0 = S[ctrl].mean(axis=0)
    c1 = S[trtd].mean(axis=0)
    Delta12 = c1 - c0
    d_traj = float(Delta12 @ u_traj)
    m12 = float(np.linalg.norm(Delta12))
    # 每状态的位移投影（只对 kept 状态）
    per_state = {}
    for s in STATES:
        cs0 = ctrl & (stt == s)
        cs1 = trtd & (stt == s)
        if cs0.sum() >= MIN_CELLS_PER_STATE and cs1.sum() >= MIN_CELLS_PER_STATE:
            per_state[s] = float((S[cs1].mean(axis=0) - S[cs0].mean(axis=0)) @ u_traj)
    row = dict(patient=p, drug=d, d_traj=d_traj, m12=m12,
               d_traj_AC=per_state.get("AC", np.nan), d_traj_OPC=per_state.get("OPC", np.nan),
               d_traj_NPC=per_state.get("NPC", np.nan), d_traj_MES=per_state.get("MES", np.nan))
    rows.append(row)
traj_proj = pd.DataFrame(rows)
if traj_proj.empty:
    traj_proj = pd.DataFrame(columns=["patient", "drug", "d_traj", "m12",
                                      "d_traj_AC", "d_traj_OPC", "d_traj_NPC", "d_traj_MES"])
traj_proj = traj_proj.merge(geom_b[["patient", "drug", "C", "M", "TV"]], on=["patient", "drug"])
print(traj_proj.round(3).to_string(index=False))
print()
print("d_traj（沿轨迹有向位移）：中位数", round(float(traj_proj["d_traj"].median()), 3),
      "| 正值比较数", int((traj_proj["d_traj"] > 0).sum()), "/", len(traj_proj),
      "（正 = 向晚期移动）")
print("m12（模块空间位移幅度）：中位数", round(float(traj_proj["m12"].median()), 3))
if len(traj_proj) >= 3:
    print("d_traj 与 B 的 M 的相关系数：", round(float(np.corrcoef(traj_proj["d_traj"].abs(), traj_proj["M"])[0, 1]), 3))
    print("d_traj 与 B 的 C 的相关系数：", round(float(np.corrcoef(traj_proj["d_traj"], traj_proj["C"])[0, 1]), 3))
else:
    print("可用比较太少（SMOKE 模式预期），跳过相关系数")
""")

code("""# --- (b) 轨迹分箱趋势：载体 TRAJ 五分位，每箱内算位移 ---
from scipy.stats import spearmanr

def bin_indices(values, edges):
    # edges：分位数切点（含 min/max）；返回 0..4
    b = np.searchsorted(edges[1:-1], values, side="right")
    return np.clip(b, 0, 4)

bin_rows = []
for _, g in geom_b.iterrows():
    p, d = g["patient"], g["drug"]
    ctrl = (tid == p) & (trt == "vehicle")
    trtd = (tid == p) & (trt == d)
    edges = np.quantile(traj[ctrl], [0, 0.2, 0.4, 0.6, 0.8, 1.0])
    edges = np.unique(edges)
    if len(edges) < 6:
        continue
    b_all = bin_indices(traj, edges)   # 对全部细胞分箱（切点取自载体），再与臂 mask 相交
    for b in range(5):
        c0, c1 = ctrl & (b_all == b), trtd & (b_all == b)
        if c0.sum() < 10 or c1.sum() < 10:
            continue
        disp = float(traj[c1].mean() - traj[c0].mean())          # 沿轨迹轴的位移（标量）
        mag = float(np.linalg.norm(S[c1].mean(axis=0) - S[c0].mean(axis=0)))  # 12 维位移幅度
        bin_rows.append(dict(patient=p, drug=d, quintile=b + 1, n0=int(c0.sum()), n1=int(c1.sum()),
                             disp_traj=disp, mag12=mag,
                             mean_traj_ctrl=float(traj[c0].mean())))
bins = pd.DataFrame(bin_rows)
print("分箱记录数:", len(bins), "（预期 43 比较 x 5 箱 ≈ 215，缺失的是细胞数 <10 的箱）")
if not SMOKE:
    assert len(bins) >= 150, "分箱记录太少，检查分箱逻辑！"
if len(bins) >= 10:
    summary = bins.groupby("quintile").agg(
        mean_disp=("disp_traj", "mean"), sem_disp=("disp_traj", lambda x: x.std(ddof=1) / np.sqrt(len(x))),
        mean_mag=("mag12", "mean"), sem_mag=("mag12", lambda x: x.std(ddof=1) / np.sqrt(len(x))),
        n_comp=("disp_traj", "size")).round(3)
    print()
    print("按轨迹五分位汇总（mean ± SEM）：")
    print(summary.to_string())
    rho_pool, p_pool = spearmanr(bins["quintile"], bins["disp_traj"])
    rho_mag, p_mag = spearmanr(bins["quintile"], bins["mag12"])
    rho_mean, p_mean = spearmanr(summary.index.values, summary["mean_disp"].values)
    print()
    print(f"趋势检验（位移方向）：合并数据 Spearman rho = {rho_pool:.3f} (p={p_pool:.2e})；"
          f"五箱均值 Spearman rho = {rho_mean:.3f} (p={p_mean:.3f})")
    print(f"趋势检验（位移幅度）：合并数据 Spearman rho = {rho_mag:.3f} (p={p_mag:.2e})")
    print("解读：rho 显著不为 0 且单调，支持“位移是轨迹位置的连续函数”；rho 接近 0 则不支持。")
else:
    summary = pd.DataFrame()
    rho_pool = p_pool = rho_mag = p_mag = rho_mean = p_mean = float("nan")
    print("分箱记录太少（SMOKE 模式预期），跳过趋势检验")
""")

code("""# --- (c) 幅度分歧是否收敛：同一空间里，离散状态分箱 vs 轨迹分箱的幅度离散度 ---
# 关键：两种分箱都在 12 维模块空间里算幅度，才可比。状态分箱用 kept 状态的每状态幅度；
# 轨迹分箱用上一步的每箱幅度（同一比较）。
from scipy.stats import wilcoxon

disp_rows = []
for _, g in geom_b.iterrows():
    p, d = g["patient"], g["drug"]
    ctrl = (tid == p) & (trt == "vehicle")
    trtd = (tid == p) & (trt == d)
    # 状态分箱幅度
    m_state = []
    for s in STATES:
        cs0, cs1 = ctrl & (stt == s), trtd & (stt == s)
        if cs0.sum() >= MIN_CELLS_PER_STATE and cs1.sum() >= MIN_CELLS_PER_STATE:
            m_state.append(float(np.linalg.norm(S[cs1].mean(axis=0) - S[cs0].mean(axis=0))))
    # 轨迹分箱幅度（复用 bins 表）
    sub = bins[(bins["patient"] == p) & (bins["drug"] == d)]
    m_bin = sub["mag12"].tolist()
    if len(m_state) >= 3 and len(m_bin) >= 3:
        disp_rows.append(dict(patient=p, drug=d,
                              sd_state=float(np.std(m_state, ddof=1)),
                              sd_traj=float(np.std(m_bin, ddof=1)),
                              cv_state=float(np.std(m_state, ddof=1) / np.mean(m_state)),
                              cv_traj=float(np.std(m_bin, ddof=1) / np.mean(m_bin)),
                              mean_state=float(np.mean(m_state)),
                              mean_traj=float(np.mean(m_bin))))
disp = pd.DataFrame(disp_rows)
print("可比较的比较数:", len(disp))
print(disp.round(3).to_string(index=False))
if len(disp) >= 5:
    w, pw = wilcoxon(disp["sd_state"], disp["sd_traj"])
    frac = float((disp["sd_traj"] < disp["sd_state"]).mean())
    print()
    print(f"幅度离散度（标准差）：状态分箱中位数 {disp['sd_state'].median():.3f} vs 轨迹分箱中位数 {disp['sd_traj'].median():.3f}")
    print(f"Wilcoxon 配对检验 p = {pw:.3e}；轨迹分箱离散度更小的比较占 {frac:.0%}")
    print(f"变异系数 CV：状态 {disp['cv_state'].median():.3f} vs 轨迹 {disp['cv_traj'].median():.3f}")
    print("解读：若轨迹分箱的离散度显著更小，说明 B 的“幅度分歧”部分是离散分箱的人为产物；")
    print("      若两者相当或轨迹更大，说明幅度分歧在轨迹坐标下并未收敛。")
else:
    pw = float("nan"); frac = float("nan")
    print("可比较的比较数太少（SMOKE 模式预期），跳过 Wilcoxon 检验")
""")

# ============================ Section 6 ============================
md("""## 6. 混杂拆分：状态内移动 vs 跨状态重分布

**强制要求（备忘）**：治疗后状态是**重新独立分配**的 —— 同一个细胞在治疗后可能被分到另一个状态。
所以质心位移 = **状态内移动**（细胞在原状态里沿轨迹走了多远）+ **组成变化**（各状态占比变了，质心被“换了一批细胞”拖动）。
两者必须量化拆开，**结论里分开写**，不能把位移直接读成“细胞沿轨迹移动了多远”。

**精确分解**（每比较、每分组方式）：
- total = Σ p1·C1 − Σ p0·C0（组内加权质心之差；p 为该组在臂内的占比，C 为组质心，12 维模块空间）
- within = Σ p1·(C1 − C0)（组成固定在治疗后占比，只看组内质心移动）
- composition = Σ (p1 − p0)·C0（质心固定在载体位置，只看占比变化）
- 恒等式：within + composition = total（本格用 `allclose` 逐比较验证，失败即停）。
- 再把三个向量分别投影到轨迹轴 u_traj 上，比较沿轨迹方向上两者各占多少。

**两种分组**：(i) B 的紧凑 panel 四状态（kept 状态，占比在 kept 细胞内重归一化，排除的细胞占比另列一栏）；
(ii) 轨迹五分位箱（同 §5(b) 的分箱）。
""")

code("""# --- 混杂拆分：B 状态分组 ---
import numpy as np

def decompose(groups_ctrl, groups_trtd, S, u):
    # groups_ctrl/trtd：list of bool masks（同一分组的各组，两臂对齐）。
    # 返回 (total, within, composition) 三个 12 维向量 + 沿轨迹投影。
    p0 = np.array([m.sum() for m in groups_ctrl], dtype=float)
    p1 = np.array([m.sum() for m in groups_trtd], dtype=float)
    p0 = p0 / p0.sum()
    p1 = p1 / p1.sum()
    C0 = np.vstack([S[m].mean(axis=0) if m.sum() > 0 else np.zeros(S.shape[1]) for m in groups_ctrl])
    C1 = np.vstack([S[m].mean(axis=0) if m.sum() > 0 else np.zeros(S.shape[1]) for m in groups_trtd])
    total = (p1[:, None] * C1).sum(axis=0) - (p0[:, None] * C0).sum(axis=0)
    within = (p1[:, None] * (C1 - C0)).sum(axis=0)
    composition = ((p1 - p0)[:, None] * C0).sum(axis=0)
    assert np.allclose(within + composition, total, atol=1e-6), "分解恒等式不成立，停下检查！"
    return total, within, composition, float(total @ u), float(within @ u), float(composition @ u)

dec_rows = []
for _, g in geom_b.iterrows():
    p, d = g["patient"], g["drug"]
    ctrl = (tid == p) & (trt == "vehicle")
    trtd = (tid == p) & (trt == d)
    kept = []
    for s in STATES:
        if ((ctrl) & (stt == s)).sum() >= MIN_CELLS_PER_STATE and ((trtd) & (stt == s)).sum() >= MIN_CELLS_PER_STATE:
            kept.append(s)
    if len(kept) < MIN_STATES_FOR_METRIC:
        continue
    gc = [ctrl & (stt == s) for s in kept]
    gt = [trtd & (stt == s) for s in kept]
    total, within, comp, t_t, t_w, t_c = decompose(gc, gt, S, u_traj)
    # 被排除的细胞占比（kept 之外的状态）
    excl0 = 1.0 - sum(m.sum() for m in gc) / ctrl.sum()
    excl1 = 1.0 - sum(m.sum() for m in gt) / trtd.sum()
    frac_comp = abs(t_c) / (abs(t_w) + abs(t_c) + 1e-12)
    dec_rows.append(dict(patient=p, drug=d, n_states=len(kept), M=g["M"],
                         total_traj=t_t, within_traj=t_w, composition_traj=t_c,
                         frac_composition=frac_comp,
                         excl_ctrl=excl0, excl_trtd=excl1, TV=g["TV"]))
decomp_state = pd.DataFrame(dec_rows)
print("状态分组分解（沿轨迹投影）：")
print(decomp_state.round(3).to_string(index=False))
print()
if len(decomp_state) > 0:
    med_w = float(decomp_state["within_traj"].median())
    med_c = float(decomp_state["composition_traj"].median())
    print("沿轨迹方向：within 中位数", round(med_w, 3),
          "| composition 中位数", round(med_c, 3),
          "| total 中位数", round(float(decomp_state["total_traj"].median()), 3))
    print("composition 占比（|comp|/(|within|+|comp|)）中位数：",
          round(float(decomp_state["frac_composition"].median()), 3))
    print("被排除的细胞占比中位数：载体", round(float(decomp_state["excl_ctrl"].median()), 3),
          "治疗", round(float(decomp_state["excl_trtd"].median()), 3))
    print("TV 中位数：", round(float(decomp_state["TV"].median()), 3))
else:
    med_w = med_c = float("nan")
    print("无可用比较（SMOKE 模式预期），跳过汇总")
""")

code("""# --- 混杂拆分：轨迹五分位分组（同一数学，换分组）---
dec_rows2 = []
for _, g in geom_b.iterrows():
    p, d = g["patient"], g["drug"]
    ctrl = (tid == p) & (trt == "vehicle")
    trtd = (tid == p) & (trt == d)
    edges = np.quantile(traj[ctrl], [0, 0.2, 0.4, 0.6, 0.8, 1.0])
    edges = np.unique(edges)
    if len(edges) < 6:
        continue
    b_all = bin_indices(traj, edges)   # 对全部细胞分箱（切点取自载体），再与臂 mask 相交
    gc = [ctrl & (b_all == b) for b in range(5)]
    gt = [trtd & (b_all == b) for b in range(5)]
    if min(m.sum() for m in gc + gt) < 10:
        continue
    total, within, comp, t_t, t_w, t_c = decompose(gc, gt, S, u_traj)
    frac_comp = abs(t_c) / (abs(t_w) + abs(t_c) + 1e-12)
    dec_rows2.append(dict(patient=p, drug=d, total_traj=t_t, within_traj=t_w,
                          composition_traj=t_c, frac_composition=frac_comp))
decomp_traj = pd.DataFrame(dec_rows2)
print("轨迹分箱分解（沿轨迹投影）：")
print(decomp_traj.round(3).to_string(index=False))
print()
if len(decomp_traj) > 0:
    print("沿轨迹方向：within 中位数", round(float(decomp_traj["within_traj"].median()), 3),
          "| composition 中位数", round(float(decomp_traj["composition_traj"].median()), 3))
    print("composition 占比中位数：", round(float(decomp_traj["frac_composition"].median()), 3))
else:
    print("无可用比较（SMOKE 模式预期），跳过汇总")
print()
print("结论写法（必须分开写）：")
print("  1) 状态内移动（within）沿轨迹方向的中位数 = X —— 这才是“细胞在原状态里沿轨迹移动”的部分。")
print("  2) 组成变化（composition）沿轨迹方向的中位数 = Y —— 这是“各状态占比变化拖动质心”的部分。")
print("  3) 若 composition 占主导，则 §5 的位移主要反映重分布，不能读成细胞沿轨迹移动；")
print("     若 within 占主导且 (b) 的趋势显著，则支持幅度是轨迹位置的连续函数。")
""")

# ============================ Section 7 ============================
md("""## 7. 输出与 verdict

**目的**：把图与结果表存盘，最后写一段 verdict —— 幅度是不是轨迹位置的连续函数？
**判断标准**（填 verdict 前逐项核对）：
- §5(0) 复现通过（C 中位数 ≈ 0.608、43/43 为正、每患者中位数与作者 r > 0.8）；
- §5(b) 趋势检验：合并数据 Spearman 是否显著（p < 0.05）且方向单调；
- §5(c) 轨迹分箱的幅度离散度是否显著小于状态分箱；
- §6 分解：沿轨迹的位移里 within 占多少、composition 占多少 —— composition 占主导时，§5 的位移不能读成细胞沿轨迹移动。
""")

code("""# --- 图 1：轨迹位置 vs 位移（五分位趋势）---
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

if len(bins) >= 10:
    plt.figure(figsize=(8, 5))
    for (p, d), sub in bins.groupby(["patient", "drug"]):
        sub = sub.sort_values("quintile")
        plt.plot(sub["quintile"], sub["disp_traj"], color="grey", alpha=0.25, linewidth=0.7)
    plt.errorbar(summary.index.values, summary["mean_disp"].values, yerr=summary["sem_disp"].values,
                 marker="o", color="tab:blue", linewidth=2, capsize=4, label="mean ± SEM (43 比较)")
    plt.axhline(0, color="black", linewidth=0.8)
    plt.xlabel("轨迹位置五分位（1 = 最早期，5 = 最晚期）")
    plt.ylabel("沿轨迹的有向位移（治疗 − 载体）")
    plt.title(f"位移随轨迹位置的变化\\nSpearman rho = {rho_pool:.2f} (p={p_pool:.1e})，五箱均值 rho = {rho_mean:.2f}")
    plt.xticks([1, 2, 3, 4, 5])
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"{OUT_DIR}/fig1_traj_trend.png", dpi=150)
    plt.close()
    print("已存 fig1_traj_trend.png")
else:
    print("跳过图 1：分箱记录太少（SMOKE 模式预期）")

# --- 图 2：B 的 C/M 复现对照 ---
if len(geom_b) >= 3 and len(merged) >= 3 and merged["median_within_drug_C"].notna().all():
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    axes[0].scatter(merged["median_within_drug_C"], merged["C_reproduced"], s=60)
    for _, r_ in merged.iterrows():
        axes[0].annotate(r_["patient"], (r_["median_within_drug_C"], r_["C_reproduced"]), fontsize=7)
    lims = [min(merged["median_within_drug_C"].min(), merged["C_reproduced"].min()) - 0.05,
            max(merged["median_within_drug_C"].max(), merged["C_reproduced"].max()) + 0.05]
    axes[0].plot(lims, lims, "k--", linewidth=0.8)
    axes[0].set_xlim(lims); axes[0].set_ylim(lims)
    axes[0].set_xlabel("作者：每患者中位数 C")
    axes[0].set_ylabel("本 notebook 复现：每患者中位数 C")
    axes[0].set_title(f"作者对照 (r = {r:.2f})")
    axes[1].hist(geom_b["C"], bins=12, edgecolor="black")
    axes[1].axvline(med_C, color="tab:red", linewidth=2, label=f"中位数 {med_C:.3f}")
    axes[1].axvline(0.608, color="tab:green", linestyle="--", linewidth=2, label="作者 0.608")
    axes[1].set_xlabel("C（方向一致性），每比较")
    axes[1].set_ylabel("比较数")
    axes[1].set_title(f"C 分布（{n_pos}/{len(geom_b)} 为正）")
    axes[1].legend()
    plt.tight_layout()
    plt.savefig(f"{OUT_DIR}/fig2_reproduce.png", dpi=150)
    plt.close()
    print("已存 fig2_reproduce.png")
else:
    print("跳过图 2：可用比较太少（SMOKE 模式预期）")

# --- 图 3：混杂拆分（状态分组）---
if len(decomp_state) >= 1:
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    axes[0].bar(["within\\n(状态内移动)", "composition\\n(组成变化)"], [med_w, med_c],
                color=["tab:blue", "tab:orange"])
    axes[0].axhline(0, color="black", linewidth=0.8)
    axes[0].set_ylabel("沿轨迹的有向位移（中位数，43 比较）")
    axes[0].set_title("质心位移的拆分（状态分组）")
    axes[1].scatter(decomp_state["M"], decomp_state["frac_composition"], s=50)
    axes[1].set_xlabel("M（B 的幅度）")
    axes[1].set_ylabel("composition 占比 |comp|/(|within|+|comp|)")
    axes[1].set_title("组成变化占比 vs 响应幅度")
    plt.tight_layout()
    plt.savefig(f"{OUT_DIR}/fig3_decomp.png", dpi=150)
    plt.close()
    print("已存 fig3_decomp.png")
else:
    print("跳过图 3：无可用比较（SMOKE 模式预期）")

# --- 图 4：B 状态 x A 模块热图 ---
mm_cols = [c for c in adata.obs.columns if c.startswith("MM")]
means = adata.obs.groupby("state_compact", observed=True)[mm_cols].mean()
zmeans = (means - means.mean(axis=0)) / means.std(axis=0)
fig, ax = plt.subplots(figsize=(10, 4))
im = ax.imshow(zmeans.values, aspect="auto", cmap="RdBu_r", vmin=-1.5, vmax=1.5)
ax.set_xticks(range(12))
ax.set_xticklabels([f"MM{m}" for m in range(12)], rotation=45, ha="right")
ax.set_yticks(range(4))
ax.set_yticklabels(STATES)
ax.set_title("每 B 状态的平均 S10 模块分（按列 z-score）")
for i in range(4):
    for j in range(12):
        ax.text(j, i, f"{zmeans.values[i, j]:.1f}", ha="center", va="center", fontsize=7)
plt.colorbar(im, ax=ax, shrink=0.8)
plt.tight_layout()
plt.savefig(f"{OUT_DIR}/fig4_state_module.png", dpi=150)
plt.close()
print("已存 fig4_state_module.png")
print()
print("图已存到", OUT_DIR)
""")


code("""# --- 结果表存 CSV ---
geom_b.to_csv(f"{OUT_DIR}/geom_b_reproduction.csv", index=False)
traj_proj.to_csv(f"{OUT_DIR}/traj_projection.csv", index=False)
bins.to_csv(f"{OUT_DIR}/quintile_bins.csv", index=False)
summary.to_csv(f"{OUT_DIR}/quintile_summary.csv")
disp.to_csv(f"{OUT_DIR}/dispersion_state_vs_traj.csv", index=False)
decomp_state.to_csv(f"{OUT_DIR}/decomposition_state_groups.csv", index=False)
decomp_traj.to_csv(f"{OUT_DIR}/decomposition_traj_groups.csv", index=False)
means.to_csv(f"{OUT_DIR}/state_module_means.csv")
print("已存 CSV：")
for name in ["geom_b_reproduction.csv", "traj_projection.csv", "quintile_bins.csv",
             "quintile_summary.csv", "dispersion_state_vs_traj.csv",
             "decomposition_state_groups.csv", "decomposition_traj_groups.csv",
             "state_module_means.csv"]:
    print("  ", name)
print()
print("=== 关键数字汇总（填 verdict 用）===")
def _med(df, col):
    return float(df[col].median()) if (len(df) > 0 and col in df.columns) else float("nan")
def _fmt(x, fmt="{:.3f}"):
    return "—" if x != x else fmt.format(x)   # nan 时显示破折号
def _pct(x):
    return "—" if x != x else f"{x:.0%}"
print(f"复现 C 中位数：{med_C:.3f}（作者 0.608）；C>0：{n_pos}/{len(geom_b)}")
print(f"每患者中位数 C 与作者 r = {r:.3f}")
print(f"(b) 位移方向趋势：合并 Spearman rho = {rho_pool:.3f} (p={p_pool:.2e})；五箱均值 rho = {rho_mean:.3f} (p={p_mean:.3f})")
print(f"(b) 位移幅度趋势：合并 Spearman rho = {rho_mag:.3f} (p={p_mag:.2e})")
print(f"(c) 幅度离散度：状态 {_fmt(_med(disp, 'sd_state'))} vs 轨迹 {_fmt(_med(disp, 'sd_traj'))}，Wilcoxon p = {_fmt(pw, '{:.3e}')}，轨迹更小占 {_pct(frac)}")
print(f"(6) 状态分组：within 中位数 {med_w:.3f}，composition 中位数 {med_c:.3f}，composition 占比中位数 {_med(decomp_state, 'frac_composition'):.2f}")
print(f"(6) 轨迹分箱：within 中位数 {_med(decomp_traj, 'within_traj'):.3f}，composition 中位数 {_med(decomp_traj, 'composition_traj'):.3f}，占比中位数 {_med(decomp_traj, 'frac_composition'):.2f}")
""")

md("""## Verdict（跑完后填）

**问题：药物响应幅度是不是轨迹位置的连续函数？**

- **结论（圈一个）：是 / 否 / 不确定**
- **证据**：
  1. 复现对照：C 中位数 ___（作者 0.608），C>0 ___/43，每患者中位数与作者 r = ___ —— 管线可信 / 不可信。
  2. 轨迹分箱趋势：位移方向合并 Spearman rho = ___（p = ___）；位移幅度 rho = ___（p = ___）。
     五箱均值从第 1 箱 ___ 到第 5 箱 ___，单调 / 不单调。
  3. 幅度分歧收敛：状态分箱离散度中位数 ___ vs 轨迹分箱 ___，Wilcoxon p = ___，轨迹更小的比较占 ___%。
  4. 混杂拆分：沿轨迹的位移里，状态内移动（within）中位数 ___，组成变化（composition）中位数 ___，
     composition 占比中位数 ___。—— 若 composition 占主导，位移主要反映重分布，不能读成细胞沿轨迹移动。
- **对虚拟肿瘤药物模块的含义（圈一个）**：
  - **连续轨迹输入**：若趋势显著、轨迹分箱离散度更小、且 within 占主导 —— 药物模块应以连续轨迹位置（TRAJ）为输入，
    响应幅度 = f(轨迹位置)，而不是离散状态。
  - **离散状态分型**：若趋势不显著或轨迹分箱并不更优 —— 保留 B 的离散状态分型，轨迹只作解释性坐标。
  - **不确定，需更多数据**：若结果混合（如趋势显著但 composition 占主导）—— 先做 §6 要求的补充分析（按药物分层、按 TissueID 分层），
    再决定；不要强行下结论。
- **红线重申**：以上结论只针对转录响应（C/M），不是疗效。任何“阻断轨迹能治病”的说法都需要功能与体内验证。

---

### 在 Colab 上跑的具体步骤（给用户）
1. 打开 https://colab.research.google.com/ → 上传本 `.ipynb` 文件（或从 Google Drive 打开）。
2. 运行时 → 更改运行时类型 → **CPU** 即可（不需要 GPU；AI Ultra 的高内存运行时更稳；若选 GPU 也能跑，但不是必需）。
3. 运行时 → **全部运行**（Run all），按顺序执行；每个下载/加载步骤后都有断言，失败会停下报错 —— 不要跳过报错继续跑。
4. 预计 60–100 分钟（主要是两次全矩阵流式读：算文库大小 + 读选定基因；重跑读缓存会快很多）。输出在 `/content/point1_outputs/`（4 张 PNG + 8 张 CSV），跑完后从左侧文件面板下载。
5. **冒烟测试**：先把顶部配置里的 `SMOKE = True` 跑一遍（只取 1 万 Tumor 细胞，约 5 分钟），确认管线通了再把 `SMOKE` 改回 `False` 跑全量。
   注意：冒烟模式的数字只是管线验证，**不是结论**。
6. 跑完后把上面 Verdict 里的空填上，再决定药物模块走连续轨迹还是离散状态。

### 显式假设（需要用户拍板）
- **轨迹位置的 late-vs-early 定义**（§4）：早期 = MM1/4/6/7/8/10/11，晚期 = MM2/3/5，MM0（Cycling）与 MM9（AC-progenitor-like）主定义排除。
  这是按模块注释 + factor_id_map 做的**人为分组**，不是论文直接给出的标签。敏感性分析（MM9 归入早期、PC1 轴）已内置；
  若用户认为 MM9 应属晚期、或 MM0 不应排除，改 `EARLY` / `LATE` 两个列表重跑即可。
- **B 状态用紧凑 panel 的 argmax**（§3）：与作者一致；margin 阈值 0.25 只用于报告，不用于过滤（主分析用全部细胞，与作者主分析一致）。
- **轨迹分箱的切点取自载体细胞**（§5b/§6）：治疗细胞按同一套切点分箱，这是“同一坐标系下比较”的前提。
- **只建模 C 与 M，不建模 K；不把转录响应读成疗效**（两条红线，全程遵守）。
""")

# ============================ Assemble & write ============================
nb.cells = cells
with open(OUT_PATH, "w") as f:
    nbf.write(nb, f)
print(f"Wrote {OUT_PATH} with {len(cells)} cells")
