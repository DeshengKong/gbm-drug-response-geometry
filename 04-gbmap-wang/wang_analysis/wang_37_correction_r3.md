# Wang (GBmap) pair #37 correction + #39 SWITCH — _r3 recompute memo

**Date:** 2026-11-09 (runs 05:36–05:47 UTC / 13:36–13:47 CST; verification 13:48–14:00 CST)
**Trigger:** Author Aaron Diaz replied by email 2026-10-09 (message 1a11d8bd6c7f87b7,
thread "GSE174554 reuse"): *"The correct assignment, based on our supplementary table and
internal tissue records, is SF9798 as primary and SF9494 as recurrent. The GEO progression
annotations appear to have been reversed during metadata preparation."* (Diaz calls it Pair 27;
our numbering is #37.) The user verified the specimen–GSM correspondence against NCBI official
filenames.
**Decision timeline for #39:** 2026-11-09 08:37 user decided KEEP GSM5319574 (SF12754) +
disclose → **2026-11-09 13:47 REVERSED: SWITCH. GSM5319563 (SF12407, FIRST recurrence) is now
the PRIMARY analysis; GSM5319574 (SF12754, SECOND recurrence) is the SENSITIVITY variant.**
This memo frames everything around the SF12407 base. It supersedes the Phase B spec's KEEP
branch and the _s39 sensitivity memo's roles.

## 0. Branch / file naming — read this first

The _r3 scripts were written while KEEP was still the decision, so the file suffixes encode
the OLD roles:

| Branch | Recurrent sample | Files | Role NOW |
|---|---|---|---|
| _s39_r3 | GSM5319563 (SF12407, first recurrence) | `wang_per_sample_s39_r3.csv`, `wang_ap1_s39_r3.csv` | **BASE / PRIMARY** |
| _r3 (main run) | GSM5319574 (SF12754, second recurrence) | `wang_per_sample_r3.csv`, `wang_ap1_r3.csv` | **SENSITIVITY** |

Both branches carry the #37 correction. Nothing in this memo renames files; downstream use
(DOCX, Fig 7) must take numbers from the **_s39_r3** files.

## 1. GSM↔specimen mapping verification (done BEFORE any change)

Verified from the on-disk GEO metadata (SOFT `GSE174554_family.soft.gz`, 113 samples parsed;
NCBI official filenames `GSM5319532_SF9494_*`, `GSM5319559_SF9798_*`; SOFT titles agree):

- GSM5319532 = SF9494, GSM5319559 = SF9798. Both age 64 (same patient). The metadata does
  NOT contradict the mapping.
- GEO progression annotations say SF9494 = Primary, SF9798 = Recurrent — exactly what Diaz
  says is reversed. Correction: **#37 = ("GSM5319559", "GSM5319532")**, i.e. primary SF9798 /
  recurrent SF9494.
- Old #37 rows (pre-correction CSVs): GSM5319532 Primary (n=3457, 82.9% tumor, TRAJ median
  tumor −0.2201); GSM5319559 Recurrent (n=4051, 96.1% tumor, TRAJ −0.1153).
  Old pair diff (recurrent − primary) = **+0.1048 (up)**; corrected = **−0.1048 (down)** — sign flip.
- #37 AP-1 core diff likewise flips: −0.1112 → **+0.1112**.

## 2. Cache handling — full regeneration, no stale labels

`wang_analysis.py:91-93` (and the _s39 twin) reuse rows from the existing output CSV via a
`done_gsm` set, so per-sample rows carry pair/timepoint labels baked in under the old
assignment. Per the user's flag:

- The four cached CSVs (`wang_per_sample.csv`, `wang_per_sample_s39.csv`, `wang_ap1.csv`,
  `wang_ap1_s39_single.csv`) were **moved aside** to
  `~/workspace/bioinfo_downloads/wang2022/wang_analysis/backup_pre_r3_20261109/` (never deleted;
  repo copies untouched, md5s re-verified).
- Four new scripts, diff-verified to change ONLY the intended lines, all `py_compile` OK:
  `wang_analysis_r3.py`, `wang_analysis_s39_r3.py`, `wang_ap1_r3.py`,
  `wang_ap1_s39_r3.py` (the _s39_r3 AP-1 script is a FULL fresh run from `wang_ap1.py`,
  not the old single-sample + merge pattern). All four carry
  `"#37": ("GSM5319559", "GSM5319532")`.
- Driver `/tmp/run_r3.sh` ran all four sequentially: **all exit=0, total ~11.5 min
  (05:36–05:47 UTC)**. No Traceback in any log.
- **Regeneration proof** (independent verification script `compute_stats_r3.py`, see §4):
  - Main CSVs (both branches): **ZERO numeric columns changed** vs the backup old files —
    only the two #37 rows' labels flipped. 58 rows old vs 58 rows new.
  - Swap AP-1 CSV: 8 numeric columns show max |Δ| ≈ 1.1e-16 — pure floating-point noise
    from the full fresh re-run (old _s39 AP-1 was single-sample + merge); no statistic affected.
  - Label check on all four _r3 CSVs: GSM5319532 = Recurrent, GSM5319559 = Primary in every
    file; **0 occurrences** of stale labels (`GSM5319532,Primary` or `GSM5319559,Recurrent`).
  - The 57 shared samples (all except the #39 recurrent pick) reproduce the submitted CSV
    bit-identically in the main runs, as in the _s39 run.
- Pair table: `wang_33pair_table_r3.csv` — diff vs _r2 shows ONLY #37 and #39 rows changed.
  #37 row now lists primary GSM5319559 / recurrent GSM5319532 with a note recording both the
  Phase A snRNA-vs-snATAC correction and the Diaz primary/recurrent swap. #39 row updated
  2026-11-09 to record the SWITCH decision and the file mapping (§0).

## 3. Diaz's other points — checked against our mapping

- **Pair 36 = our #39 (SF11873 primary; SF12407 FIRST recurrence; SF12754 SECOND
  recurrence).** CONFIRMED by the same email and by SOFT/sample records. This is the direct
  basis for the SWITCH: a primary–first-recurrence pairing uses SF12407 (GSM5319563).
  SF12754 is the second recurrence and now the sensitivity variant.
- **SF9259R/S (separate snRNA preps of the same primary specimen SF9259; recurrent SF9510).**
  SF9259 appears **0 times** in the SOFT file. Pair #6 has only recurrent SF9510 (snRNA
  GSM5319558 + snATAC GSM5319582, age 73) — no primary sample with pair# #6 exists in GEO.
  Pair #6 is a recurrent-only singleton (Phase A memo already lists #4–#7 as recurrent-only
  singletons, never candidates); it is not in our 29 and never could be formed from GEO data.
  No action.
- **SF4297/SF6621 IDH (Diaz recommends provisionally excluding that pair for a strictly
  IDH-wildtype analysis).** This is our **#8**: SF4297 primary (snATAC only, GSM5319580) /
  SF6621 recurrent (snRNA GSM5319554 + snATAC GSM5319579). **Already excluded from our 29**
  (primary has no snRNA) — the recommendation is MOOT for our analysis. Our 29-pair analysis
  applies no IDH filter at all.
- 113 samples parsed from SOFT; 6 samples lack a pair# field.

## 4. Results — full delta table (old _r2/_s39 vs new _r3, #37 corrected)

Methods identical to the _r2/_s39 runs: Wilcoxon signed-rank (scipy defaults, exact for
n=29); permutation = 10,000 within-pair label flips, `np.random.default_rng(0)`; ligands =
Spearman vs tumor TRAJ; Bonferroni threshold 0.05/12 = 0.00417. Independent verification by
`compute_stats_r3.py` recomputing every headline statistic from the per-sample CSVs
(old = backup_pre_r3_20261109, new = _r3 files). All values below are the verification
script's output, rounded as shown.

**Note on the permutation null:** the null 95% CI is IDENTICAL old vs new in both branches
(base (−0.0699, +0.0699); sensitivity (−0.0825, +0.0825)) — the null depends only on the
multiset of |pair diffs|, which a sign flip of #37 does not change. Only the observed
statistic moved.

### 4A. BASE — SF12407 (GSM5319563), files _s39_r3; old = _s39 swap

| Statistic | Old (_s39) | New (_s39_r3) | Change |
|---|---|---|---|
| TRAJ median diff (recurrent − primary) | −0.0533 | **−0.0699** | stronger negative |
| up / down (of 29) | 9 / 20 | **8 / 21** | one more down |
| Wilcoxon p | 0.01788 | **0.007456** | stronger |
| Permutation p | 0.1276 | **0.0724** | closer to 0.05, still NS |
| Perm null 95% CI | (−0.0699, +0.0699) | (−0.0699, +0.0699) | unchanged |
| Mean-TRAJ median diff / Wilcoxon p | −0.0558 / 0.01294 | **−0.0656 / 0.006452** | stronger |
| #37 pair diff | +0.1048 | **−0.1048** | sign flip |
| T-cell median diff / up-down / Wilcoxon p | 0.0 / 12-14 / 0.6752 | 0.0 / 13-13 / **0.8291** | still null |
| VEGFA rho / p (n=49) | 0.413 / 0.003147 | 0.413 / 0.003147 | unchanged |
| VCAN rho / p (n=49) | 0.489 / 0.0003627 | 0.489 / 0.0003627 | unchanged |
| ANXA1 rho / p (n=49) | 0.205 / 0.1576 | 0.205 / 0.1576 | unchanged |
| AP-1 core median diff / up-down / p | +0.0316 / 17-12 / 0.1976 | +0.0380 / 18-11 / **0.07599** | still null |
| AP-1 ext median diff / up-down / p | +0.0140 / 17-12 / 0.2383 | +0.0215 / 18-11 / **0.09633** | still null |
| AP-1 core_all median diff / up-down / p | +0.0387 / 17-12 / 0.09633 | +0.0440 / 18-11 / **0.0308** | NOW nominally significant — see §5 |
| FOS (tumor) / p | +0.0194 / 0.6391 | +0.0369 / **0.3357** | null |
| JUN (tumor) / p | +0.0267 / 0.3145 | +0.0492 / **0.2383** | null |
| FOSL1 (tumor) / p | −0.0070 / 0.01294 | −0.0066 / **0.05338** | LOST nominal significance |
| ATF3 (tumor) / p | +0.0043 / 0.8314 | +0.0077 / **0.4294** | null |
| Primary AP-1 vs TRAJ rho / p | 0.101 / 0.6022 | 0.193 / **0.3156** | null |
| #37 AP-1 core diff | −0.1112 | **+0.1112** | sign flip |

Per-module (median diff / Wilcoxon p; BONF = survives Bonferroni 0.00417):

| Module | Old (_s39) | New (_s39_r3) |
|---|---|---|
| MM0 Cycling | −0.0158 / 0.3357 ns | −0.0178 / 0.1316 ns |
| MM1 NPC-neuronal-like 1 | +0.0465 / 0.002759 BONF | +0.0465 / 0.001527 BONF |
| MM2 Gliosis-hypoxia | −0.0222 / 0.1316 ns | −0.0222 / 0.1373 ns |
| MM3 AC-gliosis-like 2 | −0.0151 / 0.3692 ns | −0.0154 / 0.2214 ns |
| MM4 NPC-OPC-like 2 | +0.0364 / 0.01211 nom | **+0.0520 / 0.00234 BONF** |
| MM5 AC-gliosis-like 1 | −0.0323 / 0.1104 ns | −0.0323 / 0.1373 ns |
| MM6 NPC-neuronal-like 2 | +0.0865 / 0.00234 BONF | +0.0867 / 0.0004988 BONF |
| MM7 NPC-OPC-like 1 | +0.0373 / 0.09633 ns | +0.0311 / 0.1689 ns |
| MM8 AC-NPC-OPC-like 1 | +0.0349 / 0.002759 BONF | **+0.0304 / 0.01058 nom** |
| MM9 AC-progenitor-like | +0.0070 / 0.8815 ns | +0.0004 / 0.8314 ns |
| MM10 NPC-neuronal-like 3 | +0.0398 / 0.03863 nom | +0.0421 / 0.007456 nom |
| MM11 OPC-like | +0.0947 / 0.006452 nom | **+0.1114 / 0.001816 BONF** |

Bonferroni survivors, BASE: {MM1, MM6, MM8} → **{MM1, MM4, MM6, MM11}**. Nominally-significant
early modules still 6 (MM1, MM4, MM6, MM8, MM10, MM11).

### 4B. SENSITIVITY — SF12754 (GSM5319574), files _r3; old = original (KEEP)

| Statistic | Old (orig.) | New (_r3) | Change |
|---|---|---|---|
| TRAJ median diff | −0.0533 | **−0.0699** | stronger negative |
| up / down (of 29) | 10 / 19 | **9 / 20** | one more down |
| Wilcoxon p | 0.03453 | **0.01679** | stronger |
| Permutation p | 0.2023 | **0.1225** | still NS |
| Perm null 95% CI | (−0.0825, +0.0825) | (−0.0825, +0.0825) | unchanged |
| Mean-TRAJ median diff / Wilcoxon p | −0.0558 / 0.02433 | **−0.0656 / 0.01383** | stronger |
| #37 pair diff | +0.1048 | **−0.1048** | sign flip |
| T-cell median diff / up-down / Wilcoxon p | 0.0 / 12-14 / 0.8489 | 0.0 / 13-13 / **0.6752** | still null |
| VEGFA rho / p (n=49) | 0.408 / 0.003598 | 0.408 / 0.003598 | unchanged |
| VCAN rho / p (n=49) | 0.470 / 0.0006633 | 0.470 / 0.0006633 | unchanged |
| ANXA1 rho / p (n=49) | 0.199 / 0.1699 | 0.199 / 0.1699 | unchanged |
| AP-1 core median diff / up-down / p | +0.0316 / 17-12 / 0.1901 | +0.0380 / 18-11 / **0.07599** | still null |
| AP-1 ext median diff / up-down / p | +0.0140 / 17-12 / 0.2132 | +0.0215 / 18-11 / **0.08774** | still null |
| AP-1 core_all median diff / up-down / p | +0.0387 / 17-12 / 0.09633 | +0.0440 / 18-11 / **0.0308** | NOW nominally significant — see §5 |
| FOS (tumor) / p | +0.0194 / 0.6545 | +0.0369 / **0.3467** | null |
| JUN (tumor) / p | +0.0492 / 0.2214 | +0.0806 / **0.1689** | null |
| FOSL1 (tumor) / p | −0.0070 / 0.03453 | −0.0066 / **0.1104** | LOST nominal significance |
| ATF3 (tumor) / p | +0.0043 / 0.7819 | +0.0077 / **0.3927** | null |
| Primary AP-1 vs TRAJ rho / p | 0.101 / 0.6022 | 0.193 / **0.3156** | null |
| #37 AP-1 core diff | −0.1112 | **+0.1112** | sign flip |

Per-module, SENSITIVITY (median diff / Wilcoxon p):

| Module | Old (orig.) | New (_r3) |
|---|---|---|
| MM0 Cycling | −0.0158 / 0.3467 ns | −0.0178 / 0.1373 ns |
| MM1 NPC-neuronal-like 1 | +0.0465 / 0.002759 BONF | +0.0465 / 0.001527 BONF |
| MM2 Gliosis-hypoxia | −0.0216 / 0.2132 ns | −0.0216 / 0.2214 ns |
| MM3 AC-gliosis-like 2 | −0.0151 / 0.4946 ns | −0.0154 / 0.325 ns |
| MM4 NPC-OPC-like 2 | +0.0364 / 0.03654 nom | +0.0520 / 0.009877 nom |
| MM5 AC-gliosis-like 1 | −0.0323 / 0.1976 ns | −0.0323 / 0.2297 ns |
| MM6 NPC-neuronal-like 2 | +0.0865 / 0.009877 nom | **+0.0867 / 0.002992 BONF** |
| MM7 NPC-OPC-like 1 | +0.0373 / 0.1056 ns | +0.0311 / 0.1901 ns |
| MM8 AC-NPC-OPC-like 1 | +0.0349 / 0.00234 BONF | **+0.0304 / 0.009216 nom** |
| MM9 AC-progenitor-like | +0.0070 / 0.7983 ns | +0.0004 / 0.9321 ns |
| MM10 NPC-neuronal-like 3 | +0.0398 / 0.04083 nom | +0.0421 / 0.008008 nom |
| MM11 OPC-like | +0.0750 / 0.009877 nom | **+0.0947 / 0.002759 BONF** |

Bonferroni survivors, SENSITIVITY: {MM1, MM8} → **{MM1, MM6, MM11}**. Nominally-significant
early modules still 6 (MM1, MM4, MM6, MM8, MM10, MM11).

## 5. Which conclusions change (framed on the SF12407 BASE)

- **Direction (Wang paired TRAJ) STANDS and is slightly stronger.** Median diff −0.053 →
  −0.070, Wilcoxon 0.018 → 0.0075, up/down 9/20 → 8/21. BUT the permutation test is still
  not significant (0.0724, was 0.1276) — the "marginal / edge" framing must be KEPT. The
  sensitivity branch agrees (perm p 0.1225). No paper conclusion changes on direction, but
  the numbers get stronger in the classical test and remain non-significant in the
  permutation test; both must be reported together.
- **Module Bonferroni composition changes.** MM8 drops out in BOTH branches (was a survivor
  in both). BASE survivors are now {MM1, MM4, MM6, MM11}; sensitivity {MM1, MM6, MM11}.
  The count of nominally-significant early modules is unchanged (6 in both branches), so the
  qualitative "early modules rise at recurrence" statement survives — only the Bonferroni
  survivor list changes.
- **FOSL1 nominal-down finding is GONE in both branches** (BASE p 0.0534, just above 0.05;
  sensitivity p 0.1104; was 0.0129 / 0.0345). Any "FOSL1 lower" wording (Fig 7 panel D
  annotation, DOCX) must be removed.
- **AP-1 core still null** (BASE p 0.0760, was 0.1976; sensitivity same 0.0760). No change.
- **NEW WRINKLE — AP-1 core_all now nominally significant in both branches** (p 0.0308, was
  0.0963). This is the all-cells AP-1 aggregate, NOT the pre-specified tumor-cell core, and
  it was not a headline claim. Flag for the author/parent decision: do NOT silently promote
  it to a finding; options are (a) leave it unreported, (b) report as an exploratory nominal
  result with the caveat. No conclusion changes by itself.
- **T-cell still null in both branches** (BASE p 0.829, sensitivity p 0.675). No change.
- **Ligands unchanged exactly** (VEGFA, VCAN, ANXA1 identical in both branches — per-sample
  values are unaffected by the #37 label flip, which only reassigns which sample is
  "primary" within pair #37; the ligand–TRAJ Spearman pools all samples). No change.
- **Primary AP-1 vs TRAJ still null** (rho 0.193, p 0.316; was 0.101/0.60). No change.

## 6. Fig 7 — re-render REQUIRED (panels C, C-right, D)

`fig7()` in `plot_figures.py` (lines 977–1130) loads `wang_per_sample.csv` + `wang_ap1.csv`.
Panels showing pair #37's paired data:

- **Panel C (Wang paired TRAJ lines):** #37's line flips (was up, now down). Hardcoded
  annotation "Median diff −0.053 (19 down / 10 up), Wilcoxon p = 0.035, permutation p = 0.20"
  must become **"Median diff −0.070 (21 down / 8 up), Wilcoxon p = 0.0075, permutation
  p = 0.072"** (BASE values).
- **Panel C-right (per-module median-diff bars):** values shift slightly (MM4/MM11 bars grow,
  MM8 bar shrinks; Bonferroni markers move: MM8 loses its marker, MM4 and MM11 gain one).
- **Panel D (AP-1 paired lines):** #37's line flips. Hardcoded annotation "Median diff +0.032
  (17/29 up, p = 0.19), FOSL1 lower (p = 0.035), AP-1 vs TRAJ ρ = 0.10 (p = 0.60)" must
  become **"Median diff +0.038 (18/29 up, p = 0.076), FOSL1 not significant (p = 0.053),
  AP-1 vs TRAJ ρ = 0.19 (p = 0.32)"**.
- Panels A, B (GBmap data) and the D-right ligand–TRAJ panel (per-sample scatter, labels
  irrelevant) are UNAFFECTED in the base variant.
- Figures were NOT re-rendered in this task (per instructions). Re-render must use the
  **_s39_r3** files as the base data source.

## 7. Files needing a correction note (DOCX numbers to replace — SF12407 BASE)

The revised DOCX (Phase B, md5 3811f850…) currently carries the KEEP branch + wrong #37
order. Its Wang numbers must be replaced with the BASE values below. Do NOT edit the DOCX
in this task; this list is for the next revision round.

| DOCX item (current, stale) | Replacement (SF12407 BASE, #37 corrected) |
|---|---|
| TRAJ median diff −0.053 | **−0.0699 (−0.070)** |
| up/down 10 up / 19 down | **8 up / 21 down** |
| Wilcoxon p = 0.035 | **p = 0.0075** |
| Permutation p = 0.20 | **p = 0.072** (still NS — keep "marginal" framing) |
| Bonferroni survivors MM1, MM8 | **MM1, MM4, MM6, MM11** |
| FOSL1 lower (p = 0.035) | **REMOVE — FOSL1 p = 0.053, not significant** |
| AP-1 core p = 0.19 | **p = 0.076** (still null) |
| AP-1 vs TRAJ ρ = 0.10 (p = 0.60) | **ρ = 0.19 (p = 0.32)** |
| Fig 7 panel C annotation (see §6) | re-render with BASE values |
| Fig 7 panel D annotation (see §6) | re-render with BASE values |
| #39 description (SF12754 as the paired recurrence) | **SF12407 (GSM5319563) is the paired first recurrence; SF12754 is the sensitivity variant** |

Other files needing a correction note:
- **Phase A memo** (`reanalysis-results-2026-11-09.md`, Wang section): its Wang numbers are the
  KEEP + wrong-#37 values; add a correction note pointing here (do not rewrite history —
  the Phase A numbers were correct for the mapping then assumed).
- **`wang_39_sensitivity_r2.md`**: roles are now FLIPPED — what it calls the "swap" is the
  primary analysis and what it calls the base/KEEP is the sensitivity. Needs a header
  correction note pointing here; its numeric tables remain valid as the sensitivity branch
  (they match §4B up to the #37 correction — see the delta table).
- **Revision package** (`revision-package-2026-11-09/`): _r2 files must NOT be overwritten.
  The _r3 files are a separate, later deposit decision (parent's call); GitHub deposit of
  _r2/_s39 is still blocked on the user's "现在推" + fresh 30-min approval.
- **Phase B spec** (`revision-spec-phaseB-2026-11-09.md`): KEEP branch superseded by the
  13:47 SWITCH; record the supersession, do not silently edit the spec.

## 8. File inventory (_r3)

All in `~/workspace/gbm-drug-geometry-code/04-gbmap-wang/wang_analysis/` (repo copies;
originals also in `~/workspace/bioinfo_downloads/wang2022/wang_analysis/`):

| File | md5 | Notes |
|---|---|---|
| `wang_per_sample_r3.csv` (59 lines) | 587037170ce31c734f92fc25e13ff526 | SENSITIVITY branch (SF12754) |
| `wang_per_sample_s39_r3.csv` (59 lines) | 2e8a39d5386f9463d3cbfe054cacfd64 | **BASE branch (SF12407)** |
| `wang_ap1_r3.csv` (59 lines) | edc241431b867b209cdc597448f0f283 | SENSITIVITY branch (SF12754) |
| `wang_ap1_s39_r3.csv` (59 lines) | b8ddeb8c55774b8d5d97b0f8c4be461d | **BASE branch (SF12407)** |
| `wang_33pair_table_r3.csv` | 98d12e07fe82148039cb855867b8229e | #37 corrected; #39 note updated to SWITCH |
| `wang_analysis_r3.log` | 6a7033a0ae676aa673d8e75e68bbcc86 | main run, exit 0 |
| `wang_analysis_s39_r3.log` | 093dc2ca1c66c960bfb8c1e8b37270d0 | swap run, exit 0 |
| `wang_ap1_r3.log` | c6e38dc014792f86b98adcc5c380f69f | AP-1 main run, exit 0 |
| `wang_ap1_s39_r3.log` | 0debbf2d8c77dc9ee6c1ceb9ccc2f208 | AP-1 swap run, exit 0 |
| `wang_37_correction_r3.md` | (this memo) | — |

Scripts (in `04-gbmap-wang/`): `wang_analysis_r3.py` (8643a58c…), `wang_analysis_s39_r3.py`
(f2ac515d…), `wang_ap1_r3.py` (9a655038…), `wang_ap1_s39_r3.py` (2e91c592…),
`compute_stats_r3.py` (0703f225…). Backup of pre-correction cached CSVs:
`~/workspace/bioinfo_downloads/wang2022/wang_analysis/backup_pre_r3_20261109/` (4 files,
submitted originals untouched). Submitted originals and `_r2` files in
`revision-package-2026-11-09/` were never modified.

## 9. What was NOT done (per instructions)

- DOCX not edited. No GitHub push. No contact with Cell Reports or bioRxiv. Figures not
  re-rendered. No new pair-table version beyond _r3 (the _r2 table in the revision package
  is untouched).
