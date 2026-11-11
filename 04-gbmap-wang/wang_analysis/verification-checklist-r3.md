# Round-2 verification checklist — 2026-11-09

Independent verification of the round-2 revision (DOCX + Fig 7 + PDF + package).
All checks run against the built artifacts, not the builder's claims.

- **R3-1 — Fig 7 re-rendered from BASE files: PASS.** `figures_final_r3/fig7.png`
  rendered by `plot_figures_r3.py` from `wang_per_sample_s39_r3.csv` +
  `wang_ap1_s39_r3.csv` only. 2100×2340 px, RGB, 300 dpi, md5
  `73530ed8d8101690206a9b31c84e915e`. Visually inspected: panels C
  ("Median diff −0.070 (21 down / 8 up), Wilcoxon p = 0.0075, permutation
  p = 0.072"), C-right (asterisks MM1*/MM4*/MM6*/MM11*; MM8 bare; title
  "* = Bonferroni"), D ("Median diff +0.038 (18/29 up, p = 0.076), FOSL1 not
  significant (p = 0.053), AP-1 vs TRAJ ρ = 0.19 (p = 0.32)"), D-right
  ("VEGFA ρ = 0.413, p = 0.0031, VCAN ρ = 0.489, p = 0.00036") all legible.
  Panels A, B unchanged. Annotation numbers independently recomputed from the
  `_s39_r3` CSVs (TRAJ median diff −0.0699, 8 up / 21 down, Wilcoxon
  p = 0.007456; modules MM1 0.001527 / MM4 0.00234 / MM6 0.0004988 /
  MM11 0.001816 BONF, MM8 0.01058; VEGFA 0.413/0.0031, VCAN 0.489/0.0004,
  n = 49) — all match the memo.
- **R3-2 — DOCX diff is exactly the 6 intended edits: PASS.** Paragraph-level
  diff r2 → r3 = 6 changed paragraphs only: P57 (Results Wang), P58
  (ligand/AP-1/#39), P73 (Limitations), P103 (Fig 7 legend), P113 (Table S7
  legend), P138 (Methods Wang). Summary, tables, references, Key Resources
  Table, ethics statement untouched. New DOCX md5
  `6d1c8f0440efe70feb4f056094b3e8d2` (matches builder's report).
- **R3-3 — Key new numbers present: PASS.** "-0.070" ×3 (P57, P103, P138);
  "21 down, 8 up" / "21 pairs down, 8 up"; "p = 0.0075" ×5; "p = 0.072" ×3;
  "MM1, MM4, MM6, and MM11"; "FOSL1 not significant" + "p = 0.053"; "+0.038";
  "18 of 29 up"; "p = 0.076"; "rho = 0.19 (p = 0.32)"; "VEGFA rho = 0.413";
  "VCAN rho = 0.489"; "GSM5319563 (SF12407), the first recurrence";
  "p = 0.83" (T-cell); "wang_per_sample_s39_r3.csv"; "wang_33pair_table_r3.csv";
  "SF9798 primary, SF9494 recurrent" (Methods #37 disclosure sentence).
- **R3-4 — Forbidden strings absent (Wang context): PASS.** Zero hits for:
  "-0.053", "19 down, 10 up", "19 pairs down, 10 up", "only MM1 and MM8
  remain significant", "FOSL1 lower", "FOSL1 was nominally lower",
  "The analysis uses GSM5319574", "VEGFA rho = 0.408", "VCAN rho = 0.470",
  "MM6 additionally survives", "core_all", "all-cells AP-1", "0.0533",
  "0.0179", "0.1276", "p = 0.035", "p = 0.20". Remaining "0.018"/"0.13"
  substrings are all non-Wang: "-0.018" (Point 4 per-drug), "0.0187"
  (Table S6 xenium fractions), "-0.136" (Fig 2/3 legends), "p = 0.134"
  (v1 LOTO) — each checked in context.
- **R3-5 — TODO/PENDING/[VERIFY] = 0: PASS.**
- **R3-6 — No AP-1 core_all content: PASS.** Zero hits for "core_all" /
  "all-cells AP-1" in the DOCX; the post-hoc p = 0.0308 finding appears
  nowhere in the manuscript or figures.
- **R3-7 — PDF: PASS.** LibreOffice render of the round-2 DOCX: 23 pages,
  Letter, md5 `4019e271f17a713faf95dbb7a1909389`.
- **R3-8 — Phase B files untouched: PASS.** Parent package DOCX md5
  `3811f850…` (Phase B value), PDF md5 `63430a8a…`, mtimes 2026-11-09 01:12–
  01:13; `figures/fig7.png` mtime unchanged. `submission-package/` untouched.
- **R3-9 — _r3 data files in round-2/data match the memo: PASS.** All five
  data-file md5s match `wang_37_correction_r3.md` §8 exactly; fig7.png copy
  matches the render byte-for-byte.
- **R3-10 — Correction-note headers: PASS.** Prepended (additive only) to the
  Phase A memo (`reanalysis-results-2026-11-09.md`, Wang section superseded)
  and `wang_39_sensitivity_r2.md` (roles flipped); both point to
  `wang_37_correction_r3.md`. No other content in either file was modified.

## Open items (carried, not resolved by this round)

- **AP-1 core_all (p = 0.0308, nominally significant, both branches):** post-hoc,
  not pre-specified; deliberately left out of the DOCX and figures. Parent's
  call: (a) leave unreported or (b) report as exploratory with caveat.
- **C-9 (Phase B):** the STAR sentence "all other Wilcoxon p-values were
  re-derived under this convention" was still not fully re-derived this round;
  only the changed Wang values were verified. The sentence is unchanged.
- **O-1b / O-6 (Phase B):** Wang-paper recurrence order and Points 1–4 design
  docs still not re-verified; neither affects the changed numbers.
- **GitHub deposit:** still blocked on the user's "现在推" + fresh 30-min
  approval. Note: the _r3 files are a separate, later deposit decision; the
  _r2/_s39 deposit the user was asked about is now superseded in content by
  _r3 (the parent should re-confirm scope before any push).
- **Journal/bioRxiv route:** still pending the user's decision.
