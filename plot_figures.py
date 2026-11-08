#!/usr/bin/env python3
"""Publication-standard redraw of Figures 1-7 for:
    "Drug-response geometry along the glioblastoma trajectory" (Kong, Cell Reports submission)

Style conventions adapted from ChenLiu-1996/figures4papers
(https://github.com/ChenLiu-1996/figures4papers), scientific-figure-making skill:
minimalist spines (top/right off), frameless legends, semantic palette
(blue = key finding, red = contrast/null, neutral = background),
in-place annotation of key statistics, tight_layout, dpi=300 export.
Adapted (not copied): font sizes rescaled for 7-inch journal figures
(the repo's 15-24 pt presets assume 14-45 inch canvases), no use_tex
(LaTeX not installed locally), no ultra-wide aspect ratios,
7-inch full-width layout throughout.

All plotted values come from the result CSVs on disk or from the manuscript
legend values (the latter are marked LEGEND in comments). No analysis is
re-run here except trivial aggregations (means/medians/Spearman) recomputed
from the CSVs to place points; every recomputed statistic was checked against
the legend before plotting (see FIGURES_README.md).

Run with the project venv (has scanpy/anndata for the Fig 1B histogram):
    ~/workspace/bioinfo_downloads/gbm-point1/venv/bin/python plot_figures.py
Outputs: figures_final/fig1.png ... fig7.png (300 dpi, RGB).
"""

import os
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import pandas as pd
from scipy import stats

warnings.filterwarnings("ignore", category=UserWarning)

BASE = os.path.dirname(os.path.abspath(__file__))
# The script lives inside figures_final/; outputs go next to it.
OUT = BASE if os.path.basename(BASE) == "figures_final" else os.path.join(BASE, "figures_final")
os.makedirs(OUT, exist_ok=True)

P1 = os.path.join(os.path.expanduser("~"), "workspace/bioinfo_downloads/gbm-point1")
MEY = os.path.join(os.path.expanduser("~"), "workspace/bioinfo_downloads/gse296577")
GSE330 = os.path.join(os.path.expanduser("~"), "workspace/bioinfo_downloads/gse330794")
GBMAP = os.path.join(os.path.expanduser("~"), "workspace/bioinfo_downloads/gbmap")
WANG = os.path.join(os.path.expanduser("~"), "workspace/bioinfo_downloads/wang2022")

# ---------------------------------------------------------------- style ---
# figures4papers house style, rescaled for 7-inch figures.
PALETTE = {
    "blue_main": "#0F4D92",
    "blue_secondary": "#3775BA",
    "blue_light": "#9BC8FA",
    "red_strong": "#B64342",
    "red_light": "#E9A6A1",
    "green": "#8BCF8B",
    "teal": "#42949E",
    "violet": "#9A4D8E",
    "neutral": "#767676",
    "dark": "#4D4D4D",
    "light_gray": "#D9D9D9",
}
mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Liberation Sans", "DejaVu Sans"],
    "font.size": 9,
    "axes.titlesize": 9,
    "axes.labelsize": 9,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "legend.fontsize": 8,
    "axes.linewidth": 1.0,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "legend.frameon": False,
    "savefig.dpi": 300,
    "figure.dpi": 300,
    "svg.fonttype": "none",
})

BLUE = PALETTE["blue_main"]
BLUE2 = PALETTE["blue_secondary"]
RED = PALETTE["red_strong"]
GRAY = PALETTE["neutral"]
DARK = PALETTE["dark"]
TEAL = PALETTE["teal"]
VIOLET = PALETTE["violet"]


def panel_letter(ax, letter):
    """Bold panel letter inside the top-left corner of a panel."""
    ax.text(0.02, 0.96, letter, transform=ax.transAxes, fontsize=11,
            fontweight="bold", va="top", ha="left",
            bbox=dict(facecolor="white", alpha=0.85, edgecolor="none", pad=1.5))


def panel_title(ax, text):
    ax.set_title(text, fontsize=9, pad=4, loc="left")


def annotate_stats(ax, text, x=0.95, y=0.97, ha="right", va="top",
                   fontsize=8, bbox=True):
    """Key statistics on the panel. Default: top-right, white backing box
    (keeps text readable over sparse data; the panel letter is top-left)."""
    kw = dict(transform=ax.transAxes, fontsize=fontsize,
              ha=ha, va=va, color=DARK)
    if bbox:
        kw["bbox"] = dict(facecolor="white", alpha=0.85, edgecolor="none",
                          pad=1.5)
    ax.text(x, y, text, **kw)


def save(fig, name):
    """Save at 300 dpi and force RGB (no alpha channel)."""
    path = os.path.join(OUT, name)
    fig.savefig(path, dpi=300, facecolor="white")
    plt.close(fig)
    # Guarantee RGB mode (matplotlib can write RGBA when alpha is present)
    # and embed 300 dpi metadata: matplotlib's PNG writer omits the pHYs
    # chunk, so re-save through PIL with dpi=(300, 300).
    from PIL import Image
    im = Image.open(path)
    if im.mode != "RGB":
        im = im.convert("RGB")
    im.save(path, dpi=(300, 300))
    im2 = Image.open(path)
    print(f"saved {name}: {im2.size[0]}x{im2.size[1]} px, mode {im2.mode}, "
          f"dpi={im2.info.get('dpi')}")


def style_ax(ax):
    ax.tick_params(length=3, width=0.8)
    for s in ("left", "bottom"):
        ax.spines[s].set_linewidth(1.0)
        ax.spines[s].set_color("black")


# ---------------------------------------------------------------- data ----
def load_csv(*parts):
    return pd.read_csv(os.path.join(*parts))


# ============================================================== FIGURE 1 ==
def fig1():
    # --- data ---
    import anndata as ad
    h5 = ad.read_h5ad(os.path.join(P1, "b-side/adata_tumor_selected_full.h5ad"),
                      backed="r")
    traj_cells = np.asarray(h5.obs["TRAJ"], dtype=float)
    qb = load_csv(P1, "colab_outputs/quintile_bins.csv")
    qs = load_csv(P1, "colab_outputs/quintile_summary.csv")
    gb = load_csv(P1, "colab_outputs/geom_b_reproduction.csv")
    # quintile edges from the transition-matrix memo (TRAJ quintiles of all cells)
    edges = np.quantile(traj_cells, [0.2, 0.4, 0.6, 0.8])

    fig = plt.figure(figsize=(7.0, 6.6), constrained_layout=True)
    gs = fig.add_gridspec(3, 2, height_ratios=[1, 1, 1], hspace=0.3, wspace=0.35)

    # (A) TRAJ definition schematic
    ax = fig.add_subplot(gs[0, 0])
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.axis("off")
    # Full-width linear layout (robust to font metrics: the module lists use
    # the whole panel width instead of being squeezed into half-width boxes)
    ax.annotate("", xy=(9.7, 8.2), xytext=(0.3, 8.2),
                arrowprops=dict(arrowstyle="-|>", lw=1.6, color=DARK))
    ax.text(1.5, 8.9, "early", ha="left", fontsize=8.5, color=BLUE)
    ax.text(9.7, 8.9, "late", ha="right", fontsize=8.5, color=RED)
    ax.text(0.3, 7.0, "EARLY modules", ha="left", fontsize=9, color=BLUE,
            fontweight="bold")
    ax.text(0.3, 6.1, "MM1 \u00b7 MM4 \u00b7 MM6 \u00b7 MM7 \u00b7 MM8 \u00b7 MM10 \u00b7 MM11",
            ha="left", fontsize=7.5)
    ax.text(0.3, 4.9, "LATE modules", ha="left", fontsize=9, color=RED,
            fontweight="bold")
    ax.text(0.3, 4.0, "MM2 \u00b7 MM3 \u00b7 MM5", ha="left", fontsize=7.5)
    ax.text(5.0, 2.6, "TRAJ = mean(LATE) \u2212 mean(EARLY)", ha="center",
            fontsize=9)
    ax.text(5.0, 1.7, "scanpy score_genes \u00b7 Table S10 gene lists",
            ha="center", fontsize=8, color=GRAY)
    ax.text(5.0, 0.8, "MM0, MM9 excluded", ha="center", fontsize=8, color=GRAY)
    panel_letter(ax, "A")
    panel_title(ax, "TRAJ definition")

    # (B) TRAJ distribution across 252,106 tumor cells
    ax = fig.add_subplot(gs[0, 1])
    ax.hist(traj_cells, bins=60, color=BLUE2, edgecolor="white", linewidth=0.3)
    med = float(np.median(traj_cells))
    ax.axvline(med, color=RED, lw=1.2, ls="--")
    for e in edges:
        ax.axvline(e, color=GRAY, lw=0.7, ls=":")
    ax.set_xlabel("TRAJ")
    ax.set_ylabel("Tumor cells")
    annotate_stats(ax, f"n = 252,106 cells\nmedian = {med:.3f}")
    panel_letter(ax, "B")
    panel_title(ax, "TRAJ distribution")
    style_ax(ax)

    # (C) direction vs trajectory quintile (all 214 quintile observations)
    ax = fig.add_subplot(gs[1, 0])
    rng = np.random.default_rng(7)
    xj = qb["quintile"].values + rng.uniform(-0.14, 0.14, size=len(qb))
    ax.scatter(xj, qb["disp_traj"].values, s=9, color=BLUE2, alpha=0.45,
               linewidths=0, zorder=2)
    ax.errorbar(qs["quintile"], qs["mean_disp"], yerr=qs["sem_disp"],
                fmt="o-", color=BLUE, ms=5, lw=1.6, capsize=3, zorder=3,
                label="Mean \u00b1 SEM")
    ax.axhline(0, color="black", lw=0.8)
    ax.set_xlabel("Trajectory quintile (1 = earliest)")
    ax.set_ylabel("Displacement (d$_{traj}$)")
    ax.set_xticks([1, 2, 3, 4, 5])
    annotate_stats(ax, "Pooled Spearman\n$\\rho$ = \u22120.242, p = 3.6$\\times$10$^{-4}$\nn = 214 quintile obs.\nBlue: mean \u00b1 SEM")
    panel_letter(ax, "C")
    panel_title(ax, "Direction vs position")
    style_ax(ax)

    # (D) five-bin means
    ax = fig.add_subplot(gs[1, 1])
    ax.errorbar(qs["quintile"], qs["mean_disp"], yerr=qs["sem_disp"],
                fmt="o-", color=BLUE, ms=6, lw=1.8, capsize=4)
    ax.axhline(0, color="black", lw=0.8)
    ax.set_xlabel("Trajectory quintile")
    ax.set_ylabel("Mean displacement")
    ax.set_xticks([1, 2, 3, 4, 5])
    annotate_stats(ax, "Five-bin means\n$\\rho$ = \u22120.894, p = 0.041\nmedian d$_{traj}$ = \u22120.019\n14/43 toward late")
    panel_letter(ax, "D")
    panel_title(ax, "Five-bin means")
    style_ax(ax)

    # (E) validation: direction C reproduced + TRAJ robustness
    gs_e = gs[2, 0:2].subgridspec(1, 2, wspace=0.35)
    ax = fig.add_subplot(gs_e[0, 0])
    rng2 = np.random.default_rng(3)
    xj = rng2.uniform(-0.18, 0.18, size=len(gb))
    ax.scatter(xj, gb["C"].values, s=10, color=BLUE2, alpha=0.55,
               linewidths=0, zorder=2)
    ax.axhline(float(gb["C"].median()), color=RED, lw=1.4, ls="--")
    ax.set_xlim(-0.6, 0.6)
    ax.set_xticks([])
    ax.set_ylabel("Direction C")
    ax.set_ylim(0, 1.0)
    annotate_stats(ax, "C reproduced exactly\nmedian = 0.608\n43/43 positive", fontsize=7.5)
    panel_letter(ax, "E")
    panel_title(ax, "Direction C")
    style_ax(ax)

    ax = fig.add_subplot(gs_e[0, 1])
    labels = ["TRAJ vs\nMM9 variant", "TRAJ vs\nmodule PC1"]
    vals = [0.998, 0.937]  # LEGEND values (robustness checks from the analysis)
    bars = ax.bar([0, 1], vals, width=0.55, color=[BLUE2, TEAL],
                  edgecolor="black", linewidth=0.8)
    for b, v in zip(bars, vals):
        ax.text(b.get_x() + b.get_width() / 2, v + 0.03, f"{v:.3f}",
                ha="center", fontsize=8)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(labels, fontsize=7.5)
    ax.set_ylabel("Spearman $\\rho$")
    ax.set_ylim(0, 1.12)
    panel_title(ax, "TRAJ robustness")
    style_ax(ax)

    save(fig, "fig1.png")


# ============================================================== FIGURE 2 ==
def fig2():
    qb = load_csv(P1, "colab_outputs/quintile_bins.csv")
    qs = load_csv(P1, "colab_outputs/quintile_summary.csv")
    disp = load_csv(P1, "colab_outputs/dispersion_state_vs_traj.csv")
    dtr = load_csv(P1, "colab_outputs/decomposition_traj_groups.csv")
    dst = load_csv(P1, "colab_outputs/decomposition_state_groups.csv")
    tm = load_csv(P1, "colab_outputs/transition_matrix.csv")
    tpd = load_csv(P1, "colab_outputs/transition_per_drug.csv")

    fig = plt.figure(figsize=(7.0, 6.6), constrained_layout=True)
    gs = fig.add_gridspec(3, 2, height_ratios=[1, 1, 1], hspace=0.3, wspace=0.35)

    # (A) amplitude vs trajectory quintile
    ax = fig.add_subplot(gs[0, 0])
    rng = np.random.default_rng(11)
    xj = qb["quintile"].values + rng.uniform(-0.14, 0.14, size=len(qb))
    ax.scatter(xj, qb["mag12"].values, s=9, color=GRAY, alpha=0.5,
               linewidths=0, zorder=2)
    ax.errorbar(qs["quintile"], qs["mean_mag"], yerr=qs["sem_mag"],
                fmt="o-", color=DARK, ms=5, lw=1.6, capsize=3, zorder=3,
                label="Mean \u00b1 SEM")
    ax.set_xlabel("Trajectory quintile (1 = earliest)")
    ax.set_ylabel("Amplitude (m12)")
    ax.set_xticks([1, 2, 3, 4, 5])
    annotate_stats(ax, "Pooled Spearman\n$\\rho$ = 0.064, p = 0.35\nn = 214 quintile obs.\nDark: mean \u00b1 SEM")
    panel_letter(ax, "A")
    panel_title(ax, "Amplitude vs position")
    style_ax(ax)

    # (B) dispersion: state vs trajectory coordinates
    ax = fig.add_subplot(gs[0, 1])
    lims = [0, 0.075]
    ax.scatter(disp["sd_state"], disp["sd_traj"], s=14, color=BLUE2,
               alpha=0.65, linewidths=0, zorder=2)
    # identity line stops just inside the corner so its cap never escapes
    ax.plot([0, 0.070], [0, 0.070], ls="--", color=RED, lw=1.1)
    ax.set_xlim(lims)
    ax.set_ylim(lims)
    ax.set_xlabel("Dispersion, state coords (SD)")
    ax.set_ylabel("Traj. dispersion (SD)")
    annotate_stats(ax, "Median 0.012 \u2192 0.008\n56% improved\nWilcoxon p = 0.031", fontsize=7.5)
    panel_letter(ax, "B")
    panel_title(ax, "Amplitude dispersion")
    style_ax(ax)

    # (C) decomposition: traj-bin vs state-group (medians)
    ax = fig.add_subplot(gs[1, 0])
    med_tr_w = float(dtr["within_traj"].median())
    med_tr_c = float(dtr["composition_traj"].median())
    med_st_w = float(dst["within_traj"].median())
    med_st_c = float(dst["composition_traj"].median())
    share_tr = float(dtr["frac_composition"].median())
    share_st = float(dst["frac_composition"].median())
    x = np.arange(2)
    w = 0.36
    ax.bar(x - w / 2, [med_tr_w, med_st_w], w, color=BLUE2,
           edgecolor="black", linewidth=0.8, label="Within")
    ax.bar(x + w / 2, [med_tr_c, med_st_c], w, color="#E8913A",
           edgecolor="black", linewidth=0.8, label="Composition")
    ax.axhline(0, color="black", lw=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels([f"Trajectory bins\n{share_tr:.0%} composition",
                        f"State groups\n{share_st:.0%} composition"], fontsize=8)
    ax.set_ylabel("Median displacement")
    for xi, v in zip(x - w / 2, [med_tr_w, med_st_w]):
        ax.text(xi, v / 2, f"{v:.3f}", ha="center", va="center",
                fontsize=7.5, color="white", fontweight="bold")
    for xi, v in zip(x + w / 2, [med_tr_c, med_st_c]):
        ax.text(xi, v / 2, f"{v:.3f}", ha="center", va="center",
                fontsize=7.5, color="white", fontweight="bold")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.30),
              handlelength=1.0, ncol=2)
    panel_letter(ax, "C")
    panel_title(ax, "Displacement decomposition")
    style_ax(ax)

    # (D) transition matrix: mean bin-share change (aggregate)
    ax = fig.add_subplot(gs[1, 1])
    agg = tm[tm["bin"] < 5].groupby("bin")["delta"].mean()
    bins = [0, 1, 2, 3, 4]
    vals = [agg[b] for b in bins]
    early = float(tm[tm["bin"] == 99]["delta"].mean())
    late = float(tm[tm["bin"] == 98]["delta"].mean())
    labels = ["Bin 0 (ns)", "Bin 1 (p = 0.0003)", "Bin 2 (ns)",
              "Bin 3 (p = 0.033)", "Bin 4 (p = 0.0011)",
              "Early 0\u20131 (p = 0.015)", "Late 3\u20134 (p = 0.0003)"]
    all_vals = vals + [early, late]
    ypos = list(range(5)) + [6, 7]
    colors = [BLUE if v > 0 else RED for v in all_vals]
    bars = ax.barh(ypos, all_vals, height=0.6, color=colors, edgecolor="black",
                   linewidth=0.8)
    # hatched edges mark the two aggregate (summary) bars
    bars[5].set_hatch("///")
    bars[6].set_hatch("///")
    for b, v in zip(bars, all_vals):
        if v > 0:
            ax.text(v + 0.004, b.get_y() + b.get_height() / 2, f"{v:+.3f}",
                    ha="left", va="center", fontsize=8)
        else:
            ax.text(v / 2, b.get_y() + b.get_height() / 2, f"{v:+.3f}",
                    ha="center", va="center", fontsize=8, color="white",
                    fontweight="bold")
    ax.axvline(0, color="black", lw=0.8)
    ax.set_yticks(ypos)
    ax.set_yticklabels(labels, fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel("Mean bin-share change")
    ax.set_xlim(-0.085, 0.085)
    panel_letter(ax, "D")
    panel_title(ax, "Transition matrix (aggregate)")
    style_ax(ax)

    # (E) per-drug transition profiles
    gs_e = gs[2, 0:2].subgridspec(1, 1)
    ax = fig.add_subplot(gs_e[0, 0])
    highlight = {"topotecan": (BLUE, 2.0), "ANA-12": (RED, 2.0)}
    for _, row in tpd.iterrows():
        drug = row["drug"]
        yv = [row[str(b)] for b in bins]
        if drug in highlight:
            c, lw = highlight[drug]
            ax.plot(bins, yv, "o-", color=c, lw=lw, ms=4, zorder=3)
        else:
            ax.plot(bins, yv, "o-", color="#BBBBBB", lw=0.9, ms=2.5, zorder=2)
    ax.axhline(0, color="black", lw=0.8)
    ax.set_xlabel("Trajectory bin (0 = earliest, 4 = latest)")
    ax.set_ylabel("Bin-share change")
    ax.set_xticks(bins)
    # direct labels for the two highlighted drugs (placed at their bin-1 peaks)
    ax.text(1.08, 0.128, "topotecan", color=BLUE, fontsize=8, va="bottom")
    ax.text(4.0, -0.150, "ANA-12", color=RED, fontsize=8, va="top", ha="right")
    annotate_stats(ax, "Topotecan: bin 1 +0.125, bin 3 \u22120.136\nANA-12: bin 4 \u22120.165\nOther drugs weak (gray)",
                   x=0.02, y=0.03, ha="left", va="bottom", fontsize=7.5)
    panel_letter(ax, "E")
    panel_title(ax, "Per-drug transition profiles")
    style_ax(ax)

    save(fig, "fig2.png")


if __name__ == "__main__":
    fig1()
    fig2()
    print("fig1, fig2 done")

# ============================================================== FIGURE 3 ==
def fig3():
    lopo_dir = load_csv(P1, "colab_outputs/lopo_direction.csv")
    lodo_sh = load_csv(P1, "colab_outputs/lodo_shape.csv")
    v2_lopo = load_csv(P1, "colab_outputs2v2/v2_lopo.csv")
    v2_lodo = load_csv(P1, "colab_outputs2v2/v2_lodo.csv")
    v3_lopo = load_csv(P1, "colab_outputs2v3/v3_lopo.csv")
    v3_unit = load_csv(P1, "colab_outputs2v3/v3_lopo_unit.csv")
    v3_strata = load_csv(P1, "colab_outputs2v3/v3_strata_summary.csv")
    p4 = load_csv(P1, "colab_outputs2v2/point4_frame.csv")

    fig = plt.figure(figsize=(7.0, 7.8), constrained_layout=True)
    gs = fig.add_gridspec(4, 2, height_ratios=[1, 1, 1, 1], hspace=0.25, wspace=0.35)

    # (A) v1: LOPO direction + LODO shape
    ax = fig.add_subplot(gs[0, 0])
    # pred/obs are stored as stringified 5-bin vectors per comparison;
    # flatten to the 215 per-bin observed-vs-predicted points.
    import ast
    obs_all, pred_all = [], []
    for _, row in lopo_dir.iterrows():
        obs_all.extend(ast.literal_eval(row["obs"]))
        pred_all.extend(ast.literal_eval(row["pred"]))
    obs_all = np.array(obs_all)
    pred_all = np.array(pred_all)
    print("Fig3A pooled Pearson (flattened obs vs pred):",
          round(float(stats.pearsonr(obs_all, pred_all).statistic), 3))
    lims = [min(obs_all.min(), pred_all.min()) * 1.1,
            max(obs_all.max(), pred_all.max()) * 1.1]
    ax.scatter(obs_all, pred_all, s=12, color=BLUE2,
               alpha=0.7, linewidths=0, zorder=2)
    ax.plot([lims[0], lims[1] * 0.95], [lims[0], lims[1] * 0.95], ls="--", color=RED, lw=1.0)
    ax.set_xlim(lims)
    ax.set_ylim(lims)
    ax.set_xlabel("Observed direction")
    ax.set_ylabel("Predicted direction")
    annotate_stats(ax, "LOPO direction\nmedian r = 0.124\npooled r = 0.003\n44% beat flat, p = 0.134",
                   x=0.98, y=0.03, ha="right", va="bottom", fontsize=7.5)
    panel_letter(ax, "A")
    panel_title(ax, "v1 LOPO direction")
    style_ax(ax)

    ax = fig.add_subplot(gs[0, 1])
    drugs = lodo_sh["drug"].tolist()
    vals = lodo_sh["pearson"].tolist()
    cols = [RED if v < 0 else BLUE2 for v in vals]
    ax.barh(np.arange(len(drugs)), vals, height=0.62, color=cols,
            edgecolor="black", linewidth=0.7)
    ax.axvline(0.004, color=RED, lw=1.0, ls="--")
    ax.set_yticks(np.arange(len(drugs)))
    ax.set_yticklabels(drugs, fontsize=7)
    ax.invert_yaxis()
    ax.set_xlabel("LODO shape Pearson r")
    ax.set_title("v1 LODO shape \u2014 median r = 0.699\nnull 0.004, p = 0.0\nRMSE 0.0024 vs 0.0027 (p = 0.652)\nMagnitude 0.0234 vs 0.0222 (p = 0.578)",
                 fontsize=7.5, pad=4, loc="left")
    style_ax(ax)

    # (B) v2: LOPO amplitude + LODO feature combos
    ax = fig.add_subplot(gs[1, 0])
    lims = [0, max(v2_lopo["y"].max(), v2_lopo["pred_model"].max()) * 1.1]
    ax.scatter(v2_lopo["y"], v2_lopo["pred_null"], s=12, color=GRAY,
               alpha=0.6, linewidths=0, label="Drug-mean null", zorder=2)
    ax.scatter(v2_lopo["y"], v2_lopo["pred_model"], s=12, color="#E8913A",
               alpha=0.8, linewidths=0, label="Ridge (baseline feats)", zorder=3)
    ax.plot([lims[0], lims[1] * 0.95], [lims[0], lims[1] * 0.95], ls="--", color=RED, lw=1.0)
    ax.set_xlim(lims)
    ax.set_ylim(lims)
    ax.set_xlabel("Observed amplitude (m12)")
    ax.set_ylabel("Predicted amplitude")
    annotate_stats(ax, "LOPO amplitude (n = 41)\nRidge RMSE 0.0687, r = 0.065\nNull RMSE 0.0615, r = 0.129\nWilcoxon p = 0.96\nGray: drug-mean null; orange: ridge",
                   x=0.98, y=0.03, ha="right", va="bottom", fontsize=7.5)
    panel_letter(ax, "B")
    panel_title(ax, "v2 LOPO amplitude")
    style_ax(ax)

    ax = fig.add_subplot(gs[1, 1])
    feats = v2_lodo["features"].tolist()
    vals = v2_lodo["rmse"].tolist()
    cols = [BLUE if f == "lig+traj" else BLUE2 for f in feats]
    ax.barh(np.arange(len(feats)), vals, height=0.62, color=cols,
            edgecolor="black", linewidth=0.7, alpha=0.85)
    ax.axvline(0.0291, color=RED, lw=1.1, ls="--")
    ax.set_yticks(np.arange(len(feats)))
    ax.set_yticklabels(feats, fontsize=7)
    ax.invert_yaxis()
    ax.set_xlabel("LODO RMSE")
    ax.set_title("v2 LODO feature combos\nbest lig+traj 0.0288\nvs global 0.0291 (= noise)",
                 fontsize=7.5, pad=4, loc="left")
    style_ax(ax)

    # (C) v3: strata strip + RMSE comparison (comparison FAILED / unit PASSED)
    ax = fig.add_subplot(gs[2, 0])
    strata = ["low", "mid", "high"]
    means = [float(v3_strata[v3_strata["stratum"] == s]["mean_m12"].iloc[0])
             for s in strata]
    cols = [BLUE2, "#E8913A", RED]
    rng = np.random.default_rng(5)
    for i, s in enumerate(strata):
        vals = v3_lopo[v3_lopo["stratum"] == s]["y"].values
        xj = i + rng.uniform(-0.16, 0.16, size=len(vals))
        ax.scatter(xj, vals, s=10, color=cols[i], alpha=0.6, linewidths=0)
        ax.plot([i - 0.28, i + 0.28], [means[i], means[i]], color="black", lw=1.8)
    ax.set_xticks([0, 1, 2])
    ax.set_xticklabels(["Low", "Mid", "High"], fontsize=8)
    ax.set_xlabel("Late-score tertile")
    ax.set_ylabel("Amplitude (m12)")
    annotate_stats(ax, "Expected m12\n0.038 / 0.071 / 0.093\nK-W p = 0.0066 (comp.)\np = 0.0285 (per-tissue)",
                   x=0.02, y=0.03, ha="left", va="bottom", fontsize=7.5)
    panel_letter(ax, "C")
    panel_title(ax, "v3 strata")
    style_ax(ax)

    ax = fig.add_subplot(gs[2, 1])
    groups = ["Comparison\nlevel", "Tissue\nlevel"]
    strat = [0.0559, 0.0348]  # LEGEND values (verified from the CSVs)
    glob = [0.0584, 0.0422]
    x = np.arange(2)
    w = 0.34
    ax.bar(x - w / 2, strat, w, color=BLUE, edgecolor="black", linewidth=0.8,
           label="Stratum prediction")
    ax.bar(x + w / 2, glob, w, color=GRAY, edgecolor="black", linewidth=0.8,
           label="Global mean")
    for xi, v in zip(x - w / 2, strat):
        ax.text(xi, v + 0.001, f"{v:.4f}", ha="center", fontsize=7.5)
    for xi, v in zip(x + w / 2, glob):
        ax.text(xi, v + 0.001, f"{v:.4f}", ha="center", fontsize=7.5)
    ax.set_xticks(x)
    ax.set_xticklabels(groups, fontsize=7.5)
    ax.set_ylabel("LOPO RMSE")
    ax.set_ylim(0, 0.085)
    ax.text(0, 0.078, "FAILED (p = 0.23)", ha="center", fontsize=7.5,
            color=RED, fontweight="bold")
    ax.text(1, 0.078, "PASSED (p = 0.030)", ha="center", fontsize=7.5,
            color=BLUE, fontweight="bold")
    ax.set_title("v3 RMSE\nblue: stratum; gray: global",
                 fontsize=8, pad=4, loc="left")
    style_ax(ax)

    # (D) late score vs amplitude
    ax = fig.add_subplot(gs[3, 0:2])
    class_colors = {"GPX4i": TEAL, "HDACi": VIOLET, "other": GRAY}
    for cls, g in p4.groupby("drug_class"):
        ax.scatter(g["late_score"], g["m12"], s=14, color=class_colors.get(cls, GRAY),
                   alpha=0.75, linewidths=0, label=cls, zorder=3)
    # pooled regression line
    m, b = np.polyfit(p4["late_score"], p4["m12"], 1)
    xs = np.array([p4["late_score"].min(), p4["late_score"].max()])
    ax.plot(xs, m * xs + b, color=RED, lw=1.2, ls="--")
    ax.set_xlabel("Late-module score (vehicle arm)")
    ax.set_ylabel("Amplitude (m12)")
    annotate_stats(ax, "Pooled $\\rho$ = 0.360, p = 0.021 (n = 41)\nPer-TissueID $\\rho$ = 0.640, p = 0.006 (n = 17)\nMM5 $\\rho$ = 0.440, MM2 0.413, MM3 0.234",
                   x=0.98, y=0.03, ha="right", va="bottom", fontsize=7.5)
    ax.legend(loc="upper left", bbox_to_anchor=(0.0, 0.88), handlelength=1.0,
              fontsize=7, title="Drug class", title_fontsize=7)
    panel_letter(ax, "D")
    panel_title(ax, "Late score vs amplitude")
    style_ax(ax)

    save(fig, "fig3.png")


# ============================================================== FIGURE 4 ==
def fig4():
    mf = load_csv(P1, "colab_outputs/myeloid_fraction.csv")
    cd = load_csv(P1, "colab_outputs/myeloid_centroid_distance.csv")
    mvr = load_csv(P1, "colab_outputs/myeloid_vs_tumor_response.csv")

    fig = plt.figure(figsize=(7.0, 6.2), constrained_layout=True)
    gs = fig.add_gridspec(3, 2, height_ratios=[1, 1, 1], hspace=0.3, wspace=0.35)

    # (A) myeloid fraction: treated vs vehicle
    ax = fig.add_subplot(gs[0, 0])
    lims = [0, max(mf["treated"].max(), mf["vehicle"].max()) * 1.08]
    ax.scatter(mf["vehicle"], mf["treated"], s=16, color=BLUE2, alpha=0.75,
               linewidths=0, zorder=2)
    ax.plot([0, lims[1] * 0.95], [0, lims[1] * 0.95], ls="--", color=RED, lw=1.0)
    ax.set_xlim(lims)
    ax.set_ylim(lims)
    ax.set_xlabel("Myeloid fraction, vehicle")
    ax.set_ylabel("Myeloid fraction, treated")
    # TB5980 is the detection-confounded outlier (vehicle 0.875, treated 0.091)
    tb = mf[mf["tis"] == "TB5980"].iloc[0]
    ax.scatter([tb["vehicle"]], [tb["treated"]], s=60, facecolor="none",
               edgecolor=RED, linewidth=1.4, zorder=3)
    ax.annotate("TB5980", (tb["vehicle"], tb["treated"]),
                xytext=(tb["vehicle"] - 0.28, tb["treated"] + 0.12),
                fontsize=7.5, color=RED)
    annotate_stats(ax, "Median $\\Delta$ = 0.000\nWilcoxon p = 0.88\nControl (oligodendrocyte) p = 0.23", fontsize=7.5)
    panel_letter(ax, "A")
    panel_title(ax, "Myeloid composition unchanged")
    style_ax(ax)

    # (B) shift magnitude vs null
    ax = fig.add_subplot(gs[0, 1])
    lims = [0, max(cd["dist"].max(), cd["null_median"].max()) * 1.08]
    ax.scatter(cd["null_median"], cd["dist"], s=16, color=BLUE2, alpha=0.75,
               linewidths=0, zorder=2)
    ax.plot([0, lims[1] * 0.95], [0, lims[1] * 0.95], ls="--", color=RED, lw=1.0)
    ax.set_xlim(lims)
    ax.set_ylim(lims)
    ax.set_xlabel("Null median shift")
    ax.set_ylabel("Observed centroid shift")
    annotate_stats(ax, "15/15 above null\nMedian ratio 3.27\u00d7\nWilcoxon p = 6.1$\\times$10$^{-5}$", fontsize=7.5)
    panel_letter(ax, "B")
    panel_title(ax, "Myeloid state shift is real")
    style_ax(ax)

    # (C) per-drug centroid ratios — LEGEND values (per-TissueID CSV has no drug column)
    ax = fig.add_subplot(gs[1, 0:2])
    drugs = ["topotecan", "panobinostat", "rsl3_ferrostatin", "givinostat",
             "ANA-12", "RO4929097"]
    ratios = [2.232, 1.998, 1.680, 1.590, 1.094, 1.040]
    cols = [BLUE if r > 1.5 else (BLUE2 if r > 1.2 else GRAY) for r in ratios]
    ax.barh(np.arange(len(drugs)), ratios, height=0.62, color=cols,
            edgecolor="black", linewidth=0.7)
    ax.axvline(1.0, color=RED, lw=1.0, ls="--")
    ax.set_yticks(np.arange(len(drugs)))
    ax.set_yticklabels(drugs, fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel("Centroid shift / null (ratio)")
    for i, v in enumerate(ratios):
        ax.text(v + 0.04, i, f"{v:.3f}", va="center", fontsize=8)
    ax.set_xlim(0, 2.6)
    annotate_stats(ax, "Direction inconsistent:\nmean pairwise Spearman 0.105",
                   x=0.98, y=0.03, ha="right", va="bottom", fontsize=7.5)
    panel_letter(ax, "C")
    panel_title(ax, "Per-drug centroid ratios")
    style_ax(ax)

    # (D) ligand vs TRAJECTORY (+ inset: shift vs amplitude)
    ax = fig.add_subplot(gs[2, 0])
    ax.scatter(mvr["mean_traj_ctrl"], mvr["lig_vehicle"], s=16, color=BLUE2,
               alpha=0.75, linewidths=0, zorder=2)
    m, b = np.polyfit(mvr["mean_traj_ctrl"], mvr["lig_vehicle"], 1)
    xs = np.array([mvr["mean_traj_ctrl"].min(), mvr["mean_traj_ctrl"].max()])
    ax.plot(xs, m * xs + b, color=RED, lw=1.2, ls="--")
    ax.set_xlabel("Tumor trajectory position")
    ax.set_ylabel("Baseline ligand")
    annotate_stats(ax, "Pooled $\\rho$ = 0.918\np = 3.07$\\times$10$^{-17}$ (n = 41)\nPer-TissueID $\\rho$ = 0.87\np = 5.63$\\times$10$^{-6}$ (n = 17)",
                   x=0.98, y=0.03, ha="right", va="bottom", fontsize=7.5)
    panel_letter(ax, "D")
    panel_title(ax, "Ligand tracks trajectory")
    style_ax(ax)

    ax = fig.add_subplot(gs[2, 1])
    ax.scatter(mvr["m12"], mvr["my_shift"], s=16, color=TEAL, alpha=0.75,
               linewidths=0, zorder=2)
    m, b = np.polyfit(mvr["m12"], mvr["my_shift"], 1)
    xs = np.array([mvr["m12"].min(), mvr["m12"].max()])
    ax.plot(xs, m * xs + b, color=RED, lw=1.2, ls="--")
    ax.set_xlabel("Tumor amplitude (m12)")
    ax.set_ylabel("Myeloid shift magnitude")
    annotate_stats(ax, "$\\rho$ = 0.411\np = 0.0075",
                   x=0.98, y=0.03, ha="right", va="bottom", fontsize=7.5)
    panel_title(ax, "Shift vs amplitude")
    style_ax(ax)

    save(fig, "fig4.png")

# ============================================================== FIGURE 5 ==
def fig5():
    q1c = load_csv(MEY, "q1_composition.csv")
    q1d = load_csv(MEY, "q1_distances.csv")
    q2a = load_csv(MEY, "q2a_ligand_by_tier2.csv")
    t2_post = load_csv(MEY, "meylan_analysis/mm2_tcell_per_section_tier2_post.csv")
    mes2_post = load_csv(MEY, "meylan_analysis/mm2_tcell_per_section_mes2_post.csv")
    proxy_post = load_csv(MEY, "meylan_analysis/mm2_tcell_per_section_post.csv")

    # strip _pre/_post suffix from patient labels before pivoting
    q1c = q1c.copy()
    q1c["patient_base"] = q1c["patient"].str.replace(r"_(pre|post)$", "", regex=True)

    fig = plt.figure(figsize=(7.0, 7.8), constrained_layout=True)
    gs = fig.add_gridspec(4, 2, height_ratios=[1, 1, 1, 1], hspace=0.25, wspace=0.35)

    # (A) T-cell influx: fraction + tumor->T distance
    ax = fig.add_subplot(gs[0, 0])
    tf = q1c.pivot(index="patient_base", columns="timepoint", values="T_frac")
    for _, row in tf.iterrows():
        ax.plot([0, 1], [row["pre"], row["post"]], "-o", color=BLUE2,
                lw=0.9, ms=3, alpha=0.8)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["Pre", "Post"])
    ax.set_ylabel("T-cell fraction")
    ax.set_ylim(0, max(tf["post"].max() * 1.15, 0.18))
    annotate_stats(ax, "Median +0.090\n8/8 up\nWilcoxon p = 0.0078\nT per tumor +0.219", x=0.02, y=0.72, ha="left", va="top", fontsize=7.5)
    panel_letter(ax, "A")
    panel_title(ax, "T-cell fraction")
    style_ax(ax)

    ax = fig.add_subplot(gs[0, 1])
    td = q1d.pivot(index="patient_base", columns="timepoint",
                   values="med_tumor_to_T")
    for _, row in td.iterrows():
        ax.plot([0, 1], [row["pre"], row["post"]], "-o", color=BLUE2,
                lw=0.9, ms=3, alpha=0.8)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["Pre", "Post"])
    ax.set_ylabel("Tumor\u2192T dist. (\u00b5m)")
    annotate_stats(ax, "Median \u221268 \u00b5m\n8/8 down, p = 0.0078\nWithin 50 \u00b5m: +0.31", fontsize=7.5)
    panel_title(ax, "Tumor\u2192T distance")
    style_ax(ax)

    # (B) ligand zoning by tumor state
    ax = fig.add_subplot(gs[1, 0:2])
    # q2a is WIDE format: tier2 rows x VEGFA_pos/VCAN_pos/ANXA1_pos columns.
    states = ["AC-like", "MES-like 1", "MES-like 2", "NPC-like 1",
              "NPC-like 2", "OPC-like"]
    tier2_full = {"AC-like": "Tumor AC-like",
                  "MES-like 1": "Tumor MES-like 1",
                  "MES-like 2": "Tumor MES-like 2 (hypoxia)",
                  "NPC-like 1": "Tumor NPC-like 1",
                  "NPC-like 2": "Tumor NPC-like 2",
                  "OPC-like": "Tumor OPC-like"}
    short = ["AC", "MES1", "MES2", "NPC1", "NPC2", "OPC"]
    ligs = ["VEGFA", "VCAN", "ANXA1"]
    lcols = {"VEGFA": RED, "VCAN": TEAL, "ANXA1": VIOLET}
    x = np.arange(len(states))
    w = 0.26
    for i, lig in enumerate(ligs):
        col = lig + "_pos"
        vals = [float(q2a.loc[q2a["tier2"] == tier2_full[s], col].iloc[0])
                for s in states]
        ax.bar(x + (i - 1) * w, vals, w, color=lcols[lig], edgecolor="black",
               linewidth=0.6, label=lig)
    ax.set_xticks(x)
    ax.set_xticklabels(short)
    ax.set_ylabel("Ligand-positive fraction")
    ax.set_ylim(0, 1.05)
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.02), handlelength=1.0,
              fontsize=7, ncol=3)
    annotate_stats(ax, "VEGFA+ in MES-like 2: 92.2% vs 22\u201359% others (p = 4.5$\\times$10$^{-7}$)\nVEGFA+ higher mes2 in 16/16 (median +0.20, p = 3.1$\\times$10$^{-5}$)\nVCAN p = 0.27, ANXA1 p = 0.46", fontsize=7)
    panel_letter(ax, "B")
    panel_title(ax, "Ligand zoning by tumor state")
    style_ax(ax)

    # (C) hypoxic exclusion: paired distance + 50-um presence
    ax = fig.add_subplot(gs[2, 0])
    for _, row in t2_post.iterrows():
        ax.plot([0, 1], [row["dist_oth"], row["dist_hyp"]], "-o", color=BLUE2,
                lw=0.9, ms=3, alpha=0.8)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["Other states", "MES-like 2"])
    ax.set_ylabel("Median dist. to T (\u00b5m)")
    annotate_stats(ax, "Pooled medians\n90.7 vs 51.2 \u00b5m\nPaired diff +31.5 \u00b5m\n8/8, p = 0.0078", x=0.02, y=0.72, ha="left", va="top", fontsize=7.5)
    panel_letter(ax, "C")
    panel_title(ax, "Distance to nearest T")
    style_ax(ax)

    ax = fig.add_subplot(gs[2, 1])
    for _, row in t2_post.iterrows():
        ax.plot([0, 1], [row["f50_oth"], row["f50_hyp"]], "-o", color=BLUE2,
                lw=0.9, ms=3, alpha=0.8)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["Other states", "MES-like 2"])
    ax.set_ylabel("Fraction with T \u226450 \u00b5m")
    ax.set_ylim(0, 0.75)
    annotate_stats(ax, "0.281 vs 0.487\nDiff \u22120.216, p = 0.0156", fontsize=7.5)
    panel_title(ax, "T within 50 \u00b5m")
    style_ax(ax)

    # (D) method control: 4-gene MM2 proxy vs author mes2 split
    ax = fig.add_subplot(gs[3, 0:2])
    proxy_diff = proxy_post["dist_diff"].values
    mes2_diff = mes2_post["dist_diff"].values
    rng = np.random.default_rng(2)
    meds = []
    for i, (vals, label) in enumerate([(proxy_diff, "4-gene proxy"),
                                       (mes2_diff, "Author mes2")]):
        xj = i + rng.uniform(-0.12, 0.12, size=len(vals))
        ax.scatter(xj, vals, s=18, color=[GRAY, BLUE][i], alpha=0.8,
                   linewidths=0, zorder=2)
        med = float(np.median(vals))
        meds.append(med)
        ax.plot([i - 0.28, i + 0.28], [med, med], color="black", lw=2.0)
    ax.axhline(0, color="black", lw=0.8)
    ax.set_xticks([0, 1])
    ax.set_xticklabels([f"4-gene proxy\n(median {meds[0]:+.1f} \u00b5m)",
                        f"Author mes2\n(median {meds[1]:+.1f} \u00b5m)"], fontsize=8)
    ax.set_ylabel("Dist. diff (\u00b5m)")
    annotate_stats(ax, "Proxy vs author mes2: $\\rho$ = 0.158; proxy alone null (p = 0.46)\nCXCL10 in the proxy is a T-cell chemoattractant\n\u2014 cancels the signal",
                   x=0.5, y=0.97, ha="center", va="top", fontsize=7.5)
    panel_letter(ax, "D")
    panel_title(ax, "Method control: hypoxia proxy")
    style_ax(ax)

    save(fig, "fig5.png")


# ============================================================== FIGURE 6 ==
def fig6():
    trb = load_csv(GSE330, "trb_analysis/trb_diversity.csv")
    ov = load_csv(GSE330, "trb_analysis/trb_overlap.csv")
    xf = load_csv(GSE330, "xenium_analysis/xenium_vegfa_mm2_fixup.csv")
    xfr = load_csv(GSE330, "xenium_analysis/xenium_fractions.csv")

    fig = plt.figure(figsize=(7.0, 7.8), constrained_layout=True)
    gs = fig.add_gridspec(4, 2, height_ratios=[1, 1, 1, 1], hspace=0.25, wspace=0.35)

    # (A) cohort composition
    ax = fig.add_subplot(gs[0, 0:2])
    items = ["Xenium sections", "TRB samples", "Hamster TMA"]
    counts = [15, 16, 2]
    ax.barh(np.arange(len(items)), counts, height=0.55,
            color=[BLUE, TEAL, GRAY], edgecolor="black", linewidth=0.8)
    for i, v in enumerate(counts):
        ax.text(v + 0.3, i, str(v), va="center", fontsize=9)
    ax.set_yticks(np.arange(len(items)))
    ax.set_yticklabels(items, fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel("Count")
    ax.set_xlim(0, 19)
    panel_letter(ax, "A")
    panel_title(ax, "Cohort composition \u2014 no outcome labels in GEO")
    style_ax(ax)

    # (B) paired TRB: top-clone fraction + Shannon
    ax = fig.add_subplot(gs[1, 0])
    tp = trb.pivot(index="patient", columns="timepoint",
                   values="top_clone").dropna()
    for _, row in tp.iterrows():
        ax.plot([0, 1], [row["pre"], row["post"]], "-o", color=BLUE2,
                lw=0.9, ms=3, alpha=0.8)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["Pre", "Post"])
    ax.set_ylabel("Top-clone fraction")
    annotate_stats(ax, "Median +0.041\n5/5 up, p = 0.0625\n(n = 5 minimum)", x=0.02, y=0.72, ha="left", va="top", fontsize=7.5)
    panel_letter(ax, "B")
    panel_title(ax, "Top-clone fraction")
    style_ax(ax)

    ax = fig.add_subplot(gs[1, 1])
    sh = trb.pivot(index="patient", columns="timepoint",
                   values="shannon").dropna()
    for _, row in sh.iterrows():
        ax.plot([0, 1], [row["pre"], row["post"]], "-o", color=BLUE2,
                lw=0.9, ms=3, alpha=0.8)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["Pre", "Post"])
    ax.set_ylabel("Shannon diversity")
    annotate_stats(ax, "Median \u22120.81, 4/5 down (p = 0.31)\nClonality up 4/5; Morisita-Horn 0.18\n21.7% of post repertoire present pre", x=0.02, y=0.03, ha="left", va="bottom", fontsize=7.5)
    panel_title(ax, "Shannon diversity")
    style_ax(ax)

    # (C) VEGFA-MM2 replication: per-section diffs + MM2+ fractions
    ax = fig.add_subplot(gs[2, 0])
    diffs = xf["mm2_mean_diff"].values
    rng = np.random.default_rng(0)
    ax.scatter(rng.uniform(-0.2, 0.2, size=len(diffs)), diffs, s=16,
               color=BLUE2, alpha=0.8, linewidths=0, zorder=2)
    med = float(np.median(diffs))
    ax.plot([-0.4, 0.4], [med, med], color=RED, lw=1.6, ls="--")
    ax.axhline(0, color="black", lw=0.8)
    ax.set_xlim(-0.5, 0.5)
    ax.set_xticks([])
    ax.set_xlabel("Section (n = 14)")
    ax.set_ylabel("MM2 score diff (VEGFA+ \u2212 VEGFA\u2212)")
    annotate_stats(ax, "14/14 sections\nMedian diff +0.365\np = 0.00012", x=0.98, y=0.03, ha="right", va="bottom", fontsize=7.5)
    panel_letter(ax, "C")
    panel_title(ax, "VEGFA\u2013MM2 per section")
    style_ax(ax)

    ax = fig.add_subplot(gs[2, 1])
    pos = xf["mm2_pos_frac_pos"].values
    neg = xf["mm2_pos_frac_neg"].values
    x = np.arange(len(xf))
    ax.plot(x, neg, "o-", color=GRAY, lw=0.9, ms=3, label="VEGFA\u2212", alpha=0.8)
    ax.plot(x, pos, "o-", color=RED, lw=1.2, ms=3, label="VEGFA+")
    ax.set_xlabel("Section (n = 14)")
    ax.set_ylabel("MM2-positive fraction")
    ax.set_ylim(0, 1.0)
    annotate_stats(ax, "Median 62.8% vs 36.5%",
                   x=0.98, y=0.03, ha="right", va="bottom", fontsize=7.5)
    ax.legend(loc="lower left", handlelength=1.0, fontsize=7)
    panel_title(ax, "MM2-positive fraction")
    style_ax(ax)

    # (D) descriptive paired T / myeloid fractions (n = 2)
    # xenium_fractions.csv is per-section: aggregate per patient-timepoint
    # as a cell-count-weighted mean (this reproduces the legend values).
    xw = xfr.copy()
    xw["Tw"] = xw["frac_T"] * xw["n_cells"]
    xw["Mw"] = xw["frac_Myeloid"] * xw["n_cells"]
    agg = (xw.groupby(["patient", "timepoint"])
             .apply(lambda d: pd.Series(
                 {"frac_T": d["Tw"].sum() / d["n_cells"].sum(),
                  "frac_Myeloid": d["Mw"].sum() / d["n_cells"].sum()}),
                    include_groups=False)
             .reset_index())
    paired = agg[agg["patient"].isin(["P17", "P20"])]
    ax = fig.add_subplot(gs[3, 0])
    for _, row in paired.pivot(index="patient", columns="timepoint",
                               values="frac_T").iterrows():
        ax.plot([0, 1], [row["pre"], row["post"]], "-o",
                color=BLUE2, lw=1.2, ms=4)
    for pat, row in paired.pivot(index="patient", columns="timepoint",
                                 values="frac_T").iterrows():
        ax.text(1.05, row["post"], pat, fontsize=7, va="center", color=DARK)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["Pre", "Post"])
    ax.set_ylabel("T-cell fraction")
    ax.set_xlim(-0.1, 1.45)
    annotate_stats(ax, "P17 0.013\u21920.024\nP20 0.006\u21920.044",
                   x=0.02, y=0.85, ha="left", va="top", fontsize=7.5)
    panel_letter(ax, "D")
    panel_title(ax, "T-cell fraction (n = 2, descriptive)")
    style_ax(ax)

    ax = fig.add_subplot(gs[3, 1])
    for pat, row in paired.pivot(index="patient", columns="timepoint",
                                 values="frac_Myeloid").iterrows():
        ax.plot([0, 1], [row["pre"], row["post"]], "-o", color=TEAL,
                lw=1.2, ms=4)
        ax.text(1.05, row["post"], pat, fontsize=7, va="center", color=DARK)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["Pre", "Post"])
    ax.set_ylabel("Myeloid fraction")
    ax.set_xlim(-0.1, 1.45)
    annotate_stats(ax, "P17 0.089\u21920.148\nP20 0.070\u21920.167",
                   x=0.02, y=0.85, ha="left", va="top", fontsize=7.5)
    panel_title(ax, "Myeloid fraction (n = 2, descriptive)")
    style_ax(ax)

    save(fig, "fig6.png")


# ============================================================== FIGURE 7 ==
def fig7():
    gdon = load_csv(GBMAP, "gbmap_analysis/gbmap_per_donor.csv")
    z = np.load(os.path.join(GBMAP, "gbmap_analysis/gbmap_per_cell.npz"),
                allow_pickle=True)
    traj = z["traj"].astype(float)
    mod = z["mod_scores"].astype(float)
    lig = z["lig_vals"].astype(float)
    ist = z["is_tumor"].astype(bool)
    ismy = z["is_myeloid"].astype(bool)
    w = load_csv(WANG, "wang_analysis/wang_per_sample.csv")
    a1 = load_csv(WANG, "wang_analysis/wang_ap1.csv")

    # PC1 of the 12 module scores among tumor cells (GBmap)
    X = mod[ist]
    Xc = X - X.mean(axis=0)
    _, s, vt = np.linalg.svd(Xc, full_matrices=False)
    pc1 = Xc @ vt[0]
    pc1_var = float(s[0] ** 2 / (s ** 2).sum())
    rho_pc1 = stats.spearmanr(traj[ist], pc1)

    fig = plt.figure(figsize=(7.0, 7.8), constrained_layout=True)
    gs = fig.add_gridspec(4, 2, height_ratios=[1, 1, 1, 1], hspace=0.25, wspace=0.35)

    # (A) GBmap: TRAJ vs PC1 (per-cell hexbin, tumor cells only)
    ax = fig.add_subplot(gs[0, 0:2])
    hb = ax.hexbin(traj[ist], pc1, gridsize=60, cmap="Blues", mincnt=1,
                   linewidths=0)
    ax.set_xlabel("TRAJ")
    ax.set_ylabel("Module PC1")
    annotate_stats(ax, f"$\\rho$ = {rho_pc1.statistic:.3f}\n"
                       f"PC1 explains {pc1_var:.1%} of module variance\n"
                       f"n = 125,136 tumor cells\nS10 coverage 1,094/1,200",
                   x=0.98, y=0.97, ha="right", va="top", fontsize=7.5)
    panel_letter(ax, "A")
    panel_title(ax, "GBmap: TRAJ vs module PC1")
    style_ax(ax)

    # (B) per-donor TRAJ vs myeloid VEGFA + tumor VEGFA MM2 split
    ax = fig.add_subplot(gs[1, 0])
    ax.scatter(gdon["traj_mean"], gdon["VEGFA_myeloid"], s=16, color=BLUE2,
               alpha=0.75, linewidths=0, zorder=2)
    m, b = np.polyfit(gdon["traj_mean"], gdon["VEGFA_myeloid"], 1)
    xs = np.array([gdon["traj_mean"].min(),
                   gdon["traj_mean"].max()])
    ax.plot(xs, m * xs + b, color=RED, lw=1.2, ls="--")
    ax.set_xlabel("Donor tumor TRAJ (mean)")
    ax.set_ylabel("Myeloid VEGFA (mean)")
    annotate_stats(ax, "$\\rho$ = 0.334, p = 0.0079 (n = 62)\nVCAN p = 0.17, ANXA1 p = 0.088",
                   x=0.98, y=0.97, ha="right", va="top", fontsize=7.5)
    panel_letter(ax, "B")
    panel_title(ax, "Per-donor ligand vs TRAJ")
    style_ax(ax)

    ax = fig.add_subplot(gs[1, 1])
    mm2 = mod[ist, 2]
    veg = lig[ist, 0]
    hi = mm2 > np.median(mm2)
    means = [float(veg[hi].mean()), float(veg[~hi].mean())]
    fracs = [float((veg[hi] > 0).mean()), float((veg[~hi] > 0).mean())]
    bars = ax.bar([0, 1], means, width=0.55, color=[RED, BLUE2],
                  edgecolor="black", linewidth=0.8)
    for i, (mv, fv) in enumerate(zip(means, fracs)):
        ax.text(i, mv + 0.03, f"{mv:.3f}", ha="center", fontsize=8)
        ax.text(i, 0.06, f"{fv:.1%} positive", ha="center", fontsize=7,
                color="white", fontweight="bold")
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["MM2-high", "MM2-low"], fontsize=8)
    ax.set_ylabel("Tumor VEGFA (mean)")
    ax.set_ylim(0, 1.15)
    panel_title(ax, "Tumor VEGFA by MM2")
    style_ax(ax)

    # (C) Wang: paired TRAJ + per-module median diffs
    ax = fig.add_subplot(gs[2, 0])
    wp = w.pivot(index="pair", columns="timepoint", values="traj_median_tumor")
    for _, row in wp.iterrows():
        ax.plot([0, 1], [row["Primary"], row["Recurrent"]], "-o", color=BLUE2,
                lw=0.7, ms=2.5, alpha=0.7)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["Primary", "Recurrent"])
    ax.set_ylabel("TRAJ (median, tumor)")
    ax.set_ylim(-0.32, 0.28)
    annotate_stats(ax, "Median diff \u22120.053 (19 down / 10 up)\nWilcoxon p = 0.035, permutation p = 0.20", x=0.98, y=0.97, ha="right", va="top", fontsize=7.5)
    panel_letter(ax, "C")
    panel_title(ax, "Wang: recurrent vs primary")
    style_ax(ax)

    ax = fig.add_subplot(gs[2, 1])
    early_idx = [1, 4, 6, 7, 8, 10, 11]
    late_idx = [2, 3, 5]
    labels = []
    diffs = []
    colors = []
    for i in range(12):
        col = f"MM{i}_tumor"
        pp = w.pivot(index="pair", columns="timepoint", values=col)
        labels.append(f"MM{i}")
        diffs.append(float((pp["Recurrent"] - pp["Primary"]).median()))
        colors.append(BLUE if i in early_idx else (RED if i in late_idx
                                                   else GRAY))
    ax.barh(np.arange(12), diffs, height=0.65, color=colors,
            edgecolor="black", linewidth=0.6)
    ax.axvline(0, color="black", lw=0.8)
    ax.set_yticks(np.arange(12))
    ax.set_yticklabels(labels, fontsize=7)
    ax.invert_yaxis()
    ax.set_xlabel("Median diff (recurrent \u2212 primary)")
    ax.set_title("Per-module median diffs\n6 early modules up",
                 fontsize=8, pad=4, loc="left")
    style_ax(ax)

    # (D) AP-1 paired lines + ligand replication
    ax = fig.add_subplot(gs[3, 0])
    ap = a1.pivot(index="pair", columns="tp", values="AP1_core")
    for _, row in ap.iterrows():
        ax.plot([0, 1], [row["Primary"], row["Recurrent"]], "-o", color=BLUE2,
                lw=0.7, ms=2.5, alpha=0.7)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["Primary", "Recurrent"])
    ax.set_ylabel("AP-1 core score")
    ax.set_ylim(-0.05, 0.45)
    annotate_stats(ax, "Median diff +0.032 (17/29 up, p = 0.19)\nFOSL1 lower (p = 0.035)\nAP-1 vs TRAJ $\\rho$ = 0.10 (p = 0.60)",
                   x=0.98, y=0.97, ha="right", va="top", fontsize=7.5)
    panel_letter(ax, "D")
    panel_title(ax, "AP-1 in recurrent tumors")
    style_ax(ax)

    ax = fig.add_subplot(gs[3, 1])
    # legend values use ALL 58 samples (primary + recurrent pooled, complete cases)
    wcc = w.dropna(subset=["lig_VEGFA", "traj_median_tumor"])
    ax.scatter(wcc["traj_median_tumor"], wcc["lig_VEGFA"], s=14, color=RED,
               alpha=0.7, linewidths=0, label="VEGFA", zorder=2)
    wcc2 = w.dropna(subset=["lig_VCAN", "traj_median_tumor"])
    ax.scatter(wcc2["traj_median_tumor"], wcc2["lig_VCAN"], s=14, color=TEAL,
               alpha=0.7, linewidths=0, label="VCAN", zorder=2)
    ax.set_xlabel("Tumor TRAJ (median)")
    ax.set_ylabel("Ligand expression")
    annotate_stats(ax, "VEGFA $\\rho$ = 0.408, p = 0.0036\nVCAN $\\rho$ = 0.470, p = 0.00066\nn = 49 complete cases (of 58)",
                   x=0.98, y=0.03, ha="right", va="bottom", fontsize=7.5)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.28), handlelength=1.0,
              fontsize=7, ncol=2)
    panel_title(ax, "Ligand\u2013TRAJ replication")
    style_ax(ax)

    save(fig, "fig7.png")


if __name__ == "__main__":
    fig1()
    fig2()
    print("fig1, fig2 done")
    fig3()
    fig4()
    print("fig3, fig4 done")
    fig5()
    fig6()
    print("fig5, fig6 done")
    fig7()
    print("fig7 done")
