"""Generate the two main figures and Supplementary Figure S1.

Usage
-----
python src/make_figures.py --results results --outdir figures
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyBboxPatch, Rectangle
from scipy import stats


COLORS = {"CK": "#8A8A8A", "DCD": "#3C78A8", "Biochar": "#D17C45"}
CSU_ORDER = [
    "Lawn", "PB", "PG", "SJ", "NZE", "GB×CH",
    "SJ×JN", "OF×RC", "AR×ND", "GYL×EF", "GYL×ND",
]


def configure_plotting() -> None:
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "DejaVu Sans", "sans-serif"],
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "font.size": 7,
            "axes.spines.right": False,
            "axes.spines.top": False,
            "axes.linewidth": 0.7,
            "xtick.major.width": 0.6,
            "ytick.major.width": 0.6,
            "legend.frameon": False,
        }
    )


def save_all(fig, stem: Path) -> None:
    fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(stem.with_suffix(".tiff"), dpi=600, bbox_inches="tight")
    fig.savefig(stem.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)


def figure_1(outdir: Path) -> None:
    fig = plt.figure(figsize=(7.2, 4.55))
    gs = fig.add_gridspec(
        2, 2, width_ratios=[0.9, 1.35], height_ratios=[0.92, 1.08],
        hspace=0.25, wspace=0.22
    )
    ax_a = fig.add_subplot(gs[:, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, 1])
    for ax in (ax_a, ax_b, ax_c):
        ax.set_axis_off()
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)

    ax_a.text(0.0, 1.02, "a", transform=ax_a.transAxes, fontsize=9, fontweight="bold")
    ax_a.text(
        0.06, 0.94, "One existing configuration-site unit (CSU)",
        fontsize=8.5, fontweight="bold", va="top"
    )
    ax_a.text(
        0.06, 0.885, "80 m × 5 m strip; 11 CSUs in total",
        fontsize=7.2, color="#555555", va="top"
    )
    ax_a.add_patch(Rectangle((0.08, 0.13), 0.84, 0.66, fill=False, ec="#333333", lw=1.0))
    orders = [
        ["Biochar", "DCD", "CK"],
        ["DCD", "CK", "Biochar"],
        ["CK", "Biochar", "DCD"],
    ]
    for i, order in enumerate(orders):
        x0 = 0.105 + i * 0.275
        ax_a.text(x0 + 0.105, 0.815, f"Block {i + 1}", ha="center", fontsize=7.2, fontweight="bold")
        ax_a.add_patch(
            FancyBboxPatch(
                (x0, 0.16), 0.21, 0.58,
                boxstyle="round,pad=0.006,rounding_size=0.01",
                fc="#F5F5F5", ec="#B0B0B0", lw=0.8
            )
        )
        for j, treatment in enumerate(order):
            y0 = 0.57 - j * 0.18
            ax_a.add_patch(
                Rectangle((x0 + 0.035, y0), 0.14, 0.11, fc=COLORS[treatment], ec="white", lw=0.8)
            )
            label = "Bio-\nchar" if treatment == "Biochar" else treatment
            ax_a.text(
                x0 + 0.105, y0 + 0.055, label, ha="center", va="center",
                color="white", fontsize=6.8, fontweight="bold", linespacing=0.85
            )
    ax_a.text(
        0.08, 0.07,
        "One CK, one DCD and one biochar plot per block;\n"
        "order randomized independently within blocks",
        fontsize=6.7, color="#555555"
    )

    ax_b.text(0.0, 1.02, "b", transform=ax_b.transAxes, fontsize=9, fontweight="bold")
    ax_b.text(0.05, 0.88, "Field timeline and repeated measurements", fontsize=8.5, fontweight="bold")
    ax_b.plot([0.12, 0.88], [0.48, 0.48], color="#888888", lw=1.4)
    timeline = [
        (0.18, "Treatment application", "6–10 May 2023", "top"),
        (0.52, "July window", "18–19 July 2023", "bottom"),
        (0.84, "October window", "19–20 October 2023", "top"),
    ]
    for x, heading, subheading, vertical in timeline:
        ax_b.scatter([x], [0.48], s=70, color="#185A9D", zorder=3)
        y = 0.34 if vertical == "top" else 0.60
        ax_b.text(x, y, heading + "\n" + subheading, ha="center", va=vertical, fontsize=6.9)
    ax_b.text(
        0.5, 0.08, "99 permanent collars × 2 windows = 198 original flux observations",
        ha="center", fontsize=7.6, fontweight="bold"
    )

    ax_c.text(0.0, 1.02, "c", transform=ax_c.transAxes, fontsize=9, fontweight="bold")
    ax_c.text(0.05, 0.90, "Paired effect construction and inference boundary", fontsize=8.5, fontweight="bold")
    boxes = [
        (0.06, 0.56, 0.24, 0.18, "CSU × block × window\nrandomization stratum"),
        (0.38, 0.56, 0.24, 0.18, "log(Treatment / CK)\nDCD and biochar"),
        (0.70, 0.56, 0.24, 0.18, "132 paired effects\n4 window × treatment means"),
    ]
    for x, y, width, height, text in boxes:
        ax_c.add_patch(
            FancyBboxPatch(
                (x, y), width, height,
                boxstyle="round,pad=0.008,rounding_size=0.015",
                fc="#EAF1F8", ec="#6389AD", lw=0.8
            )
        )
        ax_c.text(x + width / 2, y + height / 2, text, ha="center", va="center", fontsize=6.7)
    for x1, x2 in [(0.30, 0.38), (0.62, 0.70)]:
        ax_c.annotate(
            "", xy=(x2, 0.65), xytext=(x1, 0.65),
            arrowprops=dict(arrowstyle="-|>", lw=0.9, color="#555555")
        )
    ax_c.text(0.06, 0.38, "Primary inference", fontsize=7.2, fontweight="bold")
    ax_c.text(0.34, 0.38, "CSU-clustered covariance (11 clusters; t, df = 10)", fontsize=6.8)
    ax_c.text(0.06, 0.28, "Sensitivity", fontsize=7.2, fontweight="bold")
    ax_c.text(0.34, 0.28, "Randomized-block clustering (33 clusters; t, df = 32)", fontsize=6.8)
    ax_c.text(0.06, 0.15, "Inference boundary", fontsize=7.2, fontweight="bold")
    ax_c.text(
        0.36, 0.15,
        "CSUs are existing contexts, not randomized plant treatments;\n"
        "CSU-specific estimates are descriptive and are not ranked by significance.",
        fontsize=6.8, color="#8B3A3A", va="center"
    )
    fig.subplots_adjust(left=0.03, right=0.99, top=0.97, bottom=0.04)
    save_all(fig, outdir / "Figure1_revised_design_and_analysis")


def add_errorbar(ax, x, lower, upper, y, color, size=4.0) -> None:
    ax.errorbar(
        x, y, xerr=np.array([[x - lower], [upper - x]]),
        fmt="o", markersize=size, markerfacecolor=color, markeredgecolor="white",
        markeredgewidth=0.45, ecolor=color, elinewidth=1.0, capsize=2.0,
        capthick=0.8, zorder=3
    )


def figure_2(results: Path, outdir: Path) -> None:
    primary = pd.read_csv(results / "primary_csu_clustered_effects.csv")
    csu = pd.read_csv(results / "exploratory_csu_effects.csv")
    csu["CSU"] = pd.Categorical(csu["CSU"], CSU_ORDER, ordered=True)

    fig = plt.figure(figsize=(7.20, 4.35), constrained_layout=False)
    gs = fig.add_gridspec(1, 3, width_ratios=[0.92, 1.18, 1.18], wspace=0.36)
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[0, 2], sharey=ax_b)
    rows = [
        ("July", "DCD", "July · DCD"),
        ("July", "Biochar", "July · Biochar"),
        ("October", "DCD", "October · DCD"),
        ("October", "Biochar", "October · Biochar"),
    ]
    y_values = np.arange(len(rows))[::-1]
    for y, (window, treatment, label) in zip(y_values, rows):
        row = primary[
            (primary["Window"] == window) & (primary["Treatment"] == treatment)
        ].iloc[0]
        add_errorbar(
            ax_a, row.percent_change, row.ci_low_percent, row.ci_high_percent,
            y, COLORS[treatment], size=5
        )
    ax_a.axvline(0, color="#777777", lw=0.8, ls="--", zorder=0)
    ax_a.set_yticks(y_values, [row[2] for row in rows])
    ax_a.set_xlim(-32, 3)
    ax_a.set_xlabel("Change relative to CK (%)")
    ax_a.set_title("Corridor-level effects", fontsize=8, pad=6)
    ax_a.grid(axis="x", color="#E6E6E6", lw=0.55, zorder=0)

    for ax, treatment in [(ax_b, "DCD"), (ax_c, "Biochar")]:
        subset = csu[csu["Treatment"] == treatment].sort_values("CSU")
        y_values = np.arange(len(subset))[::-1]
        for y, (_, row) in zip(y_values, subset.iterrows()):
            add_errorbar(
                ax, row.percent_change, row.ci_low_percent, row.ci_high_percent,
                y, COLORS[treatment]
            )
        ax.axvline(0, color="#777777", lw=0.8, ls="--", zorder=0)
        ax.set_xlim(-45, 25)
        ax.set_xticks([-40, -20, 0, 20])
        ax.set_xlabel("Change relative to CK (%)")
        ax.set_title(f"{treatment}: exploratory CSU effects", fontsize=8, pad=6)
        ax.grid(axis="x", color="#E6E6E6", lw=0.55, zorder=0)
    subset = csu[csu["Treatment"] == "DCD"].sort_values("CSU")
    ax_b.set_yticks(np.arange(len(subset))[::-1], subset["CSU"].astype(str))
    ax_b.set_ylabel("Existing configuration-site unit (CSU)")
    ax_c.tick_params(labelleft=False)
    for label, ax in zip(["a", "b", "c"], [ax_a, ax_b, ax_c]):
        ax.text(-0.12, 1.04, label, transform=ax.transAxes, fontsize=9, fontweight="bold", va="bottom")
    fig.text(
        0.5, 0.015,
        "a, CSU-clustered 95% confidence intervals (11 CSUs). "
        "b–c, descriptive t intervals from three block-level summaries per CSU; "
        "no CSU-wise significance tests.",
        ha="center", va="bottom", fontsize=6.4
    )
    fig.subplots_adjust(left=0.115, right=0.985, top=0.90, bottom=0.19)
    save_all(fig, outdir / "Figure2_original_data_paired_effects")


def figure_s1(results: Path, outdir: Path) -> None:
    data = pd.read_csv(results / "paired_effects_132.csv")
    means = data.groupby(["Window", "Treatment"], observed=True)["LogRatio"].mean().rename("Fitted")
    data = data.join(means, on=["Window", "Treatment"])
    data["Residual"] = data["LogRatio"] - data["Fitted"]
    rng = np.random.default_rng(20260619)

    fig, axes = plt.subplots(
        1, 3, figsize=(7.2, 2.45), gridspec_kw={"width_ratios": [1.0, 1.0, 1.35]}
    )
    ax = axes[0]
    for treatment in ["DCD", "Biochar"]:
        subset = data[data["Treatment"] == treatment]
        jitter = rng.normal(0, 0.004, len(subset))
        ax.scatter(
            subset["Fitted"] + jitter, subset["Residual"], s=13, alpha=0.62,
            color=COLORS[treatment], edgecolor="none", label=treatment
        )
    ax.axhline(0, color="#777777", lw=0.8, ls="--")
    ax.set_xlabel("Fitted log ratio")
    ax.set_ylabel("Residual")
    ax.set_title("Residuals vs fitted", fontsize=8)
    ax.legend(loc="upper left")

    ax = axes[1]
    theoretical, observed = stats.probplot(data["Residual"], dist="norm", fit=False)
    slope, intercept, _ = stats.probplot(data["Residual"], dist="norm", fit=True)[1]
    ax.scatter(theoretical, observed, s=12, color="#4C78A8", alpha=0.68, edgecolor="none")
    xline = np.array([min(theoretical), max(theoretical)])
    ax.plot(xline, intercept + slope * xline, color="#777777", lw=0.9)
    ax.set_xlabel("Theoretical quantile")
    ax.set_ylabel("Residual quantile")
    ax.set_title("Normal Q–Q", fontsize=8)

    ax = axes[2]
    groups = [("July", "DCD"), ("July", "Biochar"), ("October", "DCD"), ("October", "Biochar")]
    positions = np.arange(4)
    values = [
        data[(data["Window"] == window) & (data["Treatment"] == treatment)]["LogRatio"].values
        for window, treatment in groups
    ]
    boxplot = ax.boxplot(
        values, positions=positions, widths=0.55, patch_artist=True, showfliers=False,
        medianprops={"color": "#333333", "linewidth": 1.0},
        whiskerprops={"color": "#666666"}, capprops={"color": "#666666"}
    )
    for patch, (_, treatment) in zip(boxplot["boxes"], groups):
        patch.set_facecolor(COLORS[treatment])
        patch.set_alpha(0.45)
        patch.set_edgecolor(COLORS[treatment])
    for x, ((window, treatment), group_values) in enumerate(zip(groups, values)):
        ax.scatter(
            np.full(len(group_values), x) + rng.uniform(-0.18, 0.18, len(group_values)),
            group_values, s=10, color=COLORS[treatment], alpha=0.58, edgecolor="none"
        )
    ax.axhline(0, color="#777777", lw=0.8, ls="--")
    ax.set_xticks(positions, ["Jul\nDCD", "Jul\nBiochar", "Oct\nDCD", "Oct\nBiochar"])
    ax.set_ylabel("log(Treatment / CK)")
    ax.set_title("Paired-effect distributions", fontsize=8)
    for label, ax in zip(["a", "b", "c"], axes):
        ax.text(-0.15, 1.04, label, transform=ax.transAxes, fontsize=9, fontweight="bold")
    fig.tight_layout(pad=1.0, w_pad=1.3)
    save_all(fig, outdir / "FigureS1_paired_model_diagnostics")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--results", type=Path, default=Path("results"))
    parser.add_argument("--outdir", type=Path, default=Path("figures"))
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    configure_plotting()
    figure_1(args.outdir)
    figure_2(args.results, args.outdir)
    figure_s1(args.results, args.outdir)


if __name__ == "__main__":
    main()
