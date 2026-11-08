import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.offsetbox import AnnotationBbox, DrawingArea
from matplotlib.patches import Patch, Wedge


DATA = Path("data/docking")
OUT = Path(os.environ.get("HONEYBEE_FIGURE_DIR", "figures"))
MODELS = ("ir75g_ir8a_22", "ir75g_ir25a_22", "ir75g_ir76b_ir8a_112",
          "ir75g_ir76b_ir25a_112", "ir75g_ir76b_ir8a_ir25a_1111")
COLORS = {"ir75g": "#6a3d9a", "ir8a": "#2691aa", "ir76b": "#d81b60", "ir25a": "#007f5f"}


def composition(composition_id):
    *roles, ratio = composition_id.split("_")
    return dict(zip(roles, map(int, ratio), strict=True))


def label(composition_id):
    return "-".join(f"IR{role[2:]}{n}" for role, n in composition(composition_id).items())


def diameter(delta_g):
    return 1 + 2 * np.sqrt((36 + 3 * max(0, -delta_g)) / np.pi)


def tetramer(composition_id, size):
    counts = composition(composition_id)
    roles = [r for r in COLORS for _ in range(counts.get(r, 0))]
    repeated = next((r for r in COLORS if counts.get(r) == 2), None)
    if repeated:
        other = [r for r in roles if r != repeated]
        roles = [repeated, other[0], repeated, other[1]]
    glyph = DrawingArea(size, size)
    for role, (start, stop) in zip(roles, [(90, 180), (0, 90), (270, 360), (180, 270)], strict=True):
        glyph.add_artist(Wedge((size/2, size/2), size/2-.5, start, stop,
                               facecolor=COLORS[role], edgecolor="#303030", lw=.5))
    return glyph


def composition_plot():
    table = pd.read_csv(DATA / "composition_summary.tsv", sep="\t", index_col="composition")
    table = table.loc[list(MODELS)]
    fields = ["mean_iptm", "mean_pae", "mean_delta_g"]
    if (not table.index.is_unique or not np.isfinite(table[fields]).all().all() or
            not table[["n_models", "n_delta_g"]].eq(25).all().all() or not table["n_seeds"].eq(5).all()):
        raise ValueError("Expected five complete 25-model composition summaries.")
    fig, ax = plt.subplots(figsize=(4.55, 2.65))
    fig.subplots_adjust(left=.14, right=.985, bottom=.18, top=.79)
    ax.set(xlim=(.612, .752), ylim=(25.25, 19.75), xlabel="Mean ipTM",
           ylabel="Mean inter-chain PAE (\u00c5)")
    ax.set_xticks(np.arange(.62, .75, .02))
    ax.set_yticks(range(20, 26))
    ax.grid(axis="y", color="#e3e5e7", lw=.5)
    ax.set_axisbelow(True)
    ax.spines[["right", "top"]].set_visible(False)
    labels = [(0, -13, "center"), (-6, 13, "right"), (10, 10, "left"),
              (-9, 12, "center"), (0, -19, "center")]
    for row, (dx, dy, align) in zip(table.itertuples(), labels, strict=True):
        xy = (row.mean_iptm, row.mean_pae)
        ax.add_artist(AnnotationBbox(tetramer(row.Index, diameter(row.mean_delta_g)), xy,
                                     frameon=False, pad=0, zorder=3))
        ax.annotate(label(row.Index), xy, xytext=(dx, dy), textcoords="offset points",
                    ha=align, va="center", fontsize=6.5)
    subunits = [Patch(facecolor=c, label="IR"+r[2:]) for r, c in COLORS.items()]
    fig.text(.325, .978, "Subunits", ha="center", va="top")
    fig.legend(handles=subunits, loc="upper left", bbox_to_anchor=(.12, .928), ncols=4,
               frameon=False, handlelength=.8, columnspacing=.8)
    sizes = [Line2D([], [], marker="o", linestyle="none", markerfacecolor="#e6e8eb",
                    markeredgecolor="#555555", markeredgewidth=.5, markersize=diameter(v)-1,
                    label=str(v)) for v in (-30, -60, -90)]
    fig.text(.835, .978, "Mean \u0394G (kcal/mol)", ha="center", va="top")
    fig.legend(handles=sizes, loc="upper right", bbox_to_anchor=(.988, .909), ncols=3,
               frameon=False, handleheight=2.3, handletextpad=.7, columnspacing=.8)
    # Figure 5A
    fig.savefig(OUT / "Figure5A_composition_means.pdf")
    plt.close(fig)


def ligand_plot():
    scores = pd.read_csv(DATA / "N30_display_order_matrix.tsv", sep="\t", index_col=0)
    names = pd.read_csv(DATA / "ligand_names.tsv", sep="\t", index_col="ligand_id")["name"]
    if (scores.shape != (5, 30) or set(scores.index) != set(MODELS) or not scores.columns.is_unique
            or not scores.index.is_unique or not names.index.is_unique
            or not {"HB011", "HB012", "HB013"} <= set(scores.columns)):
        raise ValueError("Expected the five compositions and the final 30-ligand order.")
    values = scores.to_numpy(float)
    if np.isinf(values).any():
        raise ValueError("Infinite ligand score.")
    labels = names.loc[scores.columns].copy()
    labels.loc["IR004"] += " [IR004]"
    fig, ax = plt.subplots(figsize=(7.2, 3.0), layout="constrained")
    image = ax.pcolormesh(values, cmap=plt.get_cmap("RdBu_r").with_extremes(bad="#d2d5d8"),
                          vmin=2.4, vmax=3.6, edgecolors="white", lw=.5)
    ax.invert_yaxis()
    ax.set_yticks(np.arange(5)+.5, [label(name) for name in scores.index])
    ax.set_xticks(np.arange(30)+.5, labels, rotation=90)
    ax.tick_params(length=0, pad=3)
    ax.spines[:].set_visible(False)
    for tick, ligand in zip(ax.get_xticklabels(), scores.columns):
        if ligand in {"HB011", "HB012", "HB013"}:
            tick.set_color("#9b5d18")
    colorbar = fig.colorbar(image, ax=ax, fraction=.028, pad=.024)
    colorbar.set_label("Size-adjusted Vina score")
    colorbar.set_ticks(np.arange(2.4, 3.61, .2))
    colorbar.outline.set_linewidth(.5)
    colorbar.solids.set_rasterized(False)
    # Figure 5E
    fig.savefig(OUT / "Figure5E_size_adjusted_vina_score.pdf")
    plt.close(fig)


def main():
    plt.rcParams.update({"font.family": ["Arial", "DejaVu Sans"], "font.size": 6.5,
                         "pdf.fonttype": 42, "axes.linewidth": .5})
    OUT.mkdir(parents=True, exist_ok=True)
    composition_plot()
    ligand_plot()


if __name__ == "__main__":
    main()
