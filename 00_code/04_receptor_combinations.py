import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch
from matplotlib.ticker import StrMethodFormatter


DATA = Path("data/expression/Figure4D_RNA_only_plot_values.csv")
OUT = Path(os.environ.get("HONEYBEE_FIGURE_DIR", "figures"))
SPECIES = {"Honeybee": ("Amel", "#b874ad"), "Mosquito": ("Aaeg", "#a59fd5"),
           "Fruit fly": ("Dmel", "#7370a0")}
GENES = {"Ir76b": "#d81b60", "Ir25a": "#007f5f", "Ir8a": "#2691aa", "Orco": "#255580"}


def main():
    table = pd.read_csv(DATA)
    groups = table["combination"].drop_duplicates().tolist()
    table = table.set_index(["combination", "species"]).unstack("species").reindex(groups)
    if table.isna().any().any() or len(groups) != 11:
        raise ValueError("Expected complete RNA-only counts for 11 combinations.")
    counts = table["clusters"]
    if not np.allclose(table["cluster_pct"], 100 * counts / counts.sum()):
        raise ValueError("Cluster percentages do not match the counts.")

    plt.rcParams.update({"font.family": ["Arial", "DejaVu Sans"], "font.size": 7,
                         "pdf.fonttype": 42, "axes.linewidth": .5})
    fig = plt.figure(figsize=(4.33, 4.35))
    x = np.arange(len(groups))
    axes = [fig.add_axes([.145, y, .845, h]) for y, h in
            [(0.725, .205), (.635, .050), (.555, .060), (.482, .053), (.234, .205), (.058, .153)]]
    separators = np.flatnonzero(np.diff([g.count("+") for g in groups])) + .5
    for ax in axes:
        ax.set_xlim(-.53, len(x)-.45)
        ax.set_xticks([])
        ax.spines[["top", "right", "bottom"]].set_visible(False)
        ax.spines["left"].set_color("#abb5ba")
        for sep in separators:
            ax.axvline(sep, color="#c9d5da", lw=.5, dashes=(4, 4), zorder=1)
    panels = [
        (axes[0], "clusters", (0, 30), [0, 10, 20, 30]),
        (axes[1], "cells", (4500, 16000), [4500, 10000, 16000]),
        (axes[2], "cells", (200, 3000), [200, 1500, 3000]),
        (axes[3], "cells", (0, 100), [0, 50, 100]),
        (axes[4], "cluster_pct", (0, 50), list(range(0, 51, 10))),
    ]
    for ax, field, limits, ticks in panels:
        ax.set_ylim(limits)
        ax.set_yticks(ticks)
        ax.yaxis.set_major_formatter(StrMethodFormatter("{x:,.0f}"))
        ax.tick_params(length=3, width=.5, pad=3)
        ax.set_axisbelow(True)
        ax.yaxis.grid(color="#e6ebed", lw=.5)
        for i, (species, (_, color)) in enumerate(SPECIES.items()):
            xpos = x + (i - 1) * .245
            values = table[field][species]
            ax.bar(xpos, values, width=.245, color=color, zorder=2)
            for xx, v in zip(xpos, values):
                if (field == "clusters" and 0 < v <= 2) or (
                    field == "cells" and limits[0] == 0 and 0 < v <= 100
                ) or (field == "cluster_pct" and 0 < v < 5):
                    label = f"{v:.1f}" if field == "cluster_pct" else str(int(v))
                    shift = -.18 if field == "cluster_pct" and i == 0 else 0
                    ax.annotate(label, (xx + shift, v), xytext=(0, 1.5), textcoords="offset points",
                                color=color, fontsize=6, ha="center", va="bottom")

    # I keep the two breaks on the cell-count axis only.
    for a, b in [(axes[1], axes[2]), (axes[2], axes[3])]:
        for ax, y in [(a, 0), (b, 1)]:
            ax.plot([-.006, .006], [y-.045, y+.045], transform=ax.transAxes, clip_on=False,
                    color="#abb5ba", lw=.6)
    for y, label in [(0.827, "# of clusters"), (.586, "# of cells"), (.335, "Cluster %")]:
        fig.text(.024, y, label, ha="center", va="center", rotation=90)
    dots = axes[-1]
    dots.spines["left"].set_visible(False)
    dots.set_ylim(-.5, 3.5)
    dots.set_yticks([3, 2, 1, 0], list(GENES))
    dots.tick_params(length=0, pad=3)
    for (gene, color), y, tick in zip(GENES.items(), [3, 2, 1, 0], dots.get_yticklabels()):
        tick.set_color(color)
        tick.set_fontstyle("italic")
        present = np.array([gene in group.split("+") for group in groups])
        dots.scatter(x[~present], np.full((~present).sum(), y), s=21, facecolors="white",
                     edgecolors="#d4dde5", linewidths=.5)
        dots.scatter(x[present], np.full(present.sum(), y), s=21, color=color, linewidths=0)
    fig.text(.57, .022, "Co-receptor combination", ha="center")
    fig.legend(handles=[Patch(color=color, label=label) for label, color in SPECIES.values()],
               loc="upper center", bbox_to_anchor=(.58, .991), ncols=3, frameon=False)

    OUT.mkdir(parents=True, exist_ok=True)
    # Figure 4D
    fig.savefig(OUT / "Figure4D_RNA_only_combinations.pdf")
    plt.close(fig)


if __name__ == "__main__":
    main()
