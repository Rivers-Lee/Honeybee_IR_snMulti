"""Draw the Figure S3 residue-contact summaries from Data S1 tables."""

import csv
import math
import os
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
from matplotlib import pyplot as plt
from matplotlib.colors import LinearSegmentedColormap


DATA = Path("data/contacts")
SPECIES = ("Amel", "Aste", "Rpro")
ROLES = ("T1", "T2", "C1", "C2")
COLORS = ("#a972a6", "#9696ca", "#739890")
FORMATS = ("pdf", "svg", "png")


def load_matrix(path):
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, [])
        rows = list(reader)
    if (header != ["", *ROLES] or any(len(row) != 5 for row in rows)
            or tuple(row[0] for row in rows) != ROLES):
        raise ValueError(f"{path}: expected a 4 x 4 T1/T2/C1/C2 contact matrix")
    values = np.full((4, 4), np.nan)
    for i, row in enumerate(rows):
        for j, cell in enumerate(row[1:]):
            if i == j:
                if cell.strip():
                    raise ValueError(f"{path}: expected a blank diagonal")
            else:
                value = float(cell)
                if value < 0 or not value.is_integer():
                    raise ValueError(f"{path}: expected a nonnegative integer count")
                values[i, j] = value
    if not np.array_equal(values, values.T, equal_nan=True):
        raise ValueError(f"{path}: reciprocal matrix entries disagree")
    return values


def load_bars(path, aste_count):
    with path.open(newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    by_species = {row["species"]: row for row in rows}
    if len(rows) != len(SPECIES) or set(by_species) != set(SPECIES):
        raise ValueError(f"{path}: expected one row each for Amel, Aste and Rpro")
    bars = []
    for species in SPECIES:
        row = by_species[species]
        value = float(row["value"])
        if not math.isfinite(value) or value < 0:
            raise ValueError(f"{path}: invalid {species} count")
        sd = float(row["sd"]) if row["sd"].strip() else None
        if species == "Aste":
            if sd is not None or value != aste_count:
                raise ValueError(f"{path}: 43EF must match its matrix count and have no SD")
        elif sd is None or not math.isfinite(sd) or sd < 0:
            raise ValueError(f"{path}: {species} requires a nonnegative seed-mean SD")
        bars.append((value, sd))
    return bars


def draw_matrices(matrices):
    width, height = 366, 124
    fig = plt.figure(figsize=(width / 72, height / 72))
    fig.text(3.5 / width, 1 - 3 / height, "A", fontsize=10, va="top")
    limit = max(1, max(np.nanmax(matrix) for matrix in matrices))
    cmap = LinearSegmentedColormap.from_list(
        "contacts", ["#fff3ee", "#fbd7d2", "#f7b9b8", "#f395aa", "#ee6f9d",
                     "#cf4795", "#ac2b85", "#822b79", "#48305e"]
    ).with_extremes(bad="white")
    for x, matrix in zip((14, 124.8, 235.6), matrices):
        ax = fig.add_axes((x / width, 15 / height, 95.5 / width, 95.5 / height))
        image = ax.pcolormesh(np.ma.masked_invalid(matrix), cmap=cmap, vmin=0, vmax=limit,
                              edgecolors="face", linewidth=.12)
        ax.set(xlim=(0, 4), ylim=(4, 0), aspect="equal")
        ax.set_xticks(np.arange(4) + .5, ROLES)
        ax.set_yticks(np.arange(4) + .5, ROLES)
        ax.tick_params(length=0, pad=3.6, labelsize=6)
        for spine in ax.spines.values():
            spine.set_visible(False)
        for (row, col), value in np.ndenumerate(matrix):
            if np.isfinite(value):
                color = "white" if value >= .55 * limit else "#231f20"
                ax.text(col + .5, row + .5, f"{value:.0f}", ha="center", va="center",
                        color=color, fontsize=6)
    cax = fig.add_axes((334.6 / width, 15 / height, 4.8 / width, 95.5 / height))
    cbar = fig.colorbar(image, cax=cax)
    cbar.set_ticks(np.arange(0, 1 + int(limit), 50) if 250 <= limit <= 350 else
                   np.linspace(0, limit, 6))
    cbar.ax.tick_params(labelsize=6, length=3, width=.7, pad=3)
    cbar.set_label("contact residue pairs", fontsize=7, labelpad=3)
    cbar.solids.set_rasterized(False)
    cbar.solids.set_edgecolor("face")
    return fig


def draw_bars(bars):
    width, height = 131, 129
    fig = plt.figure(figsize=(width / 72, height / 72))
    ax = fig.add_axes((28 / width, 24 / height, 100 / width, 98.5 / height))
    fig.text(1.5 / width, 1 - 3 / height, "B", fontsize=10, va="top")
    values, sds = zip(*bars)
    upper = max(174.2, 1.2 * max(v + (sd or 0) for v, sd in bars))
    ax.set(xlim=(-.402, 2.402), ylim=(0, upper), ylabel="C1-C2 residue pairs")
    ax.yaxis.labelpad = 2
    ax.set_yticks(np.arange(0, upper, 20 if upper < 200 else 50))
    ax.set_xticks(range(3), ("Amel\nAF3", "Aste\n43EF", "Rpro\nAF3"))
    ax.tick_params(axis="x", length=2.5, pad=2.5, width=.8)
    ax.tick_params(axis="y", length=2.5, pad=2, width=.8, labelsize=6)
    ax.bar(range(3), values, width=.55, color=COLORS)
    for i, (value, sd) in enumerate(zip(values, sds)):
        if sd is not None:
            ax.errorbar(i, value, yerr=sd, fmt="none", color="#66686a",
                        elinewidth=.5, capsize=1.45, capthick=.5)
            offset = 1.2
        else:
            ax.annotate("cryo-EM", (i, value), xytext=(0, 1.2),
                        textcoords="offset points", ha="center", va="bottom", fontsize=5)
            offset = 6.8
        ax.annotate(f"{value:.1f}", (i, value + (sd or 0)), xytext=(0, offset),
                    textcoords="offset points", ha="center", va="bottom", fontsize=6)
    return fig


def main():
    matrices = [load_matrix(DATA / f"{species}_contact_matrix.tsv") for species in SPECIES]
    bars = load_bars(DATA / "C1_C2_contact_plot_values.tsv", matrices[1][2, 3])
    plt.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"],
        "font.size": 7, "axes.spines.top": False, "axes.spines.right": False,
        "axes.linewidth": .5, "pdf.fonttype": 42, "svg.fonttype": "none",
    })
    out = Path(os.environ.get("HONEYBEE_FIGURE_DIR") or "figures")
    out.mkdir(parents=True, exist_ok=True)

    fig = draw_matrices(matrices)
    for ext in FORMATS:
        # Figure S3A
        fig.savefig(out / f"FigureS3A_residue_contacts.{ext}", dpi=300)
    plt.close(fig)

    fig = draw_bars(bars)
    for ext in FORMATS:
        # Figure S3B
        fig.savefig(out / f"FigureS3B_co_receptor_contacts.{ext}", dpi=300)
    plt.close(fig)


if __name__ == "__main__":
    main()
