import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


DATA = Path("data/structures")
OUT = Path(os.environ.get("HONEYBEE_FIGURE_DIR", "figures"))
REFS = ("43EF_cryoEM_IR75k2_IR8a2", "RproIR75b2_IR8a2")
FIG3 = ("AmelIR75g2_IR76b2_AF3", *REFS)
FIG5 = ("ir75g_ir8a_22", "ir75g_ir25a_22", "ir75g_ir76b_ir8a_112",
        "ir75g_ir76b_ir25a_112", "ir75g_ir76b_ir8a_ir25a_1111", *REFS)
LABELS3 = ("Amel\nIR75g2-IR76b2", "Aste\nIR75k2-IR8a2", "Rpro\nIR75b2-IR8a2")
ROWS5 = ("M1\nIR75g2-IR8a2", "M2\nIR75g2-IR25a2", "M3\nIR75g1-IR76b1-IR8a2",
         "M4\nIR75g1-IR76b1-IR25a2", "M5\nIR75g1-IR76b1-IR8a1-IR25a1",
         "Aste  IR75k2-IR8a2\n43EF cryo-EM", "Rpro  IR75b2-IR8a2\nAF3")
COLS5 = ("M1", "M2", "M3", "M4", "M5", "Aste", "Rpro")


def read_matrix(name, models, diagonal):
    frame = pd.read_csv(DATA / name, sep="\t", index_col=0)
    if tuple(frame.index) != models or tuple(frame.columns) != models:
        raise ValueError(f"Incorrect model order in {name}.")
    values = frame.to_numpy(float)
    if (not np.isfinite(values).all() or not np.allclose(values, values.T) or
            not np.allclose(np.diag(values), diagonal) or np.any(values < 0) or
            (diagonal == 1 and np.any(values > 1))):
        raise ValueError(f"Invalid symmetric matrix in {name}.")
    return values


def heatmap(fig, ax, values, xlabels, ylabels, title, cmap_name, limits, ticks, gray_self=False):
    colors = plt.get_cmap(cmap_name).with_extremes(bad="#eeeeee")
    n = len(values)
    display = np.ma.array(values, mask=np.eye(n, dtype=bool)) if gray_self else values
    positions = np.arange(n+1) - .5
    image = ax.pcolormesh(positions, positions, display, cmap=colors, vmin=limits[0],
                          vmax=limits[1], edgecolors="white", linewidth=.5)
    ax.set(xlim=(-.5, n-.5), ylim=(n-.5, -.5), aspect="equal", title=title)
    ax.set_xticks(range(n), xlabels)
    ax.set_yticks(range(n), ylabels)
    ax.tick_params(length=0, pad=3)
    ax.spines[:].set_visible(False)
    for y, x in np.ndindex(values.shape):
        rgb = colors(image.norm(values[y, x]))[:3]
        text = "#222222" if gray_self and x == y or np.dot(rgb, [.2126, .7152, .0722]) > .53 else "white"
        ax.text(x, y, f"{values[y, x]:.2f}", ha="center", va="center", color=text)
    colorbar = fig.colorbar(image, ax=ax, fraction=.052, pad=.045, shrink=.82, ticks=ticks)
    colorbar.outline.set_linewidth(.5)
    colorbar.ax.tick_params(width=.5, length=2)
    colorbar.solids.set_rasterized(False)


def main():
    plt.rcParams.update({"font.family": ["Arial", "DejaVu Sans"], "font.size": 6.5,
                         "pdf.fonttype": 42, "axes.linewidth": .5})
    rmsd = read_matrix("Figure3E_rmsd_A.tsv", FIG3, 0)
    tm3 = read_matrix("Figure3E_tm_score.tsv", FIG3, 1)
    tm5 = read_matrix("Figure5B_tm_score.tsv", FIG5, 1)
    OUT.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(2, 1, figsize=(2.72, 3.58), layout="constrained")
    heatmap(fig, axes[0], rmsd, ("Amel", "Aste", "Rpro"), LABELS3, "C\u03b1 RMSD (\u00c5)", "YlOrRd",
            (0, rmsd.max()), [0, 2, 4, 6, 8])
    heatmap(fig, axes[1], tm3, ("Amel", "Aste", "Rpro"), LABELS3, "TM-score", "GnBu",
            (0, 1), np.linspace(0, 1, 6), gray_self=True)
    # Figure 3E
    fig.savefig(OUT / "Figure3E_structural_similarity.pdf")
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(4.65, 2.8), layout="constrained")
    heatmap(fig, ax, tm5, COLS5, ROWS5, "TM-score", "GnBu", (.60, .95),
            np.arange(.60, .951, .05), gray_self=True)
    # Figure 5B
    fig.savefig(OUT / "Figure5B_tm_score.pdf")
    plt.close(fig)


if __name__ == "__main__":
    main()
