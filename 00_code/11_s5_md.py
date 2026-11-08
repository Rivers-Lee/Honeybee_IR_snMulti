"""Draw only the quantitative MD plots in Figure S5F-G."""

import argparse
import csv
import math
import os
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
from matplotlib import pyplot as plt
from matplotlib.colors import LinearSegmentedColormap


DATA = Path("data/md")
FORMATS = ("pdf", "svg", "png")
CHAINS = {"A": "IR75g", "B": "IR75g", "C": "IR8a", "D": "IR8a"}
PAIRS = ("A-C", "A-D", "B-C", "B-D")
SNAPSHOTS = {0: 0, 500: 50, 1000: 100}
CORE = (-16.632050700480907, 16.886418355397844)
PURPLE, CYAN, DARK = "#6a3d9a", "#2691aa", "#344b47"


def read_tsv(path, columns):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        missing = set(columns) - set(reader.fieldnames or ())
        if missing:
            raise ValueError(f"{path}: missing columns {', '.join(sorted(missing))}")
        yield from reader


def load_rmsd(path):
    series = {chain: [] for chain in CHAINS}
    metric = "TM_assembly_fit_rmsd_nm"
    for row in read_tsv(path, ("time_ns", "chain", "protein", "domain", metric)):
        if row["domain"] != "TM_helices":
            continue
        chain = row["chain"]
        if chain not in CHAINS or row["protein"].casefold() != CHAINS[chain].casefold():
            raise ValueError(f"{path}: unexpected TM chain/protein {chain}/{row['protein']}")
        t, value = float(row["time_ns"]), float(row[metric])
        if not math.isfinite(t) or not math.isfinite(value) or value < 0:
            raise ValueError(f"{path}: invalid time or RMSD for chain {chain}")
        series[chain].append((t, value))
    expected = np.arange(1001) / 10
    for chain, points in series.items():
        points.sort()
        if len(points) != 1001 or not np.allclose([p[0] for p in points], expected,
                                                 rtol=0, atol=1e-8):
            raise ValueError(f"{path}: chain {chain} needs 1001 unique times from 0 to 100 ns at 0.1 ns")
    return series


def load_contacts(path):
    matrix = np.full((len(PAIRS), 100), np.nan)
    for row in read_tsv(path, ("pair", "start_ns", "end_ns", "n_frames", "value")):
        pair = row["pair"]
        start, end, frames = (int(row[key]) for key in ("start_ns", "end_ns", "n_frames"))
        value = float(row["value"])
        if (pair not in PAIRS or not 0 <= start < 100 or end != start + 1
                or frames != (11 if end == 100 else 10) or not 0 <= value <= 1):
            raise ValueError(f"{path}: invalid bin for {pair} at {start} ns")
        i = PAIRS.index(pair)
        if not np.isnan(matrix[i, start]):
            raise ValueError(f"{path}: duplicate bin for {pair} at {start} ns")
        matrix[i, start] = value
    if np.isnan(matrix).any():
        raise ValueError(f"{path}: expected 100 populated 1 ns bins per pair")
    return matrix


def load_pore(path):
    profiles = {time: [] for time in SNAPSHOTS.values()}
    for row in read_tsv(path, ("frame", "time_ns", "channel_coordinate_A", "radius_A")):
        frame, t = int(row["frame"]), float(row["time_ns"])
        z, radius = float(row["channel_coordinate_A"]), float(row["radius_A"])
        if (frame not in SNAPSHOTS or t != SNAPSHOTS[frame] or not CORE[0] <= z <= CORE[1]
                or not math.isfinite(radius) or radius < 0):
            raise ValueError(f"{path}: invalid pore snapshot point for frame {frame}")
        profiles[t].append((z, radius))
    expected = np.arange(-133, 136) / 8
    for t, points in profiles.items():
        points.sort()
        if len(points) != len(expected) or not np.array_equal([p[0] for p in points], expected):
            raise ValueError(f"{path}: expected 269 unique core coordinates at {t} ns")
    return profiles


def style_axes(ax):
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(length=2.3, pad=2, width=.5)


def draw_contacts(rmsd, contacts):
    fig = plt.figure(figsize=(2.85, 2.25))
    fig.text(.015, .985, "F", fontsize=10, va="top")
    upper = fig.add_axes([.205, .51, .665, .405])
    lower = fig.add_axes([.205, .17, .665, .205])
    cax = fig.add_axes([.89, .17, .022, .205])

    for chain, points in rmsd.items():
        t, values = zip(*points)
        upper.plot(t, values, lw=.55, color=PURPLE if chain in "AB" else CYAN,
                   linestyle=(0, (3.7, 1.6)) if chain in "BD" else "-")
    max_rmsd = max(value for points in rmsd.values() for _, value in points)
    limit = max(.2, math.ceil(max_rmsd / .05) * .05)
    upper.set(xlim=(0, 100), ylim=(0, limit), xticks=(0, 50, 100), ylabel="TM RMSD (nm)")
    upper.set_yticks(np.linspace(0, limit, 5), [f"{v:.2f}" for v in np.linspace(0, limit, 5)])
    upper.tick_params(labelbottom=False)
    upper.text(.30, 1.025, "IR75g A/B", color=PURPLE, transform=upper.transAxes, ha="center")
    upper.text(.73, 1.025, "IR8a C/D", color=CYAN, transform=upper.transAxes, ha="center")
    style_axes(upper)

    cmap = LinearSegmentedColormap.from_list(
        "retention", ("#f6f5fb", "#d3d3eb", "#9f9fca", "#68689e", "#3e426e"))
    heatmap = lower.pcolormesh(np.arange(101), np.arange(5), contacts,
                               cmap=cmap, vmin=0, vmax=1, shading="flat")
    lower.set(xlim=(0, 100), ylim=(4, 0), xticks=(0, 50, 100), xlabel="Time (ns)")
    lower.set_yticks(np.arange(4) + .5, PAIRS)
    lower.set_title("Initial-contact retention", loc="left", pad=2)
    style_axes(lower)
    lower.tick_params(axis="y", length=0)
    bar = fig.colorbar(heatmap, cax=cax, ticks=(0, .5, 1), format="%.1f")
    bar.ax.tick_params(length=2, pad=2, width=.5)
    bar.solids.set_rasterized(False)
    return fig


def draw_pore(profiles):
    fig, axes = plt.subplots(1, 3, figsize=(3.4, 2.06), sharey=True)
    fig.subplots_adjust(left=.16, right=.983, bottom=.19, top=.88, wspace=.22)
    fig.text(.015, .985, "G", fontsize=10, va="top")
    max_radius = max(radius for points in profiles.values() for _, radius in points)
    for ax, (t, points), color in zip(axes, profiles.items(), (DARK, CYAN, PURPLE)):
        z, radius = zip(*points)
        ax.plot(radius, z, color=color, lw=.6)
        ax.set(xlim=(0, max(3.5, max_radius)), ylim=CORE,
               yticks=(-15, 0, 15), xlabel="Radius (\u00c5)")
        ax.set_xticks((0, 1.5, 3), ("0.0", "1.5", "3.0"))
        ax.set_title(f"{t} ns", pad=3)
        style_axes(ax)
    axes[0].set_ylabel("Pore coordinate (\u00c5)")
    return fig


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--panel", choices=("F", "G", "all"), default="all")
    parser.add_argument("--domain-rmsd", type=Path,
                        help="Use a retained domain_rmsd.tsv in place of inputs/md/domain_rmsd.tsv.")
    args = parser.parse_args()
    inputs = Path(os.environ.get("HONEYBEE_INPUT_DIR", "inputs"))
    out = Path(os.environ.get("HONEYBEE_FIGURE_DIR", "figures"))
    plt.rcParams.update({"font.family": ["Arial", "DejaVu Sans"], "font.size": 6,
                         "axes.titlesize": 6, "axes.linewidth": .55, "pdf.fonttype": 42,
                         "svg.fonttype": "none", "path.simplify": False})
    f = g = None
    try:
        if args.panel in ("F", "all"):
            path = args.domain_rmsd or inputs / "md/domain_rmsd.tsv"
            if not path.is_file():
                raise ValueError(f"missing {path}; provide the retained domain_rmsd.tsv via "
                                 "HONEYBEE_INPUT_DIR or --domain-rmsd (see INPUTS.md)")
            f = draw_contacts(load_rmsd(path), load_contacts(DATA / "interface_heatmap_1ns_v2.tsv"))
        if args.panel in ("G", "all"):
            g = draw_pore(load_pore(DATA / "pore_radius_snapshots.tsv"))
    except (OSError, ValueError, KeyError) as exc:
        raise SystemExit(f"Figure S5F-G: {exc}") from None

    out.mkdir(parents=True, exist_ok=True)
    for ext in FORMATS:
        if f is not None:
            # Figure S5F
            f.savefig(out / f"FigureS5F_TM_RMSD_and_contacts.{ext}", dpi=300)
        if g is not None:
            # Figure S5G
            g.savefig(out / f"FigureS5G_pore_radius_profiles.{ext}", dpi=300)
    if f is not None:
        plt.close(f)
    if g is not None:
        plt.close(g)


if __name__ == "__main__":
    main()
