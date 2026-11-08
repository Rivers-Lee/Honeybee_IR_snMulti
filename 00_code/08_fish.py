"""Draw the Figure S2 RNA-FISH summaries from section counts."""

import csv
import math
import os
from pathlib import Path
from statistics import mean

import matplotlib

matplotlib.use("Agg")
from matplotlib import pyplot as plt
from matplotlib.patches import Circle
from matplotlib.ticker import PercentFormatter


SOURCE = Path("data/fish/replicate_summary.csv")
FIELDS = ("Rep1", "Rep2", "Rep3")
ORCO_COLOR = "#6a4c93"
IR76B_COLOR = "#8b1e3f"
FORMATS = ("pdf", "svg", "png")


def load_counts(path=SOURCE):
    with path.open(newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle))
    by_field = {row["replicate"]: row for row in rows}
    if len(rows) != len(FIELDS) or set(by_field) != set(FIELDS):
        raise ValueError(f"{path}: expected one row per section (Rep1-Rep3)")
    counts = []
    for field in FIELDS:
        row = by_field[field]
        total, orco, ir76b, double = (
            int(row[key]) for key in ("total_cell_territories", "orco_positive_cells",
                                      "ir76b_positive_cells", "double_positive_cells")
        )
        if not (0 < orco <= total and 0 <= ir76b <= total
                and max(0, orco + ir76b - total) <= double <= min(orco, ir76b)):
            raise ValueError(f"{path}: inconsistent cell-territory counts for {field}")
        counts.append((orco, ir76b, double))
    return counts


def overlap_area(a, b, d):
    if d >= a + b:
        return 0.0
    if d <= abs(a - b):
        return math.pi * min(a, b) ** 2
    x = (d * d + a * a - b * b) / (2 * d)
    h = math.sqrt(max(0.0, a * a - x * x))
    return a * a * math.acos(x / a) + b * b * math.acos((d - x) / b) - d * h


def euler_geometry(orco, ir76b, double):
    a, b = (math.sqrt(n / math.pi) for n in (ir76b, orco))
    if double == 0:
        return a, b, a + b
    if double == min(orco, ir76b):
        return a, b, abs(a - b)
    lo, hi = abs(a - b), a + b
    for _ in range(60):
        mid = (lo + hi) / 2
        if overlap_area(a, b, mid) > double:
            lo = mid
        else:
            hi = mid
    return a, b, (lo + hi) / 2


def draw_euler(counts):
    orco, ir76b, double = map(mean, zip(*counts))
    a, b, distance = euler_geometry(orco, ir76b, double)
    left, right = min(-a, distance - b), max(a, distance + b)
    scale = min(math.sqrt(math.pi) * 5.175, 59.6 / max(a, b), 150 / (right - left))
    x_ir = 96 - scale * (left + right) / 2
    x_orco = x_ir + scale * distance
    y = 67.5

    fig = plt.figure(figsize=(236 / 72, 131 / 72))
    ax = fig.add_axes((0, 0, 1, 1), xlim=(0, 236), ylim=(0, 131), aspect="equal")
    ax.set_axis_off()
    ax.text(4, 126, "B", fontsize=10, va="top")
    for x, radius, color, alpha in (
        (x_ir, a, IR76B_COLOR, .28), (x_orco, b, ORCO_COLOR, .18)
    ):
        ax.add_patch(Circle((x, y), radius * scale, color=color, alpha=alpha, linewidth=0))
    for x, radius, color in ((x_ir, a, IR76B_COLOR), (x_orco, b, ORCO_COLOR)):
        ax.add_patch(Circle((x, y), radius * scale, fill=False, edgecolor=color, linewidth=1))
    ax.text(25.5, 118.5, "Ir76b+", color=IR76B_COLOR, style="italic", fontsize=8.5)
    ax.text(148.7, 118.5, "Orco+", color=ORCO_COLOR, style="italic", fontsize=8.5)
    if double > 0:
        common_left = max(x_ir - a * scale, x_orco - b * scale)
        common_right = min(x_ir + a * scale, x_orco + b * scale)
        ax.text((common_left + common_right) / 2, y, "Double+",
                ha="center", va="center", fontsize=8.5)
    ax.text(182.6, 113, "Mean cells")
    for i, (name, value, color) in enumerate((
        ("Orco+", orco, ORCO_COLOR), ("Ir76b+", ir76b, IR76B_COLOR),
        ("Double+", double, "#111111")
    )):
        ax.text(182.6, 99.6 - i * 10.35, f"{name}: {value:.1f}", color=color)
    for x, value, label in ((189.3, 1, "1 cell"), (212.3, 5, "5 cells")):
        radius = scale * math.sqrt(value / math.pi)
        ax.add_patch(Circle((x, 48.5), radius, fill=False, edgecolor="#111111", linewidth=1))
        ax.text(x, 29, label, ha="center")
    return fig


def draw_bars(counts):
    fig = plt.figure(figsize=(207 / 72, 131 / 72))
    ax = fig.add_axes((50 / 207, 16 / 131, 150.5 / 207, 106.7 / 131))
    fig.text(3 / 207, 1 - 5 / 131, "C", fontsize=10, va="top")
    values = [100 * double / orco for orco, _, double in counts]
    ax.set(xlim=(-.5, 2.5), ylim=(0, 100), yticks=range(0, 101, 25),
           ylabel=r"% $\mathit{Ir76b+}$ among $\mathit{Orco+}$ cells")
    ax.set_xticks(range(3), [f"Field{i}" for i in range(1, 4)], fontsize=8)
    ax.yaxis.set_major_formatter(PercentFormatter(100, decimals=0))
    ax.set_axisbelow(True)
    ax.grid(axis="y", color="#dddddd", linewidth=.75)
    ax.tick_params(axis="both", length=0, pad=2.4, labelsize=7)
    ax.tick_params(axis="x", labelsize=8)
    ax.bar(range(3), values, width=.4255, color=IR76B_COLOR, edgecolor="#111111", linewidth=.8)
    for i, value in enumerate(values):
        ax.annotate(f"{value:.1f}%", (i, value), xytext=(0, 1.4),
                    textcoords="offset points", ha="center", va="bottom")
    return fig


def main():
    counts = load_counts()
    plt.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"],
        "font.size": 7, "axes.spines.top": False, "axes.spines.right": False,
        "axes.linewidth": .75, "pdf.fonttype": 42, "svg.fonttype": "none",
        "mathtext.fontset": "custom", "mathtext.rm": "sans", "mathtext.it": "sans:italic",
    })
    out = Path(os.environ.get("HONEYBEE_FIGURE_DIR") or "figures")
    out.mkdir(parents=True, exist_ok=True)

    fig = draw_euler(counts)
    for ext in FORMATS:
        # Figure S2B
        fig.savefig(out / f"FigureS2B_Euler_mean_counts.{ext}", dpi=300)
    plt.close(fig)

    fig = draw_bars(counts)
    for ext in FORMATS:
        # Figure S2C
        fig.savefig(out / f"FigureS2C_FISH_conditional_bars.{ext}", dpi=300)
    plt.close(fig)


if __name__ == "__main__":
    main()
