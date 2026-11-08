"""Draw the beam alignment and ligand annotation counts for Figure S5."""

import csv
import os
import re
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
from matplotlib import pyplot as plt
from matplotlib.patches import Rectangle


DATA = Path("data/supplementary")
FORMATS = ("pdf", "svg", "png")
LAYERS = (
    ("honeybee_semiochemical", "Honeybee semiochemicals", "#2e7d32"),
    ("floral_volatile", "Floral/plant volatiles", "#8ebbe6"),
    ("ir_acid_amine", "IR-biased acids/amines", "#ffe297"),
    ("reactive_odorant", "Reactive odorants", "#e99195"),
)
GROUPS = ("carboxylic acid", "aldehyde", "ester", "amine", "ketone", "thiol")
BEAM_COLORS = {"Amel_IR25a": "#007f5f", "Amel_IR8a": "#2691aa"}


def read_tsv(path):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def load_beam():
    path = DATA / "beam_alignment.tsv"
    rows = read_tsv(path)
    by_id = {r["sequence_id"]: r for r in rows}
    visible = sorted((r for r in rows if int(r["display_order"]) > 0),
                     key=lambda r: int(r["display_order"]))
    orders = [int(r["display_order"]) for r in visible]
    if len(rows) != 19 or len(by_id) != 19 or orders != list(range(1, 15)):
        raise ValueError(f"{path}: expected 19 unique scaffold rows and display orders 1-14")
    for row in rows:
        seq = row["aligned_sequence"]
        if (not re.fullmatch(r"[A-Z-]{80}", seq)
                or len(seq.replace("-", "")) != int(row["end"]) - int(row["start"]) + 1):
            raise ValueError(f"{path}: invalid alignment or residue range for {row['sequence_id']}")

    path = DATA / "beam_features.tsv"
    highlights = set()
    helices = defaultdict(list)
    for row in read_tsv(path):
        key, feature = row["sequence"], row["feature"]
        if key not in by_id:
            raise ValueError(f"{path}: unknown sequence {key}")
        if feature == "published_red_box_homologous_region":
            highlights.add(tuple(int(row[f"alignment_column_{end}_1based"]) - 1
                                 for end in ("start", "end")))
        elif feature.startswith("model_DSSP_"):
            helices[key].append((feature.removeprefix("model_DSSP_"),
                                  int(row["start"]), int(row["end"])))
        else:
            raise ValueError(f"{path}: unknown feature {feature}")
    if (len(highlights) != 1 or set(helices) != set(BEAM_COLORS)
            or any([h[0] for h in helices[key]] != [f"H{i}" for i in range(1, n + 1)]
                   for key, n in (("Amel_IR25a", 4), ("Amel_IR8a", 3)))):
        raise ValueError(f"{path}: expected a shared beam interval and both Amel helix tracks")
    highlight = highlights.pop()
    if not 0 <= highlight[0] <= highlight[1] < 80:
        raise ValueError(f"{path}: beam interval is outside the alignment")
    return rows, visible, highlight, helices


def load_counts():
    path = DATA / "ligand_annotation_counts.tsv"
    rows = read_tsv(path)
    by_group = {r["functional_group"]: r for r in rows}
    if len(rows) != len(GROUPS) or set(by_group) != set(GROUPS):
        raise ValueError(f"{path}: expected one row per displayed functional group")
    counts = [[int(by_group[g][layer]) for layer, _, _ in LAYERS] for g in GROUPS]
    if any(n < 0 for row in counts for n in row):
        raise ValueError(f"{path}: annotation counts must be nonnegative")
    return counts


def draw_beam(rows, visible, highlight, helices):
    shades = []
    for col in zip(*(r["aligned_sequence"] for r in rows)):
        residues = Counter(aa for aa in col if aa != "-")
        frequency = max(residues.values()) / residues.total() if residues else 0
        shades.append("#d8d8df" if frequency >= .8 else "#ededf1" if frequency >= .5 else None)

    fig, ax = plt.subplots(figsize=(5.35, 1.95))
    fig.subplots_adjust(left=.02, right=.98, bottom=.025, top=.98)
    ax.set(xlim=(-14, 92), ylim=(17.5, -2.3))
    ax.set_axis_off()
    ax.text(-14, -2.15, "A", va="top", fontsize=10)
    lo, hi = highlight
    ax.plot([lo - .35, hi + .35], [-.95, -.95], color="#c33b2b", lw=.6)
    ax.text((lo + hi) / 2, -1.15, "beam-like region", ha="center",
            va="bottom", color="#c33b2b")
    y = .1
    for row in visible:
        key = row["sequence_id"]
        color = BEAM_COLORS.get(key, "#222222")
        ax.text(-13.5, y, key.replace("_", " "), va="center", color=color)
        ax.text(82, y, f"{row['start']}-{row['end']}", va="center")
        residue = int(row["start"])
        positions = {}
        for x, aa in enumerate(row["aligned_sequence"]):
            if aa != "-":
                positions[residue] = x
                residue += 1
                if shades[x]:
                    ax.add_patch(Rectangle((x - .49, y - .39), .98, .78,
                                           facecolor=shades[x], edgecolor="none"))
            ax.text(x, y, aa, ha="center", va="center", color=color,
                    fontfamily="DejaVu Sans Mono", fontsize=5.65)
        for label, start, end in helices.get(key, ()):
            occupied = [x for pos, x in positions.items() if start <= pos <= end]
            if not occupied:
                raise ValueError(f"No aligned residues for {key} {label}")
            a, b = occupied[0], occupied[-1]
            ax.plot([a - .35, b + .35], [y + .58] * 2, lw=.6, color=color)
            ax.text((a + b) / 2, y + .68, label, ha="center", va="top", color=color, fontsize=5.4)
        y += 2.2 if key in helices else 1.05
    return fig


def draw_ligands(counts):
    fig, ax = plt.subplots(figsize=(2.5, 1.88))
    fig.subplots_adjust(left=.23, right=.985, top=.97, bottom=.185)
    fig.text(.018, .98, "B", va="top", fontsize=10)
    y = list(range(len(GROUPS)))
    left = [0] * len(GROUPS)
    for column, (_, label, color) in zip(zip(*counts), LAYERS):
        ax.barh(y, column, left=left, height=.72, color=color,
                edgecolor="#555555", linewidth=.35, label=label)
        left = [a + b for a, b in zip(left, column)]
    ax.set(xlim=(0, max(21.2, max(left) + .2)), ylim=(5.65, -.65),
           xlabel="Number of annotated ligands")
    ax.set_yticks(y, ("carboxylic\nacid", *GROUPS[1:]))
    ax.set_xticks(range(0, 21, 5), [f"{x:.1f}" for x in range(0, 21, 5)])
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(axis="both", length=2, pad=2)
    ax.set_axisbelow(True)
    ax.grid(axis="x", color="#dcdcdc", linewidth=.4)
    ax.legend(loc="lower right", fontsize=5.5, handlelength=1.6,
              borderpad=.4, labelspacing=.25, framealpha=1, edgecolor="#dddddd")
    return fig


def main():
    plt.rcParams.update({"font.family": ["Arial", "DejaVu Sans"], "font.size": 6,
                         "axes.linewidth": .55, "pdf.fonttype": 42, "svg.fonttype": "none"})
    try:
        beam = draw_beam(*load_beam())
        ligands = draw_ligands(load_counts())
    except (OSError, ValueError, KeyError) as exc:
        raise SystemExit(f"Figure S5A-B: {exc}") from None

    out = Path(os.environ.get("HONEYBEE_FIGURE_DIR", "figures"))
    out.mkdir(parents=True, exist_ok=True)
    for ext in FORMATS:
        # Figure S5A
        beam.savefig(out / f"FigureS5A_beam_sequence.{ext}", dpi=300)
        # Figure S5B
        ligands.savefig(out / f"FigureS5B_ligand_annotation_counts.{ext}", dpi=300)
    plt.close(beam)
    plt.close(ligands)


if __name__ == "__main__":
    main()
