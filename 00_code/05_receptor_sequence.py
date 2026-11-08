import hashlib
import itertools
import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import Align, SeqIO
from Bio.Align import substitution_matrices
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle


DATA = Path("data/sequence")
OUT = Path(os.environ.get("HONEYBEE_FIGURE_DIR", "figures"))
IDS = ("Amel_Ir75g", "Aste_Ir75k", "Rpro_Ir75b")
LABELS = ("Amel75g", "Aste75k", "Rpro75b")


def percent_identity(sequences):
    aligner = Align.PairwiseAligner(mode="global")
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -.5
    identity = np.eye(len(IDS)) * 100
    for i, j in itertools.combinations(range(len(IDS)), 2):
        a, b = aligner.align(sequences[IDS[i]], sequences[IDS[j]])[0]
        paired = [(x, y) for x, y in zip(a, b) if x != "-" and y != "-"]
        identity[i, j] = identity[j, i] = 100 * sum(x == y for x, y in paired) / len(paired)
    return identity


def plot_identity(values):
    fig, ax = plt.subplots(figsize=(2.2, 1.82), layout="constrained")
    palette = plt.get_cmap("Purples")
    for y, x in np.ndindex(values.shape):
        value = values[y, x]
        color = "#e6e6e6" if x == y else palette(np.clip((value-15)/20, 0, 1))
        ax.add_patch(Rectangle((x-.5, y-.5), 1, 1, facecolor=color, edgecolor="white", lw=.5))
        ax.text(x, y, f"{value:.1f}", ha="center", va="center",
                color="white" if x != y and value > 29 else "#222222")
    ax.set(xlim=(-.5, 2.5), ylim=(2.5, -.5), aspect="equal", title="Sequence identity (%)")
    ax.set_xticks(range(3), LABELS)
    ax.set_yticks(range(3), LABELS)
    ax.tick_params(length=0, pad=3)
    ax.spines[:].set_visible(False)
    # Figure 3A
    fig.savefig(OUT / "Figure3A_sequence_identity.pdf")
    plt.close(fig)


def plot_signatures(sequences, hits):
    palette = {"PF10613": "#f06a10", "PF10613_like_homology_mapped": "#f06a10",
               "PF00060": "#cf370c"}
    if set(hits["sequence_id"]) != set(IDS) or not set(hits["signature"]) <= set(palette):
        raise ValueError("Unexpected signature annotations.")
    fig = plt.figure(figsize=(2.35, 1.62))
    ax = fig.add_axes([.255, .30, .70, .62])
    ax.set(xlim=(0, 725), ylim=(-.4, 2.45))
    for y, protein in zip([2, 1, 0], IDS):
        seq = sequences[protein]
        length = len(seq)
        ax.plot([0, length], [y, y], color="#999999", lw=.5, zorder=1)
        ax.text(length+12, y, str(length), va="center")
        for row in hits[hits["sequence_id"].eq(protein)].itertuples():
            start, end = row.start_aa_1based, row.end_aa_1based
            if row.sequence_sha256 != hashlib.sha256(seq.encode("ascii")).hexdigest():
                raise ValueError(f"Signature sequence differs from FASTA: {protein}")
            if not 1 <= start <= end <= length:
                raise ValueError(f"Invalid signature interval: {protein}, {start}-{end}")
            ax.add_patch(Rectangle((start-.5, y-.13), end-start+1, .26,
                                   facecolor=palette[row.signature], edgecolor="none"))
    ax.set_yticks([2, 1, 0], LABELS)
    ax.set_xticks([0, 200, 400, 600])
    ax.set_xlabel("Residue position", labelpad=1)
    ax.tick_params(axis="y", length=0, pad=4)
    ax.tick_params(axis="x", length=2, width=.5, pad=2)
    ax.spines[["left", "top", "right"]].set_visible(False)
    handles = [Line2D([0], [0], color=palette[key], lw=3, label=label)
               for key, label in (("PF10613", "Lig_chan-Glu_bd"), ("PF00060", "Lig_chan"))]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(.57, .0),
               ncols=2, frameon=False, fontsize=6, handlelength=1.3, columnspacing=.8)
    # Figure 3B
    fig.savefig(OUT / "Figure3B_sequence_signatures.pdf")
    plt.close(fig)


def main():
    plt.rcParams.update({"font.family": ["Arial", "DejaVu Sans"], "font.size": 6.5,
                         "pdf.fonttype": 42, "axes.linewidth": .5})
    records = SeqIO.to_dict(SeqIO.parse(DATA / "Figure3_sequences.fasta", "fasta"))
    sequences = {key: str(records[key].seq) for key in IDS}
    hits = pd.read_csv(DATA / "Figure3_signature_hits_with_homology.tsv", sep="\t")
    OUT.mkdir(parents=True, exist_ok=True)
    plot_identity(percent_identity(sequences))
    plot_signatures(sequences, hits)


if __name__ == "__main__":
    main()
