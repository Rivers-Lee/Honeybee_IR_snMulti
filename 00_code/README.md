# Honeybee IR Figure Redraws

Run the scripts from this directory with Python 3.10+ or R as indicated.
Start with the numbered file for the panel you want. The scripts redraw data
panels, not the final Illustrator page layouts.

```bash
python -m pip install -r requirements.txt
python 05_receptor_sequence.py
Rscript 01_cluster_expression.R evidence
```

`data/` contains small tables and sequences supplied with this package.
`inputs/` is for intermediates prepared by the reader; see
[Inputs](INPUTS.md). `figures/` is the default output directory.
Set `HONEYBEE_INPUT_DIR` and `HONEYBEE_FIGURE_DIR` to use other directories.
Python scripts use Arial when available and DejaVu Sans otherwise.

## Find a Panel

| Panel | Script | Input |
| --- | --- | --- |
| Figure 1C-E, Figure 4A | `00_osn_embedding.R` | Final OSN RNA/t-SNE data |
| Figure 1F, Figure 4B, Figure S4A | `01_cluster_expression.R` | Cluster calls and final per-nucleus RNA; `evidence` needs only the calls |
| Figure 4C (supporting counts) | `01_cluster_expression.R evidence` | Cluster assignments and complement counts; the schematic is separate |
| Figure 2A, Figure S4B-F | `02_locus_tracks.R` | Prepared normalized accessibility, exon coordinates and per-nucleus RNA for S4D-F |
| Figure S1A (bars), Figure S1C | `03_repertoire_genome.R` | Repertoire counts, 200-kb density and IR loci |
| Figure 4D | `04_receptor_combinations.py` | RNA-only comparison table |
| Figure 3A-B | `05_receptor_sequence.py` | Protein FASTA and signature intervals |
| Figure 3E, Figure 5B | `06_structural_similarity.py` | C-alpha RMSD and symmetric TM-score matrices |
| Figure 5A, Figure 5E | `07_compositions_and_ligands.py` | AF3/PRODIGY summary and displayed 30-ligand score matrix |
| Figure S2B-C | `08_fish.py` | Three-section cell territory summary |
| Figure S3A-B | `09_contacts.py` | Residue-contact matrices and co-receptor summary |
| Figure S5A-B | `10_s5_beam_ligand.py` | Alignment, helix annotations and ligand group counts |
| Figure S5F-G (plots) | `11_s5_md.py` | TMD RMSD, contact retention and HOLE profiles |

For panels made in ChimeraX, MEGA11, PoseEdit, ImageJ or Illustrator, see
[Software Panels](SOFTWARE_PANELS.md). No plotting stubs are included.
