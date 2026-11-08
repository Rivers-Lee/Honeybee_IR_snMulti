## RNA and Accessibility

Download the honeybee multiome FASTQs (GEO accession GSE248097). Follow
Zhang et al. (2025) and the manuscript's *Generation of gene expression and
chromatin accessibility matrices* and *OSN subtype classification* to
prepare the final OSN object. Keep the RNA normalized-data layer, final
cluster identities (0-59) and RNA t-SNE coordinates for the same 4,175
nuclei.
```bash
Rscript 00_osn_embedding.R export /path/to/final_OSN.rds
Rscript 01_cluster_expression.R
```

The first command exports `inputs/osn/osn_cells.tsv` and draws the RNA maps.
It reads `Annotation`, the existing `tsne.rna` reduction, and the RNA
`data` and `counts` layers. If you already have the TSV, omit the
`export` arguments. Keep unique `cell` IDs, `cluster`, `rnaTSNE_1`,
`rnaTSNE_2`, and normalized `Orco`, `Ir76b`, `Ir75g`, `Ir25a`, `Ir8a`,
`Ir75f`, `Ir76a`, `Ir75e` columns. Add corresponding `_UMI` columns for
the violins. The export requires SeuratObject; the maps and violins use
ggplot2.

`data/osn/figure4b_calls.tsv` contains one eight-symbol string per
cluster, in final display order and the gene order above. `L` and `S`
are the two RNA categories, `A` is accessibility-only, and `.` is below
the criteria. These are the published calls, not per-nucleus ATAC data.
To regenerate them, use the cluster RNA and promoter accessibility rules
in STAR Methods. Run `Rscript 01_cluster_expression.R evidence` to draw
Figure 4B and export the class tables for Figure 4C without the RNA TSV.

## Accessibility Tracks

Use the final OSN clusters and indexed fragments from the Amel_HAv3.1
Signac workflow. For each displayed locus, retain the normalized coverage
values from the plotted Signac tracks, after the smoothing and scaling
used for that panel. Do not substitute promoter-accessibility calls for
the continuous coverage.

- `inputs/osn/locus_coverage.tsv`: `gene`, `group`, `chr`,
  `position_bp`, `signal`, `y_min`, `y_max`. Use one row per genomic
  sample and group. Retain the plotted coordinate window and the lower
  and upper displayed signal limits for each group. The script clips
  display heights at these limits without changing the values.
- `inputs/osn/locus_exons.tsv`: `gene`, `start_bp`, `end_bp`, `strand`.
  Use the exon intervals overlapping the same window, in 1-based
  coordinates, and `+` or `-` for the transcript orientation.
- Figure 2A: supply Orco (Group1), Ir76b (Group14) and Ir75g (Group2),
  each with groups `0` through `59`. Run `Rscript 02_locus_tracks.R fig2`.
- Figure S4B-F: supply Ir25a, Ir8a, Ir76a, Ir75e and Ir75f with groups
  `positive` and `negative`. In the retained Ir25a/Ir8a grouping,
  clusters with RNA or accessibility support form the positive group.
  For the specific genes, positive
  clusters are 0 and 5 for Ir76a, 0 for Ir75e and 5 for Ir75f. The
  accompanying RNA panels use `inputs/osn/osn_cells.tsv` from script 00.
  Run `Rscript 02_locus_tracks.R s4`.

The original export settings are in the coverage sections of `02_recall_peak_cov_plot_final.R`
and `Fig3_S3.R`; use the window retained for the final panel.

## Repertoire and Chromosomes

- Figure S1A: `inputs/repertoire/receptor_counts.tsv` must contain
  `species`, `family`, `total_annotated`, and `olfactory_detected`. Use
  one row for each combination of `A.mellifera`/`D.melanogaster` and
  `IR`/`OR`. Use the final repertoire and expression counts, not gene
  density or the number of tips in the tree. Run
  `Rscript 03_repertoire_genome.R counts`.
- Figure S1C: use Amel_HAv3.1 chromosome lengths and the curated
  annotation. Count GTF `transcript` start coordinates in the 200-kb
  windows used by `06_RIdeogram.R`. Supply
  `inputs/repertoire/transcript_density.tsv` with integer `chr` (1-16),
  `start`, `end`, `count`; starts are 0, 200000, 400000, and so on,
  and the final end is clipped to chromosome length. Run
  `Rscript 03_repertoire_genome.R ideogram`. The 11 study loci are in
  `data/repertoire/ir_loci.tsv`. To redraw the displayed bands without
  the numeric density table, run `Rscript 03_repertoire_genome.R ideogram-display`;
  the supplied band colors are display values, not gene counts.

## Bundled Figure Tables

| Script | Files under `data/` | How to obtain the intermediate |
| --- | --- | --- |
| 04 | `expression/Figure4D_RNA_only_plot_values.csv` | STAR Methods, RNA-only comparison. One row per species and receptor combination, with `species`, `combination`, `clusters`, `cells`, `cluster_pct`. Percentages use denominators 60, 47 and 43, respectively. Do not use the joint RNA/accessibility calls. |
| 05 | `sequence/Figure3_sequences.fasta`, `sequence/Figure3_signature_hits_with_homology.tsv` | Supplementary Data S1, sequence features. Use canonical AmelIR75g, AsteIR75k and RproIR75b. The TSV contains inclusive 1-based Pfam or homology-mapped intervals. The script computes identity from global BLOSUM62 alignments, excluding gap pairs. |
| 06 | `structures/Figure3E_rmsd_A.tsv`, `structures/Figure3E_tm_score.tsv`, `structures/Figure5B_tm_score.tsv` | Supplementary Data S1, structural alignments. One square, labeled TSV per metric; use the symmetric scores obtained by rescoring the retained USalign pairs with the pair-mean observed C-alpha count. Do not average the two directional scores. |
| 07 | `docking/composition_summary.tsv`, `docking/N30_display_order_matrix.tsv`, `docking/ligand_names.tsv` | Supplementary Data S1, composition metrics and ligand scores. Use the five 25-model means (no baseline) and the final ordered 5 x 30 score matrix from the 48-ligand panel. The matrix contains size-adjusted scores, not raw affinities. |
| 08 | `fish/replicate_summary.csv` | Supplementary Data S1, three sections. Use the CellProfiler territory summary; each section is one field. |
| 09 | `contacts/*_contact_matrix.tsv`, `contacts/C1_C2_contact_plot_values.tsv` | Supplementary Data S1, residue contacts. AF3 matrix entries are medians across 25 predictions; AF3 bars are means and SDs of five seed means. The 43EF value is from one experimental structure. |
| 10 | `supplementary/beam_alignment.tsv`, `supplementary/beam_features.tsv`, `supplementary/ligand_annotation_counts.tsv` | Retained beam alignment and annotation coordinates; final functional-group counts by curation layer. The counts are the plotted values, not compound-level records. Compounds may contribute to more than one group. |
| 11 | `md/interface_heatmap_1ns_v2.tsv`, `md/pore_radius_snapshots.tsv` | Retained one-nanosecond contact-retention bins and 0, 50 and 100 ns HOLE profiles. The only external input is the retained TMD RMSD table for Figure S5F. |

## Figure S5F TMD RMSD

Follow the MD trajectory analysis in STAR Methods for the ligand-free
IR75g2-IR8a2 model. After a collective TMD-core fit, export the retained
`04_domains/domain_rmsd.tsv` to `inputs/md/domain_rmsd.tsv`.
Required columns: `time_ns`, `chain`, `protein`, `domain`,
`TM_assembly_fit_rmsd_nm`. The script selects `TM_helices` and expects
chains A/B (IR75g), C/D (IR8a) with 1,001 time points each, 0 to 100 ns at
0.1-ns spacing. The contact heatmap is already supplied.

```bash
python 11_s5_md.py --domain-rmsd /path/to/domain_rmsd.tsv
python 11_s5_md.py --panel G
```