# Data dictionary

`data/figure_source_data/README_source_data.md` provides a file-by-file inventory. The most frequently used files are summarized below.

| File or group | Contents |
| --- | --- |
| `main139_Figure1B_*` | Pairwise ANI matrix and type-strain ANI/coverage values for the 139-genome cohort. |
| `main139_144tip_rooted.nwk`, `main139_139tip_midpoint.nwk` | Rooted outgroup-inclusive and midpoint-rooted ingroup phylogenetic trees. |
| `main139_iTOL_*` | Plain-text iTOL labels, color strips, and legends for Figure 2. |
| `main139_postprocess.*`, `main139_pav_matrix.tsv.gz` | PGAP2 partition statistics, rarefaction results, fitted curves, and the compressed PAV matrix. |
| `main139_pcoa_coordinates.csv`, `main139_pam_cluster_assignments.csv` | Accessory-genome PCoA coordinates and PAM assignments. |
| `main139_pam_associated_genes.csv`, `main139_cog_*` | PAM-associated clusters and COG composition/enrichment results. |
| `main139_amr_*`, `main139_tetm_*` | Per-sample AMR calls and tet(M) summaries. |
| `main139_cluster_mge_status.tsv`, `main139_mge_regions_per_sample.tsv`, `main139_genomic_islands_per_sample.tsv` | Predicted MGE and genomic-island summaries. |
| `main139_defense_*` | Per-sample and group-level DefenseFinder summaries. |
| `data/script_inputs/` | Direct, low-volume inputs retained for the released Figure 3, 5, and 6 scripts. |
| `metadata/Supplementary_Table_S1...S6.xlsx` | Submission-ready supplementary tables. |

All tabular text files are UTF-8, tab- or comma-delimited. The PAV matrix is gzip-compressed and can be read directly in R with `readr::read_tsv()`.
