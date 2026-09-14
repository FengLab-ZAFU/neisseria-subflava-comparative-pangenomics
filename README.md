# Comparative pangenomics of the Neisseria subflava-associated group

This repository contains the code and supporting data for our study, *Comparative pangenomics of the Neisseria subflava-associated group reveals structured accessory diversity and lineage-associated tet(M) carriage*.

The study cohort comprises 139 genomes: 47 *N. subflava*, 82 *N. perflava*, and 10 *N. flavescens*. Genospecies assignments are based on reference ANI comparisons.

## Files

- `code/`: R scripts used for the analyses and figures.
- `code/figure_scripts/`: scripts for Figures 1B and 3-6 and Supplementary Figure S2.
- `data/figure_source_data/`: source data for the manuscript figures, including the ANI matrix, trees, pangenome summaries, PCoA coordinates, PAM assignments, and AMR/MGE/defense results.
- `data/script_inputs/`: additional input tables for Figures 3, 5, and 6.
- `metadata/`: cohort identifiers and Supplementary Tables S1-S6.
- `figure_assets/`: editable source for the Figure 1A workflow.
- `docs/`: data dictionary and notes on reproducing the figures.

Genome assemblies can be downloaded from NCBI GenBank using the accessions in `metadata/Supplementary_Table_S1_main_cohort_metadata.xlsx`.

The files in `data/figure_source_data/` contain the values used in the figures. See `docs/REPRODUCIBILITY.md` for the main input files.

## Software

The analyses were run in R 4.5.2. The main software versions are listed in the manuscript. R packages can be installed with `code/install_packages.R`.

## Citation

Citation information will be added after publication.
