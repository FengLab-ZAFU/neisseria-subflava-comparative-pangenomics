# Reproducing figures

Install the R packages listed in `code/install_packages.R`. The figure scripts and their input tables are organized by manuscript figure. `data/script_inputs/` contains the additional metadata and per-sample tables used by Figures 3, 5, and 6.

The principal inputs are:

- Figure 1B ANI values: `data/figure_source_data/main139_Figure1B_pairwise_ANI_matrix.tsv` and `main139_Figure1B_type_strain_ANI_coverage.tsv`.
- Figure 2 tree and annotations: `main139_144tip_rooted.nwk` plus the `main139_iTOL_*_144tip.txt` files.
- Figure 3 pangenome and PCoA data: `main139_postprocess.*`, `main139_pcoa_coordinates.csv`, and `main139_pav_matrix.tsv.gz`.
- Figures 4-6: the corresponding PAM, COG, AMR, MGE, and defense tables in the same directory.
- Supplementary Tables S1-S6: `metadata/`.

The deposited scripts retain the directory structure used for the analysis. When running them in another environment, set the input paths at the beginning of each script to the corresponding files in this repository.
