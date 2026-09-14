# Figure source data — main139 cohort (139 genomes)

All files in this directory are **copies** (never moved) of analysis outputs from
`results_main139/` (and the iTOL annotation files from `outputs_main139/figures/iTOL/`),
renamed with a `main139_` prefix. The 6 excluded boundary samples
(SAMN18012341, SAMN18012350, SAMN24604773, SAMN51436610, SAMN51436611, SAMN51436612)
do not appear in any of these files. Sample IDs are BioSample accessions.

The two files `pcoa_coordinates.csv` and `pam_cluster_assignments.csv` (without the
`main139_` prefix) are identical copies placed here by the population-structure step;
`main139_pcoa_coordinates.csv` / `main139_pam_cluster_assignments.csv` are the same
content with the standard prefix.

## Figure 1: ANI cohort definition

| File | Content | Origin |
|---|---|---|
| `main139_Figure1B_pairwise_ANI_matrix.tsv` | Uncapped symmetric 139 x 139 pairwise ANI matrix used for average-linkage clustering in Figure 1B; the first column contains row sample IDs | `results/fastani_associated_group/pairwise_ani.txt.matrix`, subset to the frozen main139 cohort |
| `main139_Figure1B_type_strain_ANI_coverage.tsv` | Per-genome ANI and conservative bidirectional minimum coverage against the *N. subflava* ATCC 49275, *N. perflava* CCUG 17915 and *N. flavescens* ATCC 13120 references, plus the selected genospecies and ANI margin | `results/type_strain_reclassification/reclassification.tsv`, subset to the frozen main139 cohort |

## Phylogenetic tree and iTOL annotations

| File | Content | Origin |
|---|---|---|
| `main139_144tip_rooted.nwk` | 144-tip ML tree (139 ingroup + 5 N. flava outgroup), outgroup-rooted; GTR+F+I+R10, UFBoot+SH-aLRT support | `results_main139/04_recombination/iqtree2/main139_GTRFIR10.144tip.rooted.tree` |
| `main139_139tip_midpoint.nwk` | 139-tip ingroup-only ML tree, midpoint-rooted (diameter 0.0957) | `results_main139/04_recombination/iqtree2/main139_GTRFIR10.139tip.midpoint.tree` |
| `main139_iTOL_labels_144tip.txt` | iTOL LABELS annotation (tip labels = BioSample IDs) | `outputs_main139/figures/iTOL/iTOL_labels_144tip.txt` |
| `main139_iTOL_labels_strain_names_144tip.txt` | iTOL LABELS annotation for the 144-tip rooted tree; display labels prioritise strain/isolate names and use BioSample IDs only as fallbacks or duplicate-name disambiguators | `outputs_main139/figures/iTOL/iTOL_labels_strain_names_144tip.txt` |
| `main139_iTOL_colorstrip_group_144tip.txt` | iTOL COLORSTRIP annotation: ANI-defined genospecies per tip (N_subflava / N_perflava / N_flavescens / N. flava outgroup) | `outputs_main139/figures/iTOL/iTOL_colorstrip_group_144tip.txt` |
| `main139_iTOL_colorstrip_collection_region_144tip.txt` | iTOL COLORSTRIP annotation for descriptive collection region: United States (66), East Asia (37), Europe (19), and Unknown/not applicable (22) | `outputs_main139/figures/iTOL/iTOL_colorstrip_collection_region_144tip.txt` |
| `main139_iTOL_colorstrip_PAM_cluster_144tip.txt` | iTOL COLORSTRIP annotation for accessory-genome PAM clusters C1-C6 (139 ingroup) plus the rooted outgroup category | `outputs_main139/figures/iTOL/iTOL_colorstrip_PAM_cluster_144tip.txt` |
| `main139_iTOL_colorstrip_tetM_144tip.txt` | Narrow iTOL COLORSTRIP annotation for tet(M): 44 positive main-cohort genomes highlighted; negatives left white | `outputs_main139/figures/iTOL/iTOL_colorstrip_tetM_144tip.txt` |
| `main139_iTOL_group_legend_144tip.txt` | Legend: group colours and tip counts (47 / 82 / 10 / 5) | `outputs_main139/figures/iTOL/iTOL_group_legend_144tip.txt` |

## PCoA / PAM cluster structure

| File | Content | Origin |
|---|---|---|
| `main139_pcoa_coordinates.csv` | PCoA coordinates (PCo1–PCo10, cmdscale on Jaccard distances of accessory PAV) + variance explained per axis + `pam_cluster` / `operational_group` columns; 139 rows. Here, accessory PAV is defined as the 9,167 Shell + Cloud PGAP2 clusters; the full PAV matrix is not used for these coordinates. (The raw `pcoa_coordinates.csv` in `results_main139/05_population_structure/` lacks the two annotation columns; this copy is the PAM-labelled plotting version.) | `results_main139/05_population_structure/pcoa_coordinates.csv` + `pam_cluster_assignments.csv` |
| `main139_pam_cluster_assignments.csv` | Columns: `sample`, `pam_cluster` (C1–C6), `operational_group`, `silhouette_width`; 139 rows | `results_main139/05_population_structure/pam_cluster_assignments.csv` |
| `main139_permdisp_results.tsv` | PERMDISP tests of multivariate dispersion for PAM clusters and ANI-defined genospecies, calculated from the same binary Jaccard distances used for Figure 3C-D | `results_main139/05_population_structure/permdisp_results.tsv` |
| `main139_Supplementary_Table_PAM_by_BioProject.tsv` | One row per PAM x BioProject combination, including within-PAM project fraction, genospecies composition, tet(M) prevalence, assembly composition and sample IDs | Main139 metadata, PAM assignments and tet(M)-positive sample list |

## Pangenome (PGAP2) composition, rarefaction, Heap law

| File | Content | Origin |
|---|---|---|
| `main139_pgap2_partition_summary_statistics.txt` | PGAP2 partition counts: strict core 1203 / core 187 / soft core 145 / shell 872 / cloud 8295 / total 10702, with % | `results_main139/01_pgap2/pgap2.partition.summary_statistics.txt` |
| `main139_postprocess.pan_group_stat.tsv` | Same partition counts in TSV (Group / Count / Proportion) | `results_main139/01_pgap2/postprocess/profile/postprocess.pan_group_stat.tsv` |
| `main139_postprocess.rarefaction.tsv` | Core- and pan-genome rarefaction curves (columns: Type / Strain / Sampling); 10 repetitions per strain count | `results_main139/01_pgap2/postprocess/profile/postprocess.rarefaction.tsv` |
| `main139_postprocess.curve_fit.txt` | Heap/power-law fits: pan fit `y = 2009.668 * x^0.340` (α = 0.340, 95% CI [0.618, 0.644] for the new-gene decay exponent 0.631 → open pangenome) | `results_main139/01_pgap2/postprocess/profile/postprocess.curve_fit.txt` |
| `main139_pav_matrix.tsv.gz` | Gene presence/absence matrix, 10702 clusters × 139 samples (1/0; gzip of the `.pav` file) | `results_main139/01_pgap2/pgap2.partition.gene_content.pav` |

## Lineage-associated genes

| File | Content | Origin |
|---|---|---|
| `main139_pam_associated_genes.csv` | Significant PAM-cluster-associated gene clusters (one-vs-rest Fisher, BH-FDR < 0.05, in-group prevalence > 80%, out-group < 20%); 935 associations / 914 unique clusters | `results_main139/06_lineage_association/pam_associated_genes.csv` |
| `main139_taxon_associated_genes.csv` | Significant operational-group-associated clusters; 99 associations / 99 unique clusters | `results_main139/06_lineage_association/taxon_associated_genes.csv` |

Columns (both files): `scheme, clust, group, present_in_group, group_size,
present_in_rest, rest_size, odds_ratio, p_value, prev_in_group, prev_in_rest, fdr,
product, gene_name, pgap2_group`.

## COG functional composition / enrichment

| File | Content | Origin |
|---|---|---|
| `main139_cog_class_composition.tsv` | COG class counts per partition (core/shell/cloud columns) | `results_main139/07_functional_enrichment/cog_class_composition.tsv` |
| `main139_cog_by_class_enrichment.csv` | COG-class enrichment per partition (core/shell/cloud): count, odds ratio, p, FDR | `results_main139/07_functional_enrichment/cog_by_class_enrichment.csv` |
| `main139_cog_by_pam_group_enrichment.csv` | COG-class enrichment among PAM-associated genes per PAM group (C6 × X significant: OR=2.21, FDR=3.44e-03) | `results_main139/07_functional_enrichment/cog_by_pam_group_enrichment.csv` |

## AMR / tet(M) / AMR–MGE co-localization

| File | Content | Origin |
|---|---|---|
| `main139_amr_per_sample_summary.tsv` | Per-sample RGI and abritamr hits (sample_id, rgi_n_hits, rgi_aros, abritamr_n_genes, abritamr_hits); 139 rows | `results_main139/08_amr_mge/amr_summary/amr_per_sample_summary.tsv` |
| `main139_amr_concordance_summary.tsv` | Tool concordance: both=48, only abritamr=8, only RGI=6, neither=77 (total 139) | `results_main139/08_amr_mge/amr_summary/concordance_summary.tsv` |
| `main139_tetm_by_pam_group.tsv` | tet(M) carrier counts per PAM group (C1 21/72 … C5 10/10, C6 0/9) with sample lists | `results_main139/08_amr_mge/tetM_lineage_association/tetm_by_pam_group.tsv` |
| `main139_tetm_by_operational_group.tsv` | tet(M) carriers per operational group (N_flavescens 0/10, N_perflava 30/82, N_subflava 14/47) | `results_main139/08_amr_mge/tetM_lineage_association/tetm_by_operational_group.tsv` |
| `main139_key_gene_summary.tsv` | AMR–MGE co-localization per key gene hit (sample, gene, tools, contigs, MGE types, best evidence tier: sequence_only / mge_colocated / verified_complete) | `results_main139/08_amr_mge/amr_mge_colocalization/key_gene_summary.tsv` |

## MGE / genomic islands

| File | Content | Origin |
|---|---|---|
| `main139_cluster_mge_status.tsv` | Per-cluster MGE overlap status (pgap2_group, is_mge_associated, n_mge_members, MGE type counts); 10702 rows | `results_main139/08_amr_mge/mge_pav_overlap/cluster_mge_status.tsv` |
| `main139_mge_regions_per_sample.tsv` | All MGE regions per sample (source_tool, mge_type, contig, start, end, detail); 3562 rows | `results_main139/08_amr_mge/mge_pav_overlap/mge_regions_per_sample.tsv` |
| `main139_genomic_islands_per_sample.tsv` | IslandPath-DIMOB islands per sample (n_islands, gi_bp, genome_bp, gi_fraction, assembly info, PAM cluster, operational group); 1031 islands over 139 samples | `results_main139/11_genomic_islands/genomic_islands_per_sample.tsv` |

## Defense systems (DefenseFinder)

| File | Content | Origin |
|---|---|---|
| `main139_defense_finder_per_sample_summary.tsv` | Per-sample complete defense systems (n_hmmer_hits, n_complete_systems, n_system_types, system types list); 139 rows | `results_main139/09_defense_systems/defense_finder_per_sample_summary.tsv` |
| `main139_defense_by_pam_group.tsv` | Per PAM group: total/mean/median/min–max complete systems and per-type counts | `results_main139/09_defense_systems/defense_by_pam_group.tsv` |
| `main139_defense_by_operational_group.tsv` | Same per operational group | `results_main139/09_defense_systems/defense_by_operational_group.tsv` |

## Supplementary figures

- The Panstripe accumulation analysis remains a supplementary figure in
  `outputs_main139/figures/figure3/` with its plotting source and TSV source data.
- The per-genome VFDB candidate-homolog profile and all figure-level source tables
  are collected in `outputs_main139/figures/supplementary_vfdb/`. These are
  sequence-level homology calls and do not establish virulence or phenotype.

Generated: 2026-09-07. All numbers in these files are derived from the 139-genome
manuscript main cohort; per-module details and old-145-cohort comparisons are in
`results_main139/*/` summary markdown files.
