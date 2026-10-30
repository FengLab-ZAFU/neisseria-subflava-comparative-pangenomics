# Analysis stages

Paths below are relative to the prepared workspace. The launcher runs only the selected stage. Supplied intermediate tables allow most downstream stages to be run independently; the order below describes a sequential reanalysis.

| Launcher stage | Calculation and inputs |
|---|---|
| `summary` | Count reconciliation from the supplied matrices, cohort lists and supplementary tables |
| `revision2-checks` | Check current table snapshots, single-reference membership from 912 directional ANI records, and 38-carrier/57-locus evidence |
| `pangenome` | Group repertoires and accumulation curves from the shared PGAP2 dictionary |
| `structure` | PCoA, PERMANOVA and residual-permutation PERMDISP from the supplied distance matrices |
| `permdisp` | Recalculate PERMDISP using the saved permutation draws and update the current distance-test table |
| `association` | Gene associations and COG tests from the binary matrix and 144-genome assignments |
| `functional` | Exact protein-sequence mapping, GO ancestry and GO/KEGG tests; run after `association` |
| `amr` | Group AMR and tet(M) summaries from per-genome evidence |
| `defense` | Defense/VFDB tests and descriptive summaries using the supplied overlap observations |
| `overlap` | Interval randomization and scope annotations; requires external normalized interval tables and the `defense` stage |
| `tetm-sequences` | Reference-site comparison, flanking sequences, two new tblastn searches and marker summaries |
| `tetm-cached-summary` | Marker summaries from the supplied search results, without new searches |
| `mge-gi` | Six-endpoint MGE/GI count tests, assembly-adjusted models and representative-subset checks |
| `figure1` | Cohort flow from the supplied screening and group metadata |
| `figure2` | Verify and package the supplied 157-tip outgroup-rooted tree and iTOL tracks; no ML inference is executed |
| `figure3` | Repertoire, ordination and accumulation panels |
| `figure4` | Associated-cluster matrix, COG profiles and KEGG mapping; uses `association` and `functional` outputs |
| `figure5` | AMR, tet(M) and seven-type sequence panels |
| `figure6` | MGE/GI/defense panels; uses the count-model outputs |
| `supplementary1` | ANI-margin distribution for all 169 dereplicated candidates |
| `supplementary2` | Current four-endpoint assembly-evidence plot from S4 and the original 38-carrier column order; exports the accession-to-column mapping |
| `supplementary3` | VFDB homolog prevalence; verifies raw hit calls and the 42 eligible group tests |

The final comparison scripts are under `results_revision_20261019_boundary8_analysis/scripts/`. The launcher uses that directory as their working directory and the workspace root for the other stages. The Figure 1 script receives `manuscript/Figures/new/` as its first argument and `--fig1-only` as its second. The `figure2` stage writes a checked package to `analysis_outputs/figure2_itol/`; render its supplied tree and annotation tracks in iTOL as an outgroup-rooted cladogram. The `supplementary2` stage uses `revision2/plot_supplementary2.R`.

Panstripe and tet(M) acquisition/loss reconstruction are withdrawn analyses and have no current launcher stage. Older sources are retained for provenance, as recorded in [SCRIPT_SCOPE.tsv](SCRIPT_SCOPE.tsv). The coordinate-preserving outgroup extraction sources are included under `results_review2_20261028_v2/02_rooting/outgroup_orthologs/scripts/`; see [revision notes](REVISION_20261030.md) for their dependencies.

For individual upstream scripts, use the workspace root and set `NEISSERIA_PROJECT_ROOT` to its absolute path:

| Upstream step | Script directory |
|---|---|
| Retrieval, filtering, downloads, ANI and QC | `scripts/revision_20261001/a2_*` |
| Pairwise SNP distances | `scripts/revision_20261001/a3_*` |
| Source correction and canonical units | `results_revision_20261010_sourcefix_v7/delivery/scripts/` |
| Annotation and cache binding | `scripts/revision_20261001/b1_bakta.sh`, `scripts_phaseB_v7/p1_*` |
| PGAP2 and strict single-copy core tree | `scripts_main152/m2_*`, `scripts_main152/m3*` |
| AMR, MGE and defense prediction | `scripts/revision_20261001/d*`, `scripts_main152/m7c_vfdb_run.sh` |

Source-correction scripts retain version identifiers because later steps read specific earlier tables. They are upstream records, rather than a self-contained automated pipeline. Required external inputs are listed in [REPRODUCIBILITY.md](REPRODUCIBILITY.md). Figure and data directory names are kept stable for dependency resolution; supplementary table numbering follows the six workbooks in `data/supplementary_tables/`.
