# Neisseria subflava-associated group pangenome analysis

Analysis code and data for **Pangenome analysis of the Neisseria subflava-associated group reveals variation in accessory gene content and tet(M) carriage**.

The study includes 152 genomes. Comparisons among the three reference-assigned groups use 144 genomes: 23 subflava-, 114 perflava- and 7 flavescens-reference-assigned genomes. Eight species-unresolved genomes are included in the overall analyses. Group repertoires are derived from a shared PGAP2 gene-cluster dictionary.

## Getting started

From the repository root:

```bash
python code/verify_release.py
python code/reproduce.py prepare
python code/reproduce.py run summary
python code/reproduce.py run revision2-checks
python code/reproduce.py run figure2
python code/reproduce.py list
```

These commands require Python 3.10 or later and use only its standard library. Preparation restores the supplied files and copies the analysis scripts to `work/project/`. Allow at least 1 GB of free space. The summary checks cohort membership, pangenome partitions, single-copy core selection, tet(M) carriage and associated-cluster counts against the supplied tables.

Install the dependencies for the selected analysis before running further stages:

```bash
python -m pip install -r requirements.txt
python code/reproduce.py run association
python code/reproduce.py run mge-gi --rscript Rscript
python code/reproduce.py run figure3 --rscript Rscript
```

Each command runs one stage in the prepared workspace. See [analysis stages](docs/ANALYSIS_STAGES.md) for their inputs and execution order, and [reproduction instructions](docs/REPRODUCIBILITY.md) for software requirements.

## Repository contents

| Directory | Contents |
|---|---|
| `code/` | Preparation, validation and analysis scripts |
| `data/frozen/` | Analytical inputs, result tables, trees, model objects and cached sequence records |
| `data/supplementary_tables/` | Supplementary Tables S1–S6 |
| `data/tables/` | TSV exports of the supplementary worksheets |
| `metadata/` | Cohort lists, screening decisions and study composition |
| `docs/` | Methods, software records, file checksums and script provenance |

This revision was synchronized with Supplementary Tables S1–S6 on 30 October 2026. Figure 2 uses the supplied 157-tip tree with five N. flava outgroups. The 152-genome study cohort and 144-genome primary comparisons are unchanged. `revision2-checks` reconciles directional ANI scores and current locus evidence; `figure2` verifies and packages the supplied iTOL tree and tracks. It does not infer a new tree.

Panstripe inference and tet(M) gain/loss reconstruction have been withdrawn from the manuscript and removed from the current launcher. Historical implementations and frozen outputs remain available as provenance; their scope is recorded in [SCRIPT_SCOPE.tsv](docs/SCRIPT_SCOPE.tsv). They do not support current conclusions. The retained phylogenies and descriptive patristic-distance comparisons are based on alignments without recombination filtering.

The release supports downstream reanalysis from the supplied inputs. Full reconstruction from assemblies also requires external genome files, tool databases and several historical ANI/SNP intermediates. Repeating positional randomization requires the original normalized interval tables. Details are given in [retrieval and curation](docs/RETRIEVAL.md) and the [revision notes](docs/REVISION_20261030.md).
