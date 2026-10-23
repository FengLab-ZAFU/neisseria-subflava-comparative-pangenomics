# Reproduction instructions

## Prepare the workspace

Use Python 3.10 or later. Run the following from the repository root:

```bash
python code/verify_release.py
python code/reproduce.py prepare
python code/reproduce.py run summary
```

Preparation checks the stored and decompressed SHA256 values in `data/frozen_manifest.tsv`. Large compressed files are split into ordered parts and joined before decompression. Files are restored atomically into `work/project/`; the inputs under `data/` are preserved. Repeating `prepare` restores modified analytical files to their supplied versions and refreshes the scripts, so use a separate `--project /absolute/workspace/path` for each independent rerun.

`summary` recalculates counts from copy-number matrices and curated tables. Screening totals are reconciled from the recorded dispositions; this step does not repeat ANI screening or SNP clustering. [VALIDATION.md](VALIDATION.md) describes the checks performed for this release.

## Software

Python analysis dependencies are listed in `requirements.txt`. R packages vary by stage:

| Stage | Main R packages |
|---|---|
| Structure and PERMDISP | vegan, permute |
| Association, functional tests and defense | stats (included with R) |
| Panstripe | panstripe 0.4.0, ape, phangorn, jsonlite, ragg |
| tet(M) parsimony | ape, phangorn |
| MGE/GI count models | dplyr, tidyr, jsonlite |
| Figures and sequence displays | ggplot2, patchwork, dplyr, tidyr, ape, phangorn, Biostrings, ggtext, svglite, ragg, xml2, jsonlite, pdftools |
| Additional figure checks | pdftools, png; pwalign for sequence alignment |

Panstripe was run at commit `7843a8d6c25269dc1ae30d52480d45274fa51092`. The original main workflow recorded R 4.3.3; the corrected PERMDISP stage recorded R 4.5.3, vegan 2.7.5 and permute 0.9.10. Local Panstripe, tet(M) and final figure analyses used R 4.5.2. Software records are listed in [SOFTWARE.tsv](SOFTWARE.tsv), with session records among the supplied data.

R and Python analyses can be selected separately. Use `--rscript /path/to/Rscript` to select an R installation. BLAST+ is required for `tetm-sequences`, which runs both translated searches before summarizing their outputs. `tetm-cached-summary` only summarizes the supplied search tables. Figure rendering also requires suitable fonts and graphics devices. Upstream tool stages use Linux shell utilities and their recorded conda environments; structure and PERMDISP use `sha256sum`.

## Run analyses

`python code/reproduce.py list` lists the available stages. `run STAGE --dry-run` displays commands without executing them. All stages use a prepared workspace; `NEISSERIA_PROJECT_ROOT`, `PYTHON` and `RSCRIPT` are set by the launcher. For a native script, set these variables explicitly and use the working directory specified in [ANALYSIS_STAGES.md](ANALYSIS_STAGES.md).

Current group comparisons and plots share `results_revision_20261019_boundary8_analysis/`. The `permdisp` stage uses the archived permutation matrices, writes detailed residual-permutation results separately, and updates the current distance-test table only after successful completion. Panstripe reruns write to `analysis_outputs/panstripe/`; supplied model objects remain under `results_panstripe_minimal_20261020/models/`.

The upstream main152 scripts contain the initial 31/114/7 ANI grouping. Their tree, dictionary and annotations are shared inputs. Final comparisons use the boundary8 scripts and the 23/114/7 assignments; the initial grouping is used only for the specified sensitivity analyses.

## External inputs

Raw genome FASTA and annotation files, large databases, original genus-query JSON pages and some corrected raw ANI/SNP files are not bundled. Source-correction code also reads intermediate ledgers from earlier stages. The corrected source, canonical-unit and cohort tables are supplied, but a full upstream rerun requires recovery or regeneration of these intermediate inputs. Versioned assembly accessions are provided for genome retrieval.

Set the following variables when running the relevant upstream scripts:

| Variable | Required input |
|---|---|
| `CONDA_SH` | Conda initialization script; default `$HOME/miniconda3/etc/profile.d/conda.sh` |
| `BAKTA_DB` | Bakta v6.0 full database |
| `CHECKM2_DB` | CheckM2 `uniref100.KO.1.dmnd` file |
| `GENOMAD_DB` | geNomad database |
| `ISESCAN_SCRIPT` | ISEScan `isescan.py` |
| `EGGNOG_SIF`, `EGGNOG_DB` | emapper 3.0.0-beta6 image and eggNOG 7.0.0 database |
| `APPTAINER`, `SKANI`, `DEFENSE_FINDER` | Optional executable overrides; otherwise resolved on PATH |
| `ANNOTATION_PENDING_TSV` | Original annotation-pending accession table |
| `NEISSERIA_SOURCE_ROOT` | Optional root for recovered external genome/annotation files |
| `BAKTA_RECORDED_DB_PATH` | Database path recorded in the original logs, for the cache-reuse check |
| `DIMOB_LIST`, `GFF3_TO_GBK` | DIMOB accession/FASTA list and original GFF3-to-per-contig-GenBank conversion helper |

The original annotation-pending table and DIMOB conversion helper were unavailable in the retained files and must be recovered before rerunning those upstream stages. This does not affect the supplied GI predictions or downstream interval analyses. Database versions and records should be checked before comparing new tool predictions with the published results.

## Positional randomization inputs

The supplied interval-overlap tables, eligibility records and null-test results support checking the published summaries. A new positional-randomization run additionally needs five normalized tables under `results_revision_20261015_main152_species/09_interval_overlap/tables/`: `contigs.tsv`, `genes.tsv`, `mge_local.tsv`, `gi.tsv` and `whole_contig.tsv`. These files were absent from the retained source directories and returned archive. Recover them from the original analysis workspace, or regenerate them with `scripts_main152/m8_build_intervals.py` after restoring the raw genomes and prediction outputs. The launcher checks these requirements before starting `overlap`.

The `defense` stage recalculates defense and VFDB tests from the supplied per-genome tables. Its descriptive overlap summaries reuse the supplied gene/GI overlap observations; they do not rerun interval detection or positional randomization.

Historical environment, preflight and source-audit records are retained as provenance. They may name analyses omitted from the current manuscript; the current clustering comparison is PGAP2 orth_id 0.98 versus 0.95, as reported in Supplementary Table S2.
