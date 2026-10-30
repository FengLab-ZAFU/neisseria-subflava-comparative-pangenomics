# Release validation

The current revision is dated 30 October 2026. Fresh checks for this revision are recorded below. The following execution record is the earlier 23 October 2026 validation, performed using an independently prepared workspace; it is not evidence that withdrawn inference was restored.

- All 86 Python, 37 R and 25 shell sources passed syntax checks.
- Stored and restored checksums were checked during input preparation, including the split compressed protein FASTA. The release validator also checks the complete file inventory, cohort membership, manifest paths and the 100 MiB file-size limit.
- The count summary reconciled the 152/144/8 cohort, pangenome partitions, 962 single-copy core clusters, 38 tet(M) carriers, seven sequence types and 180 distinct associated clusters.
- Pangenome summaries, structure/PERMDISP, gene associations, AMR, defense/VFDB tests, GO/KEGG mapping and tests, MGE/GI count models and cached sequence summaries were run through the public entry points. GO/KEGG Fisher and BH calculations were also checked independently. tet(M) parsimony was exercised in that historical validation but has since been withdrawn from the manuscript and current launcher.
- Figures 1 and 3–6 and Supplementary Figures S1–S3 were regenerated from the prepared workspace. Figure 2 rendering uses the supplied tree and iTOL tracks.

The validation environment used Python 3.12, Biopython 1.88, SciPy 1.18.1, statsmodels 0.15.0 and R 4.5.2. Published software records are listed separately in SOFTWARE.tsv. Comparisons with supplied results are listed in [NUMERICAL_CHECKS.tsv](NUMERICAL_CHECKS.tsv). Rows were aligned by identifiers where needed. Numeric differences were checked against a tolerance of `1e-6 × (1 + |reference value|)`; text and identifiers were compared exactly. Execution checks and numerical comparisons are reported separately.

The raw-genome prediction pipeline, new BLAST searches, 1,000-resample Panstripe fits and positional randomization were not repeated for this check. Full upstream reconstruction and new interval randomizations require the external inputs described in REPRODUCIBILITY.md. The sequence test used cached BLAST tables and is labeled accordingly.

## Checks on 30 October 2026

The current input manifest restored 1,106 frozen files with stored and decompressed SHA256 verification into an independent workspace. The summary entry point reconciled the cohort, screening flow, pangenome partitions, carrier and association counts. The revision check recalculated the 152 retained genomes' single-reference eligibility from all 912 directional ANI records and checked the 38-carrier/57-locus classifications against the current S4 export. Figure 2 packaging preserved the original 157-tip tree's unrooted splits, branch lengths and support labels and verified its root split and tip annotations.

The updated Figure 1 and Supplementary Figure S2 sources ran successfully through the public launcher using R 4.5.2, and their rendered output was inspected. S2 retained 3 chromosome-verified carriers, 1 assembly-supported association, 22 same-contig MGE co-occurrences and 0 direct overlaps. The plot exports the accession-to-column mapping. New Python sources and the updated launcher were checked by the release validator; new or changed R sources were parsed and exercised. The complete release checksum inventory and source-provenance inventory were refreshed.

This validation did not repeat new recombination analyses, Panstripe fitting, tet(M) acquisition/loss inference, ML inference or participant dereplication. The extraction implementations are supplied, but their complete upstream workflow was not rerun. Minor local locale/font messages occurred during plotting; the exported PDF and image content was inspected.
