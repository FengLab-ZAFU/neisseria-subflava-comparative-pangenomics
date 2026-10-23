# Release validation

Checked on 23 October 2026 using an independently prepared workspace.

- All 86 Python, 37 R and 25 shell sources passed syntax checks.
- Stored and restored checksums were checked during input preparation, including the split compressed protein FASTA. The release validator also checks the complete file inventory, cohort membership, manifest paths and the 100 MiB file-size limit.
- The count summary reconciled the 152/144/8 cohort, pangenome partitions, 962 single-copy core clusters, 38 tet(M) carriers, seven sequence types and 180 distinct associated clusters.
- Pangenome summaries, structure/PERMDISP, gene associations, AMR, defense/VFDB tests, GO/KEGG mapping and tests, MGE/GI count models, tet(M) parsimony and cached sequence summaries were run through the public entry points. GO/KEGG Fisher and BH calculations were also checked independently.
- Figures 1 and 3–6 and Supplementary Figures S1–S3 were regenerated from the prepared workspace. Figure 2 rendering uses the supplied tree and iTOL tracks.

The validation environment used Python 3.12, Biopython 1.88, SciPy 1.18.1, statsmodels 0.15.0 and R 4.5.2. Published software records are listed separately in SOFTWARE.tsv. Comparisons with supplied results are listed in [NUMERICAL_CHECKS.tsv](NUMERICAL_CHECKS.tsv). Rows were aligned by identifiers where needed. Numeric differences were checked against a tolerance of `1e-6 × (1 + |reference value|)`; text and identifiers were compared exactly. Execution checks and numerical comparisons are reported separately.

The raw-genome prediction pipeline, new BLAST searches, 1,000-resample Panstripe fits and positional randomization were not repeated for this check. Full upstream reconstruction and new interval randomizations require the external inputs described in REPRODUCIBILITY.md. The sequence test used cached BLAST tables and is labeled accordingly.
