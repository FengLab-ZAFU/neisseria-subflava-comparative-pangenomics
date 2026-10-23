# Data dictionary

`metadata/overall152.tsv`: 152 representative genomes, reference assignments, unresolved status, source studies and sensitivity labels. `metadata/comparison144.tsv`: the 23/114/7 primary comparison. `metadata/unresolved8.txt`: genomes excluded from reference-group comparisons. `metadata/screening1303.tsv`: one row per screened assembly and its final disposition. Study/participant labels follow frozen provenance decisions; genomes are not assumed fully independent.

`data/tables/S1_*`: metadata, screening ledger, boundary dDDH, clone thresholds and source composition.
`S2_*`: dictionary sensitivity, repertoire summaries, distance tests and Panstripe models/contrasts.
`S3_*`: association-eligible universe, all gene tests, candidates, COG and exploratory KEGG results.
`S4_*`: per-genome AMR, mapped tet(M) loci, seven sequence types, reference sites and conditional gain/loss summaries.
`S5_*`: per-genome MGE/GI/defense observations, interval overlap/null tests and group comparisons.
`S6_*`: candidate VFDB homolog presence and group tests.

The overall PGAP2 dictionary has 13,215 clusters. The ≥99% core (1,347 clusters) differs from the strict all-genome core (1,118) and the single-copy alignment set (962). Shell/cloud features are defined in all 152 genomes even when association tests use only 144. There are 181 group–cluster associations representing 180 distinct clusters.

tet(M) has 38 positive genomes, 57 mapped loci and 52 retained full-length sequences representing seven sequence types. Genomes may carry multiple loci/types; counts are not interchangeable. GI/CDS-instance proportions, cluster-level MGE overlap and randomization-eligible denominators likewise use different units, as labeled in the tables.

CSV/TSV files use UTF-8 and retain identifiers. Empty/NA values are not silently converted to biological absence. Frozen input data preserve the original project table names for dependency resolution; versioned metadata and current supplement exports define the final manuscript scopes.

Functional mapping tables distinguish an exact annotation match, a sequence input without an annotation row, and an unverified source sequence. Missing annotations are treated as missing information. GO/KEGG tests use their stated annotated universes, rather than all clusters as an automatically annotated background.
