# Outgroup ortholog extraction report (module B, 5 N. flava outgroups)

- BLAST: blastn: 2.17.0+
- blastn params: -task blastn -word_size 11 -dust no -evalue 1e-20, outfmt 6 std + qseq sseq
- Acceptance: merged (same contig+strand) query coverage >= 50% of gap-stripped consensus, weighted identity >= 70%
- Uniqueness: second independent hit (different contig or non-overlapping subject interval) with bitscore >= 90% of best -> ambiguous, excluded
- Consensus: per-column majority over A/C/G/T only, ties A>C>G>T; 20 all-gap/all-N columns skipped over 962 genes

## Per-outgroup summary

| outgroup | present_unique | % | covered cols | % of 906074 | mean ident | ident range | ambiguous | low_coverage | absent | discarded bases |
|---|---|---|---|---|---|---|---|---|---|---|
| SAMN20865411 | 894 | 92.9 | 821421 | 90.7 | 83.37 | 70.08-98.37 | 0 | 36 | 32 | 1804 |
| SAMEA122247875 | 876 | 91.1 | 802482 | 88.6 | 83.34 | 70.16-98.37 | 1 | 36 | 49 | 1800 |
| SAMEA122247926 | 885 | 92.0 | 814693 | 89.9 | 83.34 | 70.16-98.37 | 0 | 39 | 38 | 1803 |
| SAMEA122247929 | 889 | 92.4 | 820039 | 90.5 | 83.38 | 70.16-98.37 | 0 | 37 | 36 | 1819 |
| SAMEA122247940 | 889 | 92.4 | 820039 | 90.5 | 83.38 | 70.16-98.37 | 0 | 37 | 36 | 1819 |

## Outgroups available per gene (genes count by #present_unique outgroups)

0 outgroups: 68 genes, 1 outgroups: 4 genes, 2 outgroups: 1 genes, 3 outgroups: 3 genes, 4 outgroups: 12 genes, 5 outgroups: 874 genes

## Verification

- All 5 concatenated rows have length 906074.
- Concat-vs-per-gene slice check for 3 genes (clust_1001, clust_3682, clust_9978) x 5 outgroups: ALL MATCH
- Manual spot checks (logs/spotcheck_*.txt): clust_1001 x SAMN20865411 (full 1269-bp single minus-strand HSP, 1269/1269 bases re-verified against subject FASTA), clust_1103 x SAMEA122247875 (partial HSP q156-527, 372/372 bases verified, row gaps outside HSP), clust_1017 x SAMN20865411 (best coverage 48.57% -> low_coverage, correctly no row; second hit bitscore 569 < 90% of 1637 -> not ambiguous). Coordinate mapping confirmed in all three.

## Observations

- SAMEA122247929 and SAMEA122247940 yield identical acceptance/coverage/identity statistics on every gene; their blast outputs differ only in contig names/coordinates and line order (same HSP query intervals and identities). Input FASTAs are distinct files (distinct GCA accessions); the two assemblies are near-identical over all core-gene regions. This is a property of the provided inputs, not a pipeline artifact.
- SAMN20865411 shows 36 overlapping-HSP placement conflicts (first-placed base kept); all other outgroups 0.
- Stop rule NOT triggered: all outgroups have mean identity > 60% and > 90% of columns covered.

Full machine-readable details in extraction_report.json; per-gene/outgroup decisions in availability_matrix.tsv.
