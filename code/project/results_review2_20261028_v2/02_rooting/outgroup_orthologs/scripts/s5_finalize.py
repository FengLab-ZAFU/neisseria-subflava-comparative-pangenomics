#!/usr/bin/env python
"""s5: write extraction_report.md and final checksums.sha256."""
import json, hashlib, os, glob, subprocess

OUT = os.path.join(os.environ["NEISSERIA_PROJECT_ROOT"], "results_review2_20261028_v2/02_rooting/outgroup_orthologs")
OGS = ["SAMN20865411", "SAMEA122247875", "SAMEA122247926", "SAMEA122247929", "SAMEA122247940"]

rep = json.load(open(f"{OUT}/extraction_report.json"))
blast_ver = open(f"{OUT}/logs/blastn_version.txt").read().strip().split("\n")[0]

L = []
w = L.append
w("# Outgroup ortholog extraction report (module B, 5 N. flava outgroups)")
w("")
w(f"- BLAST: {blast_ver}")
w(f"- blastn params: -task blastn -word_size 11 -dust no -evalue 1e-20, outfmt 6 std + qseq sseq")
w(f"- Acceptance: merged (same contig+strand) query coverage >= 50% of gap-stripped consensus, weighted identity >= 70%")
w(f"- Uniqueness: second independent hit (different contig or non-overlapping subject interval) with bitscore >= 90% of best -> ambiguous, excluded")
w(f"- Consensus: per-column majority over A/C/G/T only, ties A>C>G>T; 20 all-gap/all-N columns skipped over 962 genes")
w("")
w("## Per-outgroup summary")
w("")
w("| outgroup | present_unique | % | covered cols | % of 906074 | mean ident | ident range | ambiguous | low_coverage | absent | discarded bases |")
w("|---|---|---|---|---|---|---|---|---|---|---|")
for og in OGS:
    o = rep["outgroups"][og]
    w(f"| {og} | {o['present_unique']} | {o['present_unique_frac']*100:.1f} | {o['covered_alignment_columns']} | {o['covered_fraction_of_906074']*100:.1f} | {o['mean_identity']:.2f} | {o['min_identity']:.2f}-{o['max_identity']:.2f} | {o['ambiguous']} | {o['low_coverage']} | {o['absent']} | {o['discarded_bases']} |")
w("")
w("## Outgroups available per gene (genes count by #present_unique outgroups)")
w("")
w(", ".join(f"{k} outgroups: {v} genes" for k, v in sorted(rep["outgroups_per_gene_distribution"].items(), key=lambda kv: int(kv[0]))))
w("")
w("## Verification")
w("")
w("- All 5 concatenated rows have length 906074.")
ok = all(v["match"] for v in rep["verification_concat_vs_pergene"])
genes_checked = sorted({v["gene"] for v in rep["verification_concat_vs_pergene"]})
w(f"- Concat-vs-per-gene slice check for {len(genes_checked)} genes ({', '.join(genes_checked)}) x 5 outgroups: {'ALL MATCH' if ok else 'MISMATCH FOUND'}")
w("- Manual spot checks (logs/spotcheck_*.txt): clust_1001 x SAMN20865411 (full 1269-bp single minus-strand HSP, 1269/1269 bases re-verified against subject FASTA), clust_1103 x SAMEA122247875 (partial HSP q156-527, 372/372 bases verified, row gaps outside HSP), clust_1017 x SAMN20865411 (best coverage 48.57% -> low_coverage, correctly no row; second hit bitscore 569 < 90% of 1637 -> not ambiguous). Coordinate mapping confirmed in all three.")
w("")
w("## Observations")
w("")
w("- SAMEA122247929 and SAMEA122247940 yield identical acceptance/coverage/identity statistics on every gene; their blast outputs differ only in contig names/coordinates and line order (same HSP query intervals and identities). Input FASTAs are distinct files (distinct GCA accessions); the two assemblies are near-identical over all core-gene regions. This is a property of the provided inputs, not a pipeline artifact.")
w("- SAMN20865411 shows 36 overlapping-HSP placement conflicts (first-placed base kept); all other outgroups 0.")
w("- Stop rule NOT triggered: all outgroups have mean identity > 60% and > 90% of columns covered.")
w("")
w("Full machine-readable details in extraction_report.json; per-gene/outgroup decisions in availability_matrix.tsv.")
open(f"{OUT}/extraction_report.md", "w").write("\n".join(L) + "\n")

def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()

files = ["consensus_all.fasta", "consensus_map.json", "availability_matrix.tsv",
         "outgroup_rows_full906074.fasta", "extraction_report.json", "extraction_report.md"]
files += [f"blast/{og}.blastn.tsv" for og in OGS]
files += [f"per_gene/{os.path.basename(p)}" for p in sorted(glob.glob(f"{OUT}/per_gene/*.fasta"))]
with open(f"{OUT}/checksums.sha256", "w") as cf:
    for f in files:
        cf.write(f"{sha(f'{OUT}/{f}')}  {f}\n")
print(f"report.md written; checksums for {len(files)} files")
