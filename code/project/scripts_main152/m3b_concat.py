#!/usr/bin/env python3
"""M3a (cont.): concatenate per-gene MAFFT alignments -> core alignment.

- cluster order: lexicographic cluster ID; sequences per gene in roster order
  (MAFFT --inputorder preserves input order; verified explicitly).
- all-gap columns: reported; removed only if they exist, with before/after
  lengths and column mapping . Invariant sites retained.
- outputs: core_alignment.fasta, core_alignment.partition.nex (1-based closed
  intervals), core_gene_lengths.tsv, alignment_QC.tsv (per-genome ACGT/N/gap),
  site_summary.tsv (constant/variable/parsimony-informative, identical-seq groups)
"""
import json
import os
import sys
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import RESULTS, read_roster, write_tsv

CORE = os.path.join(RESULTS, "03_phylogeny/strict_singlecopy_core")


def read_aln(path):
    seqs, name, buf = [], None, []
    for ln in open(path):
        if ln.startswith(">"):
            if name is not None:
                seqs.append((name, "".join(buf).upper()))
            name, buf = ln[1:].strip().split()[0], []
        else:
            buf.append(ln.strip())
    if name is not None:
        seqs.append((name, "".join(buf).upper()))
    return seqs


def main():
    roster = read_roster()
    genes = sorted(fn[:-4] for fn in os.listdir(os.path.join(CORE, "per_gene_fastas")) if fn.endswith(".fna"))
    aln_dir = os.path.join(CORE, "per_gene_alignments")
    missing = [g for g in genes if not os.path.isfile(os.path.join(aln_dir, g + ".aln.fna"))]
    if missing:
        sys.exit(f"FATAL: {len(missing)} alignments missing, e.g. {missing[:5]}")

    per_gene = {}
    problems = []
    for g in genes:
        seqs = read_aln(os.path.join(aln_dir, g + ".aln.fna"))
        names = [n for n, _ in seqs]
        if names != roster:
            problems.append(f"{g}: sequence order/content != roster ({len(names)} seqs)")
            continue
        Ls = {len(s) for _, s in seqs}
        if len(Ls) != 1:
            problems.append(f"{g}: unequal alignment lengths {sorted(Ls)[:5]}")
            continue
        per_gene[g] = (seqs, Ls.pop())
    if problems:
        print("\n".join(problems[:20]))
        sys.exit(1)

    # all-gap column check (per gene)
    allgap = {}
    for g, (seqs, L) in per_gene.items():
        idx = [i for i in range(L) if all(s[i] == "-" for _, s in seqs)]
        if idx:
            allgap[g] = idx
    if allgap:
        with open(os.path.join(CORE, "allgap_column_map.tsv"), "w") as fh:
            fh.write("cluster\tlength_before\tlength_after\tremoved_0based_columns\n")
            for g, idx in allgap.items():
                fh.write(f"{g}\t{per_gene[g][1]}\t{per_gene[g][1] - len(idx)}\t{';'.join(map(str, idx))}\n")
        for g, idx in allgap.items():
            seqs, L = per_gene[g]
            keep = [i for i in range(L) if i not in set(idx)]
            per_gene[g] = ([(n, "".join(s[i] for i in keep)) for n, s in seqs], L - len(idx))
        print(f"all-gap columns removed in {len(allgap)} genes (map written)")
    else:
        print("no all-gap columns")

    # concatenate
    gene_order = sorted(per_gene)
    concat = {acc: [] for acc in roster}
    bounds = []
    pos = 1
    for g in gene_order:
        seqs, L = per_gene[g]
        for n, s in seqs:
            concat[n].append(s)
        bounds.append((g, pos, pos + L - 1, L))
        pos += L
    total_len = pos - 1

    with open(os.path.join(CORE, "core_alignment.fasta"), "w") as fh:
        for acc in roster:
            s = "".join(concat[acc])
            assert len(s) == total_len
            fh.write(f">{acc}\n")
            for i in range(0, len(s), 60):
                fh.write(s[i:i + 60] + "\n")

    with open(os.path.join(CORE, "core_alignment.partition.nex"), "w") as fh:
        fh.write("#nexus\nbegin sets;\n")
        for g, a, b, L in bounds:
            fh.write(f"  charset {g} = {a}-{b};\n")
        fh.write("end;\n")

    write_tsv(os.path.join(CORE, "core_gene_lengths.tsv"),
              ["cluster", "alignment_length_bp", "concat_start_1based", "concat_end_1based"],
              [dict(cluster=g, alignment_length_bp=L, concat_start_1based=a, concat_end_1based=b)
               for g, a, b, L in bounds])

    # per-genome QC
    qc = []
    for acc in roster:
        s = "".join(concat[acc])
        c = Counter(s)
        qc.append(dict(accession=acc, length=len(s), ACGT=c["A"] + c["C"] + c["G"] + c["T"],
                       N=c["N"], gap=c["-"],
                       ACGT_frac=f"{(c['A'] + c['C'] + c['G'] + c['T']) / len(s):.6f}",
                       missing_frac=f"{(c['N'] + c['-']) / len(s):.6f}"))
    write_tsv(os.path.join(CORE, "alignment_QC.tsv"),
              ["accession", "length", "ACGT", "N", "gap", "ACGT_frac", "missing_frac"], qc)

    # site summary
    n_const = n_var = n_pi = 0
    col_count = Counter()
    seqs_all = ["".join(concat[acc]) for acc in roster]
    for i in range(total_len):
        col = [s[i] for s in seqs_all]
        bases = Counter(ch for ch in col if ch in "ACGT")
        if len(bases) <= 1:
            n_const += 1
        else:
            n_var += 1
            if sum(1 for v in bases.values() if v >= 2) >= 2:
                n_pi += 1
    # identical sequence groups
    groups = defaultdict(list)
    for acc, s in zip(roster, seqs_all):
        groups[s].append(acc)
    ident = sorted((sorted(v) for v in groups.values() if len(v) > 1), key=lambda v: (-len(v), v))

    write_tsv(os.path.join(CORE, "site_summary.tsv"),
              ["metric", "value"],
              [dict(metric="n_genes", value=len(gene_order)),
               dict(metric="total_length_bp", value=total_len),
               dict(metric="constant_sites", value=n_const),
               dict(metric="variable_sites", value=n_var),
               dict(metric="parsimony_informative_sites", value=n_pi),
               dict(metric="n_identical_sequence_groups", value=len(ident)),
               dict(metric="invariant_sites_retained", value="yes"),
               dict(metric="trimming", value="none"),
               dict(metric="allgap_columns_removed", value=sum(len(v) for v in allgap.values()))])
    write_tsv(os.path.join(CORE, "identical_sequence_groups.tsv"),
              ["group_index", "n_members", "accessions"],
              [dict(group_index=i + 1, n_members=len(v), accessions=";".join(v))
               for i, v in enumerate(ident)])

    print(json.dumps(dict(status="PASS", n_genes=len(gene_order), total_length_bp=total_len,
                          constant=n_const, variable=n_var, parsimony_informative=n_pi,
                          identical_groups=len(ident)), indent=1))


if __name__ == "__main__":
    main()
