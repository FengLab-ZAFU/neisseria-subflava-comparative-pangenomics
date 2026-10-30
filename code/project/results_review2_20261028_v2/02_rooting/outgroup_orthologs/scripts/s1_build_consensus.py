#!/usr/bin/env python
"""s1: per-gene majority consensus of the 152 ingroup aligned sequences.

Rules (pre-registered):
- per alignment column count only A/C/G/T (ignore N and gap)
- majority wins; ties broken in order A>C>G>T
- columns with zero ACGT do not enter the consensus
- save map consensus_pos(1-based) -> alignment column(1-based)
Outputs: consensus_all.fasta, consensus_map.json, consensus_stats.tsv
"""
import json, sys, os
import numpy as np

ROOT = os.environ["NEISSERIA_PROJECT_ROOT"]
ALN_DIR = f"{ROOT}/results_review2_20261028_v2/01_recombination/per_gene_alignments"
GENES_TSV = f"{ROOT}/task_unpack_20261028/remote_revision2_task_20261028_v2/inputs/phylogeny/core_gene_lengths.tsv"
OUT = f"{ROOT}/results_review2_20261028_v2/02_rooting/outgroup_orthologs"

BASE2INT = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
INT2BASE = "ACGT"

def read_fasta(path):
    names, seqs = [], []
    name, buf = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    names.append(name); seqs.append("".join(buf))
                name = line[1:].split()[0]; buf = []
            else:
                buf.append(line.strip())
    if name is not None:
        names.append(name); seqs.append("".join(buf))
    return names, seqs

def main():
    genes = []
    with open(GENES_TSV) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            p = line.rstrip("\n").split("\t")
            genes.append((p[0], int(p[1]), int(p[2]), int(p[3])))
    assert len(genes) == 962, len(genes)

    cmap = {}            # gene -> list of alignment columns (1-based), index = consensus_pos-1
    stats = []
    with open(f"{OUT}/consensus_all.fasta", "w") as out:
        for gene, alen, cstart, cend in genes:
            names, seqs = read_fasta(f"{ALN_DIR}/{gene}.fasta")
            assert len(names) == 152, (gene, len(names))
            L = len(seqs[0])
            assert L == alen == cend - cstart + 1, (gene, L, alen, cstart, cend)
            assert all(len(s) == L for s in seqs), gene
            arr = np.zeros((152, L), dtype=np.int8)  # 0-3 ACGT, -1 other
            arr[:] = -1
            for i, s in enumerate(seqs):
                b = np.frombuffer(s.upper().encode("ascii"), dtype=np.uint8)
                for j, ch in enumerate(b"ACGT"):
                    arr[i, b == ch] = j
            counts = np.stack([(arr == j).sum(axis=0) for j in range(4)], axis=0)  # 4 x L
            total = counts.sum(axis=0)
            keep = total > 0
            # np.argmax returns first maximum -> tie order A>C>G>T
            cons = counts.argmax(axis=0)
            cols = (np.nonzero(keep)[0] + 1).tolist()  # 1-based alignment columns
            cmap[gene] = cols
            cseq = "".join(INT2BASE[cons[k]] for k in np.nonzero(keep)[0])
            out.write(f">{gene}\n")
            for i in range(0, len(cseq), 120):
                out.write(cseq[i:i+120] + "\n")
            stats.append((gene, L, len(cols), L - len(cols)))
    with open(f"{OUT}/consensus_map.json", "w") as fh:
        json.dump(cmap, fh)
    with open(f"{OUT}/logs/consensus_stats.tsv", "w") as fh:
        fh.write("gene\talignment_len\tconsensus_len\tskipped_columns\n")
        for g, L, cl, sk in stats:
            fh.write(f"{g}\t{L}\t{cl}\t{sk}\n")
    print(f"genes={len(genes)} total_consensus_bases={sum(s[2] for s in stats)} "
          f"skipped_columns={sum(s[3] for s in stats)}")

if __name__ == "__main__":
    main()
