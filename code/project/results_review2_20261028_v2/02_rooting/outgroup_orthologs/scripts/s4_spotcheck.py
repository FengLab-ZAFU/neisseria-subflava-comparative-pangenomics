import os
#!/usr/bin/env python
"""s4: manual QC spot checks for 3 genes. For each gene x chosen outgroups:
- raw blastn hit lines
- availability decision
- accepted HSPs with q/s coords
- consensus vs generated row: first/last 200 alignment columns with rulers
- independent re-verification of base placement: re-walk accepted HSPs and
  check every placed base against subject FASTA directly
"""
import json, sys

ROOT = os.environ["NEISSERIA_PROJECT_ROOT"]
OUT = f"{ROOT}/results_review2_20261028_v2/02_rooting/outgroup_orthologs"
TASK = f"{ROOT}/task_unpack_20261028/remote_revision2_task_20261028_v2"
OUTGROUPS = ["SAMN20865411", "SAMEA122247875", "SAMEA122247926",
             "SAMEA122247929", "SAMEA122247940"]

def read_fasta_dict(path):
    d = {}
    name = None
    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith(">"):
            name = line[1:].split()[0]; d[name] = []
        else:
            d[name].append(line.strip())
    return {k: "".join(v) for k, v in d.items()}

def revcomp(s):
    return s.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]

def main(gene_og_pairs):
    cmap = json.load(open(f"{OUT}/consensus_map.json"))
    cons = read_fasta_dict(f"{OUT}/consensus_all.fasta")
    per_gene = {}
    for gene, _ in gene_og_pairs:
        per_gene[gene] = read_fasta_dict(f"{OUT}/per_gene/{gene}.outgroups.fasta")
    og_seqs = {og: read_fasta_dict(f"{TASK}/inputs/outgroups/{og}.fna") for og in OUTGROUPS}

    for gene, og in gene_og_pairs:
        lines = []
        w = lines.append
        w(f"########## SPOT CHECK gene={gene} outgroup={og} ##########")
        raw = [l for l in open(f"{OUT}/blast/{og}.blastn.tsv")
               if l.startswith(gene + "\t")]
        w(f"--- raw blastn hits ({len(raw)} lines) ---")
        for l in raw:
            p = l.rstrip("\n").split("\t")
            w("  " + "\t".join(p[:12]))
        # decision from availability matrix
        for l in open(f"{OUT}/availability_matrix.tsv"):
            if l.startswith(gene + "\t" + og + "\t"):
                w("--- availability: " + l.rstrip("\n"))
        if gene not in per_gene or og not in per_gene.get(gene, {}):
            w("--- no row generated for this pair ---")
        else:
            row = per_gene[gene][og]
            cseq = cons[gene]; cols = cmap[gene]
            # rebuild consensus-position row in alignment columns for display
            w(f"--- row length={len(row)} consensus_len={len(cseq)} ---")
            # verify every placed base against subject FASTA
            n_checked = n_bad = 0
            for l in raw:
                p = l.rstrip("\n").split("\t")
                sid = p[1]; qs, qe = int(p[6]), int(p[7])
                ss, se = int(p[8]), int(p[9])
                strand = "+" if se >= ss else "-"
                subj = og_seqs[og][sid]
                qseq, sseq = p[12], p[13]
                spos = ss if strand == "+" else ss  # walk from sstart
                qpos = qs
                for qc, sc in zip(qseq, sseq):
                    if qc != "-":
                        col = cols[qpos - 1]
                        # subject base in genome orientation
                        if sc != "-":
                            genome_base = subj[spos - 1] if strand == "+" else revcomp(subj)[len(subj) - spos]
                            n_checked += 1
                            if genome_base.upper() != sc.upper():
                                n_bad += 1
                            if row[col - 1] != sc:
                                w(f"  MISMATCH row[col {col}]={row[col-1]} hsp_base={sc} qpos={qpos}")
                        qpos += 1
                        spos += 1 if strand == "+" else -1
                    else:
                        if sc != "-":
                            spos += 1 if strand == "+" else -1
            w(f"--- base-level verification vs subject FASTA: {n_checked} checked, {n_bad} mismatches vs genome")
            # display first/last 200 columns
            for label, lo, hi in [("FIRST", 0, 200), ("LAST", max(0, len(row) - 200), len(row))]:
                w(f"--- {label} 200 alignment columns ({lo+1}-{hi}) ---")
                aln_cons = ["-"] * len(row)
                for i, c in enumerate(cols):
                    aln_cons[c - 1] = cseq[i]
                aln_cons = "".join(aln_cons)
                w("col1k : " + "".join(str((i + 1) // 100 % 10) for i in range(lo, hi)))
                w("col10 : " + "".join(str((i + 1) // 10 % 10) for i in range(lo, hi)))
                w("col1  : " + "".join(str((i + 1) % 10) for i in range(lo, hi)))
                w("cons  : " + aln_cons[lo:hi])
                w(f"{og[:9]}: " + row[lo:hi])
        out = f"{OUT}/logs/spotcheck_{gene}_{og}.txt"
        with open(out, "w") as fh:
            fh.write("\n".join(lines) + "\n")
        print(out)

if __name__ == "__main__":
    pairs = [tuple(a.split(":")) for a in sys.argv[1:]]
    main(pairs)
