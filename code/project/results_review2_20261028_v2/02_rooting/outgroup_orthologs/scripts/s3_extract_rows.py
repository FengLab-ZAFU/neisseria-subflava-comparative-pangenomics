#!/usr/bin/env python
"""s3: hit processing, outgroup row construction, concatenation, report.

Pre-registered rules:
- best hit per gene x outgroup by bitscore; HSPs on same contig+strand are
  merged (chained by query coords, union coverage; group bitscore = sum of
  HSP bitscores -- documented implementation choice)
- accept if merged query coverage >= 50% of gap-stripped consensus length
  AND weighted identity >= 70% (weighted by HSP alignment length)
- uniqueness: a second independent hit (different contig, or same contig
  with subject interval not overlapping the best hit's) with bitscore
  >= 90% of the best hit's -> ambiguous, exclude
- row construction: walk merged HSP alignments; query base + subject base
  -> subject base at map[query_pos]; query base + subject gap -> '-';
  query gap + subject base -> discarded (counted); uncovered columns '-'
"""
import json, os, hashlib, subprocess, sys
from collections import defaultdict

ROOT = os.environ["NEISSERIA_PROJECT_ROOT"]
OUT = f"{ROOT}/results_review2_20261028_v2/02_rooting/outgroup_orthologs"
GENES_TSV = f"{ROOT}/task_unpack_20261028/remote_revision2_task_20261028_v2/inputs/phylogeny/core_gene_lengths.tsv"
OUTGROUPS = ["SAMN20865411", "SAMEA122247875", "SAMEA122247926",
             "SAMEA122247929", "SAMEA122247940"]
FULL_LEN = 906074

def union_len(intervals):
    if not intervals:
        return 0
    iv = sorted(intervals)
    total = 0; cs, ce = iv[0]
    for s, e in iv[1:]:
        if s <= ce + 1:
            ce = max(ce, e)
        else:
            total += ce - cs + 1; cs, ce = s, e
    total += ce - cs + 1
    return total

def union_positions(intervals):
    pos = set()
    for s, e in intervals:
        pos.update(range(s, e + 1))
    return pos

def main():
    genes = []
    with open(GENES_TSV) as fh:
        fh.readline()
        for line in fh:
            p = line.rstrip("\n").split("\t")
            genes.append((p[0], int(p[1]), int(p[2]), int(p[3])))
    assert len(genes) == 962
    assert sum(g[1] for g in genes) == FULL_LEN
    cmap = json.load(open(f"{OUT}/consensus_map.json"))
    cons_len = {g: len(cmap[g]) for g, *_ in genes}

    # ---- parse blast outputs, group HSPs ----
    avail = {}          # (gene, og) -> record dict
    rows = {}           # og -> {gene: row string}
    discarded_total = defaultdict(int)
    conflict_total = defaultdict(int)
    covered_cols = defaultdict(int)   # og -> covered alignment columns (present genes)
    for og in OUTGROUPS:
        hsps = defaultdict(list)
        with open(f"{OUT}/blast/{og}.blastn.tsv") as fh:
            for line in fh:
                p = line.rstrip("\n").split("\t")
                (qid, sid, pident, length, mismatch, gapopen, qstart, qend,
                 sstart, send, evalue, bitscore, qseq, sseq) = p
                strand = "+" if int(send) >= int(sstart) else "-"
                hsps[qid].append(dict(
                    sid=sid, pident=float(pident), length=int(length),
                    qstart=int(qstart), qend=int(qend),
                    sstart=int(sstart), send=int(send),
                    evalue=float(evalue), bitscore=float(bitscore),
                    qseq=qseq, sseq=sseq, strand=strand))
        rows[og] = {}
        for gene, alen, cs, ce in genes:
            L = alen
            hlist = hsps.get(gene, [])
            if not hlist:
                avail[(gene, og)] = dict(status="absent")
                continue
            # group by contig+strand
            groups = defaultdict(list)
            for h in hlist:
                groups[(h["sid"], h["strand"])].append(h)
            ginfo = []
            for (sid, strand), gh in groups.items():
                qs = [(min(h["qstart"], h["qend"]), max(h["qstart"], h["qend"])) for h in gh]
                ss = [(min(h["sstart"], h["send"]), max(h["sstart"], h["send"])) for h in gh]
                ginfo.append(dict(
                    sid=sid, strand=strand, hsps=gh,
                    bitscore=sum(h["bitscore"] for h in gh),
                    qcov=union_len(qs) / cons_len[gene],
                    qintervals=qs,
                    s_min=min(s for s, _ in ss), s_max=max(e for _, e in ss),
                    wident=sum(h["pident"] * h["length"] for h in gh) /
                           sum(h["length"] for h in gh),
                    n_hsps=len(gh)))
            ginfo.sort(key=lambda g: -g["bitscore"])
            best = ginfo[0]
            # uniqueness: second independent hit
            second = None
            for g in ginfo[1:]:
                independent = (g["sid"] != best["sid"] or
                               g["s_max"] < best["s_min"] or g["s_min"] > best["s_max"])
                if independent:
                    second = g
                    break
            accepted = best["qcov"] >= 0.50 and best["wident"] >= 70.0
            if not accepted:
                avail[(gene, og)] = dict(
                    status="low_coverage", identity=best["wident"],
                    coverage=best["qcov"], contig=best["sid"],
                    s_start=best["s_min"], s_end=best["s_max"],
                    n_hsps=best["n_hsps"], bitscore=best["bitscore"])
                continue
            if second is not None and second["bitscore"] >= 0.90 * best["bitscore"]:
                avail[(gene, og)] = dict(
                    status="ambiguous", identity=best["wident"],
                    coverage=best["qcov"], contig=best["sid"],
                    s_start=best["s_min"], s_end=best["s_max"],
                    n_hsps=best["n_hsps"], bitscore=best["bitscore"],
                    second_contig=second["sid"], second_bitscore=second["bitscore"],
                    second_s_start=second["s_min"], second_s_end=second["s_max"])
                continue
            # ---- build row ----
            row = ["-"] * L
            cmapg = cmap[gene]
            placed_cols = set()
            for h in sorted(best["hsps"], key=lambda h: h["qstart"]):
                qpos = h["qstart"]
                for qc, sc in zip(h["qseq"], h["sseq"]):
                    if qc != "-":
                        col = cmapg[qpos - 1]
                        placed_cols.add(col)
                        if sc != "-":
                            if row[col - 1] != "-" and row[col - 1] != sc:
                                conflict_total[og] += 1
                            elif row[col - 1] == "-":
                                row[col - 1] = sc
                        qpos += 1
                    else:
                        if sc != "-":
                            discarded_total[og] += 1
            assert qpos - 1 == h["qend"], (gene, og, qpos, h["qend"])
            rows[og][gene] = "".join(row)
            covered_cols[og] += len(placed_cols)
            avail[(gene, og)] = dict(
                status="present_unique", identity=best["wident"],
                coverage=best["qcov"], contig=best["sid"],
                s_start=best["s_min"], s_end=best["s_max"],
                n_hsps=best["n_hsps"], bitscore=best["bitscore"])

    # ---- stop rule: anomalous outgroup ----
    per_og = {}
    problems = []
    for og in OUTGROUPS:
        ids = [avail[(g, og)]["identity"] for g, *_ in genes
               if avail[(g, og)]["status"] == "present_unique"]
        npresent = len(ids)
        mean_id = sum(ids) / len(ids) if ids else 0.0
        per_og[og] = (npresent, mean_id)
        if npresent == 0 or mean_id < 60.0:
            problems.append((og, npresent, mean_id))
    if problems:
        print("STOP RULE TRIGGERED:", problems)
        sys.exit(2)

    # ---- per-gene outgroup fasta + availability matrix ----
    with open(f"{OUT}/availability_matrix.tsv", "w") as av:
        av.write("gene\toutgroup\tstatus\tidentity\tcoverage\tcontig\ts_start\ts_end"
                 "\tn_hsps\tbitscore\tsecond_contig\tsecond_bitscore\n")
        for gene, alen, cs, ce in genes:
            with open(f"{OUT}/per_gene/{gene}.outgroups.fasta", "w") as pf:
                for og in OUTGROUPS:
                    r = avail[(gene, og)]
                    av.write("\t".join(map(str, [
                        gene, og, r["status"],
                        f"{r.get('identity', 0):.2f}" if "identity" in r else "NA",
                        f"{r.get('coverage', 0):.4f}" if "coverage" in r else "NA",
                        r.get("contig", "NA"), r.get("s_start", "NA"), r.get("s_end", "NA"),
                        r.get("n_hsps", "NA"), r.get("bitscore", "NA"),
                        r.get("second_contig", "NA"), r.get("second_bitscore", "NA")])) + "\n")
                    if r["status"] == "present_unique":
                        s = rows[og][gene]
                        pf.write(f">{og}\n")
                        for i in range(0, len(s), 120):
                            pf.write(s[i:i+120] + "\n")

    # ---- concatenation into full 906074 columns ----
    full = {}
    for og in OUTGROUPS:
        buf = ["-"] * FULL_LEN
        for gene, alen, cs, ce in genes:
            if gene in rows[og]:
                s = rows[og][gene]
                assert len(s) == ce - cs + 1
                buf[cs-1:ce] = s
        full[og] = "".join(buf)
    with open(f"{OUT}/outgroup_rows_full906074.fasta", "w") as fh:
        for og in OUTGROUPS:
            assert len(full[og]) == FULL_LEN, (og, len(full[og]))
            fh.write(f">{og}\n")
            for i in range(0, FULL_LEN, 120):
                fh.write(full[og][i:i+120] + "\n")

    # ---- verification: 3 genes, concat slice == per-gene row ----
    check_genes = [genes[0][0], genes[481][0], genes[961][0]]
    verif = []
    for gene in check_genes:
        g = next(g for g in genes if g[0] == gene)
        for og in OUTGROUPS:
            slice_ = full[og][g[2]-1:g[3]]
            per = rows[og].get(gene, "-" * g[1])
            verif.append((gene, og, slice_ == per))
    assert all(v[2] for v in verif), [v for v in verif if not v[2]]

    # ---- stats + report ----
    report = {"parameters": {
        "blastn_task": "blastn", "word_size": 11, "dust": "no", "evalue": 1e-20,
        "outfmt": "6 std + qseq sseq (aligned strings needed for step-4 row construction)",
        "coverage_threshold": 0.50, "identity_threshold": 70.0,
        "ambiguity_bitscore_fraction": 0.90,
        "hit_grouping": "contig+strand; group bitscore = sum of HSP bitscores",
        "weighted_identity": "sum(pident*length)/sum(length) over merged HSPs",
        "consensus": "majority over A/C/G/T only; ties A>C>G>T; columns with no ACGT skipped",
        "queries_batched_per_outgroup": True,
        "note": "e-values are per query vs database; batching 962 queries in one "
                "run per outgroup does not change per-gene results"},
        "full_length": FULL_LEN, "n_genes": 962, "outgroups": {}}
    dist = defaultdict(int)
    for og in OUTGROUPS:
        stats = [avail[(g, og)]["status"] for g, *_ in genes]
        ids = [avail[(g, og)]["identity"] for g, *_ in genes
               if avail[(g, og)]["status"] == "present_unique"]
        amb = [g for g, *_ in genes if avail[(g, og)]["status"] == "ambiguous"]
        low = [g for g, *_ in genes if avail[(g, og)]["status"] == "low_coverage"]
        ab = [g for g, *_ in genes if avail[(g, og)]["status"] == "absent"]
        report["outgroups"][og] = {
            "present_unique": stats.count("present_unique"),
            "present_unique_frac": stats.count("present_unique") / 962,
            "ambiguous": len(amb), "low_coverage": len(low), "absent": len(ab),
            "covered_alignment_columns": covered_cols[og],
            "covered_fraction_of_906074": covered_cols[og] / FULL_LEN,
            "mean_identity": sum(ids) / len(ids),
            "min_identity": min(ids), "max_identity": max(ids),
            "discarded_bases": discarded_total[og],
            "placement_conflicts": conflict_total[og],
            "ambiguous_genes": amb, "low_coverage_genes": low, "absent_genes": ab}
    for g, *_ in genes:
        n = sum(1 for og in OUTGROUPS if avail[(g, og)]["status"] == "present_unique")
        dist[n] += 1
    report["outgroups_per_gene_distribution"] = {str(k): dist[k] for k in range(6)}
    report["verification_concat_vs_pergene"] = [
        {"gene": g, "outgroup": o, "match": m} for g, o, m in verif]

    with open(f"{OUT}/extraction_report.json", "w") as fh:
        json.dump(report, fh, indent=2)

    # sha256 of outputs
    def sha(p):
        h = hashlib.sha256()
        with open(p, "rb") as fh:
            for chunk in iter(lambda: fh.read(1 << 20), b""):
                h.update(chunk)
        return h.hexdigest()
    with open(f"{OUT}/checksums.sha256", "w") as cf:
        for f in ["consensus_all.fasta", "consensus_map.json", "availability_matrix.tsv",
                  "outgroup_rows_full906074.fasta", "extraction_report.json"]:
            cf.write(f"{sha(f'{OUT}/{f}')}  {f}\n")
        for og in OUTGROUPS:
            cf.write(f"{sha(f'{OUT}/blast/{og}.blastn.tsv')}  blast/{og}.blastn.tsv\n")

    # ---- console summary ----
    for og in OUTGROUPS:
        o = report["outgroups"][og]
        print(f"{og}: present={o['present_unique']}/962 ({o['present_unique_frac']:.3f}) "
              f"cov_cols={o['covered_alignment_columns']} ({o['covered_fraction_of_906074']:.3f}) "
              f"ident mean={o['mean_identity']:.2f} range=[{o['min_identity']:.2f},{o['max_identity']:.2f}] "
              f"amb={o['ambiguous']} low={o['low_coverage']} absent={o['absent']} "
              f"discarded={o['discarded_bases']} conflicts={o['placement_conflicts']}")
    print("outgroups_per_gene:", dict(sorted(dist.items())))
    print("verification:", all(v[2] for v in verif))

if __name__ == "__main__":
    main()
