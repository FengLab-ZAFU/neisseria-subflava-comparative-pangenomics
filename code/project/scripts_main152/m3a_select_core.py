#!/usr/bin/env python3
"""M3a: strict single-copy core selection + per-gene nucleotide FASTAs.

Selection rule: clusters of the main152 PGAP2-0.98
dictionary with copy==1 in ALL 152 genomes. Cluster order lexicographic by
cluster ID; sequences written in frozen roster order. DNA uppercase;
IUPAC ambiguity RYSWKMBDHV -> N with a per-gene/per-sample replacement ledger.
Sequences come from the M1-bound Bakta FFN keyed by (accession, gene_id).
"""
import csv
import json
import os
import sys
from multiprocessing import Pool

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import RESULTS, read_roster, read_tsv, write_tsv

PAN = os.path.join(RESULTS, "02_pangenome/pgap2_main152/pav_tables")
OUT = os.path.join(RESULTS, "03_phylogeny/strict_singlecopy_core")
ANN = os.path.join(RESULTS, "01_annotation")

IUPAC_AMB = set("RYSWKMBDHV")


def read_ffn(path):
    seqs, name, buf = {}, None, []
    with open(path) as fh:
        for ln in fh:
            if ln.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf).upper()
                name, buf = ln[1:].split()[0], []
            else:
                buf.append(ln.strip())
    if name is not None:
        seqs[name] = "".join(buf).upper()
    return seqs


def main():
    os.makedirs(os.path.join(OUT, "per_gene_fastas"), exist_ok=True)
    roster = read_roster()

    # copy matrix
    with open(os.path.join(PAN, "pav_copy.tsv"), newline="") as fh:
        rd = csv.reader(fh, delimiter="\t")
        hdr = next(rd)[1:]
        assert hdr == roster, "PAV column order != roster order"
        clusters = [(r[0], [int(x) for x in r[1:]]) for r in rd]

    singlecopy = sorted(cid for cid, cp in clusters if all(c == 1 for c in cp))
    sc_set = set(singlecopy)
    strict_multi = sorted(cid for cid, cp in clusters
                          if all(c >= 1 for c in cp) and any(c > 1 for c in cp))

    write_tsv(os.path.join(OUT, "core_gene_list.tsv"),
              ["cluster", "selection_rule"],
              [dict(cluster=c, selection_rule="copy==1 in all 152 genomes (new PGAP2-0.98 main152 dictionary)")
               for c in singlecopy])
    write_tsv(os.path.join(OUT, "excluded_gene_clusters.tsv"),
              ["cluster", "reason"],
              [dict(cluster=c, reason="Strict_core (present in all 152) but multi-copy in >=1 genome; "
                    "not eligible for single-copy core tree") for c in strict_multi])
    write_tsv(os.path.join(OUT, "cluster_selection_log.tsv"),
              ["metric", "value"],
              [dict(metric="total_clusters_main152_pgap2_orth98", value=len(clusters)),
               dict(metric="strict_core_clusters_freq_1", value=sum(1 for _, cp in clusters if all(c >= 1 for c in cp))),
               dict(metric="strict_core_multicopy_excluded", value=len(strict_multi)),
               dict(metric="singlecopy_core_selected", value=len(singlecopy)),
               dict(metric="selection_date_utc", value=__import__("datetime").datetime.utcnow().isoformat()),
               dict(metric="dictionary_id", value="pgap2_main152_orth0.98")])

    # occurrence map -> (cluster, accession) -> gene_id
    occ = {}
    for r in read_tsv(os.path.join(PAN, "occurrence_map.tsv")):
        if r["cluster"] in sc_set:
            key = (r["cluster"], r["accession"])
            if key in occ:
                sys.exit(f"FATAL: duplicate occurrence for {key}")
            occ[key] = r["gene_id"]

    # FFN paths from M1 manifest
    ffn = {r["accession"]: r["ffn_path_resolved"]
           for r in read_tsv(os.path.join(ANN, "annotation_manifest_main152.tsv"))}

    ledger = []
    missing = []
    per_gene_len = []
    for cid in singlecopy:
        out_fa = os.path.join(OUT, "per_gene_fastas", cid + ".fna")
        n_amb_total = 0
        with open(out_fa, "w") as out:
            for acc in roster:
                gid = occ.get((cid, acc))
                if gid is None:
                    missing.append((cid, acc, "no occurrence")); continue
                seqs = getattr(read_ffn, "_cache", {}).get(ffn[acc])
                if seqs is None:
                    seqs = read_ffn(ffn[acc])
                    read_ffn._cache = getattr(read_ffn, "_cache", {})
                    read_ffn._cache[ffn[acc]] = seqs
                s = seqs.get(gid)
                if s is None:
                    missing.append((cid, acc, f"gene_id {gid} not in FFN")); continue
                n_amb = sum(1 for ch in s if ch in IUPAC_AMB)
                if n_amb:
                    s = "".join("N" if ch in IUPAC_AMB else ch for ch in s)
                bad = set(s) - set("ACGTN")
                if bad:
                    missing.append((cid, acc, f"illegal chars {bad}")); continue
                if not s:
                    missing.append((cid, acc, "empty sequence")); continue
                n_amb_total += n_amb
                ledger.append(dict(cluster=cid, accession=acc, gene_id=gid,
                                   length_bp=len(s), iupac_amb_to_N=n_amb))
                out.write(f">{acc}\n")
                for i in range(0, len(s), 60):
                    out.write(s[i:i + 60] + "\n")
        per_gene_len.append((cid, os.path.getsize(out_fa)))

    write_tsv(os.path.join(OUT, "sequence_extraction_ledger.tsv"),
              ["cluster", "accession", "gene_id", "length_bp", "iupac_amb_to_N"], ledger)
    write_tsv(os.path.join(OUT, "per_gene_fasta_sizes.tsv"),
              ["cluster", "bytes"], [dict(cluster=c, bytes=b) for c, b in per_gene_len])

    summary = dict(module="M3a", n_singlecopy_clusters=len(singlecopy),
                   n_strict_core_multicopy_excluded=len(strict_multi),
                   n_sequences_written=len(ledger),
                   expected_sequences=len(singlecopy) * 152,
                   n_iupac_replaced=sum(int(r["iupac_amb_to_N"]) for r in ledger),
                   missing_or_invalid=missing,
                   status="PASS" if not missing and len(ledger) == len(singlecopy) * 152 else "FAILED")
    with open(os.path.join(OUT, "m3a_selection_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=1)
    print(json.dumps({k: v for k, v in summary.items() if k != "missing_or_invalid"}, indent=1))
    if missing:
        print("MISSING:", missing[:20])
        sys.exit(1)


if __name__ == "__main__":
    main()
