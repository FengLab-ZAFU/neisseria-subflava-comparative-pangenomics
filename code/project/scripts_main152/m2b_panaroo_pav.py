#!/usr/bin/env python3
"""M2b: Panaroo strict PAV tables + three-way method sensitivity summary.

Panaroo 1.8.0 strict run on the same 152 staged GFF3 files. Cluster IDs are
NOT comparable with PGAP2 dictionaries; comparison is at count/frequency level
only (plan M2b).
"""
import csv
import json
import os
import sys
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import RESULTS, read_roster

OUT = os.path.join(RESULTS, "02_pangenome")
PAN = os.path.join(OUT, "panaroo_main152_strict")


def read_rtab(path):
    with open(path, newline="") as fh:
        rd = csv.reader(fh, delimiter="\t")
        hdr = next(rd)
        rows = [(r[0], r[1:]) for r in rd if r and r[0]]
    return hdr[1:], rows


def pan_judger(n_present, n):
    f = n_present / n
    if f == 1:
        return "Strict_core"
    if f >= 0.99:
        return "Core"
    if f >= 0.95:
        return "Soft_core"
    if f >= 0.15:
        return "Shell"
    if f > 0:
        return "Cloud"
    raise ValueError


def main():
    roster = read_roster()
    tbl = os.path.join(PAN, "pav_tables")
    os.makedirs(tbl, exist_ok=True)
    samples, rows = read_rtab(os.path.join(PAN, "gene_presence_absence.Rtab"))
    errors = []
    if set(samples) != set(roster):
        errors.append(f"sample set mismatch: only_rtab={sorted(set(samples)-set(roster))[:5]}")
    n = len(samples)
    order = [samples.index(a) for a in roster]

    bin_rows, freq_rows = [], []
    part_count = Counter()
    n_singletons = 0
    gene_totals = [0] * n
    for gid, vals in rows:
        binv = [1 if v not in ("", "0") else 0 for v in vals]
        np_ = sum(binv)
        if np_ == 0:
            errors.append(f"{gid}: all-zero row")
        f = np_ / n
        p = pan_judger(np_, n)
        part_count[p] += 1
        if np_ == 1:
            n_singletons += 1
        for i, b in enumerate(binv):
            gene_totals[i] += b
        bin_rows.append((gid, [binv[i] for i in order]))
        freq_rows.append((gid, np_, f"{f:.6f}", p))

    with open(os.path.join(tbl, "pav_binary.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["Gene"] + roster)
        for gid, v in bin_rows:
            w.writerow([gid] + v)
    with open(os.path.join(tbl, "cluster_frequency.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["cluster", "n_genomes_present", "frequency",
                    "partition_unified_frequency_rule"])
        w.writerows(freq_rows)
    with open(os.path.join(tbl, "checks.txt"), "w") as fh:
        fh.write(f"set=panaroo_main152_strict tool=panaroo clean_mode=strict remove_invalid_genes=yes\n")
        fh.write(f"n_samples={n} n_clusters={len(rows)} singletons={n_singletons}\n")
        fh.write(f"partitions_unified_rule={dict(part_count)}\n")
        fh.write(f"assertions={'ALL_PASS' if not errors else 'FAIL'}\n")
        for e in errors[:50]:
            fh.write(f"ERROR {e}\n")

    # three-way sensitivity summary (counts only; cluster IDs not joinable)
    summ = json.load(open(os.path.join(OUT, "m2a_pav_summary.json")))
    lines = [["metric", "pgap2_orth98", "pgap2_orth95", "panaroo_strict", "note"]]
    m98, m95 = summ["main152"], summ["main152_orth95"]
    pa = dict(n_samples=n, total_clusters=len(rows), strict_core=part_count["Strict_core"],
              core=part_count["Core"], soft_core=part_count["Soft_core"],
              shell=part_count["Shell"], cloud=part_count["Cloud"],
              singletons=n_singletons,
              mean_genes_per_genome=f"{sum(gene_totals) / n:.1f}")
    for k in ["n_samples", "total_clusters", "strict_core", "core", "soft_core",
              "shell", "cloud", "singletons", "mean_genes_per_genome"]:
        lines.append([k, m98[k], m95[k], pa[k],
                      "cluster IDs not joinable across dictionaries" if k != "n_samples" else ""])
    with open(os.path.join(OUT, "method_sensitivity_summary.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerows(lines)
    print(f"panaroo strict: clusters={pa['total_clusters']} strict_core={pa['strict_core']} "
          f"shell={pa['shell']} cloud={pa['cloud']} errors={len(errors)}")
    sys.exit(1 if errors else 0)


if __name__ == "__main__":
    main()
