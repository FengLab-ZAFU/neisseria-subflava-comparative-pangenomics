#!/usr/bin/env python3
"""M2: PAV / partition / occurrence tables for the main152 PGAP2 sets.

Adapted from the historically verified parser
(scripts/revision_20261010_phaseB_v7/p2_pav_tables.py) — same PGAP2 2.4.1
native outputs, same assertions. Partition by genome frequency (plan M2a):
Strict_core f==1; Core 0.99<=f<1; Soft_core 0.95<=f<0.99; Shell 0.15<=f<0.95;
Cloud 0<f<0.15 (n=152: 152 / 151 / 145-150 / 23-144 / 1-22).
"""
import csv
import os
import re
import sys
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import RESULTS, read_roster

OUT = os.path.join(RESULTS, "02_pangenome")


def sniff_delim(path):
    with open(path, newline="") as fh:
        first = fh.readline()
    return "\t" if "\t" in first else ","


def read_matrix(path):
    d = sniff_delim(path)
    with open(path, newline="") as fh:
        rd = csv.reader(fh, delimiter=d)
        hdr = next(rd)
        while hdr and hdr[-1] == "":
            hdr.pop()
        rows = []
        for row in rd:
            if not row or row[0] == "":
                continue
            while len(row) < len(hdr) + 1:
                row.append("")
            rows.append((row[0], row[1:len(hdr)]))
    return hdr[1:], rows


def pan_judger(n_present, n_total):
    f = n_present / n_total
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
    raise ValueError(f"invalid freq {n_present}/{n_total}")


def write_matrix(path, row_header, samples, rows):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow([row_header] + samples)
        for rid, vals in rows:
            w.writerow([rid] + list(vals))


def process_pgap2(set_name, orth_id, roster):
    d = os.path.join(OUT, f"pgap2_{set_name}")
    tbl = os.path.join(d, "pav_tables")
    os.makedirs(tbl, exist_ok=True)

    samples, copy_rows = read_matrix(os.path.join(d, "pgap2.partition.gene_content.pav"))
    errors = []
    if len(samples) != len(set(samples)):
        errors.append("duplicate columns in pav")
    if set(samples) != set(roster):
        errors.append(f"sample set != roster: only_pav={sorted(set(samples) - set(roster))[:5]} "
                      f"only_roster={sorted(set(roster) - set(samples))[:5]}")
    n = len(samples)

    detail = {}
    with open(os.path.join(d, "pgap2.partition.gene_content.detail.tsv"), newline="") as fh:
        rd = csv.reader(fh, delimiter="\t")
        next(rd)
        for row in rd:
            if len(row) < 13:
                continue
            detail[row[0]] = {"group": row[3], "involved_strain": int(row[9]),
                              "involved_gene": int(row[11])}

    freq_rows, part_rows = [], []
    copy_out, bin_out = [], []
    n_clust = 0
    part_count = Counter()
    n_singletons = 0
    gene_totals = [0] * n
    for cid, vals in copy_rows:
        copies = [int(v) if v not in ("", None) else 0 for v in vals]
        n_clust += 1
        if cid not in detail:
            errors.append(f"{cid} missing in detail.tsv")
            continue
        det = detail[cid]
        if sum(copies) != det["involved_gene"]:
            errors.append(f"{cid}: copy rowsum {sum(copies)} != involved_gene {det['involved_gene']}")
        binv = [1 if c > 0 else 0 for c in copies]
        n_present = sum(binv)
        if n_present == 0:
            errors.append(f"{cid}: all-zero row (not allowed)")
        if n_present != det["involved_strain"]:
            errors.append(f"{cid}: binary rowsum {n_present} != involved_strain {det['involved_strain']}")
        recomputed = pan_judger(n_present, n)
        if recomputed != det["group"]:
            errors.append(f"{cid}: recomputed {recomputed} != pgap2 group {det['group']}")
        part_count[det["group"]] += 1
        if n_present == 1:
            n_singletons += 1
        for i, c in enumerate(copies):
            gene_totals[i] += c
        copy_out.append((cid, copies))
        bin_out.append((cid, binv))
        freq_rows.append((cid, n_present, f"{n_present / n:.6f}", det["group"]))
        part_rows.append((cid, det["group"], recomputed,
                          "MATCH" if recomputed == det["group"] else "MISMATCH"))
    if sum(part_count.values()) != n_clust:
        errors.append(f"partition sum {sum(part_count.values())} != total clusters {n_clust}")

    # order columns by frozen roster order for all exported matrices
    order = [samples.index(a) for a in roster]
    write_matrix(os.path.join(tbl, "pav_copy.tsv"), "#Clust", roster,
                 [(cid, [v[i] for i in order]) for cid, v in copy_out])
    write_matrix(os.path.join(tbl, "pav_binary.tsv"), "#Clust", roster,
                 [(cid, [v[i] for i in order]) for cid, v in bin_out])
    with open(os.path.join(tbl, "cluster_frequency.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["cluster", "n_genomes_present", "frequency", "partition"])
        w.writerows(freq_rows)
    with open(os.path.join(tbl, "partition_table.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["cluster", "partition_pgap2", "partition_recomputed", "check"])
        w.writerows(part_rows)

    # per-genome gene burden (copy sums per sample, roster order)
    with open(os.path.join(tbl, "per_genome_gene_burden.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["accession", "n_CDS_copies", "n_clusters_present",
                    "n_accessory_clusters", "n_singleton_clusters"])
        for i, acc in enumerate(roster):
            copies_i = [v[order[i]] for _, v in copy_out]
            n_acc_clusters = 0
            w_acc = 0
            n_single = 0
            for (cid, binv), (cid2, grp, _, _) in zip(bin_out, part_rows):
                if binv[order[i]]:
                    n_acc_clusters += 1
                    if grp in ("Shell", "Cloud"):
                        w_acc += 1
                    if sum(binv) == 1:
                        n_single += 1
            w.writerow([acc, sum(copies_i), n_acc_clusters, w_acc, n_single])

    # occurrence map: gene content csv x involved_annot
    annot = {}
    with open(os.path.join(d, "total.involved_annot.tsv"), newline="") as fh:
        rd = csv.reader(fh, delimiter="\t")
        next(rd)
        for row in rd:
            if len(row) < 6:
                continue
            m = re.match(r"\[(\d+):(\d+)\]\((.)\)", row[3])
            annot[row[5]] = (row[1], row[2], m.group(1) if m else "", m.group(2) if m else "",
                             m.group(3) if m else "")
    g_samples, g_rows = read_matrix(os.path.join(d, "pgap2.partition.gene_content.csv"))
    n_occ = 0
    with open(os.path.join(tbl, "occurrence_map.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["cluster", "accession", "gene_id", "contig", "start", "end", "strand"])
        for cid, cells in g_rows:
            for sample, cell in zip(g_samples, cells):
                if not cell:
                    continue
                for gid in cell.split(";"):  # multi-copy cells join gene IDs with ';'
                    gid = gid.strip().strip('"')
                    if not gid:
                        continue
                    a = annot.get(gid)
                    if a is None:
                        errors.append(f"gene {gid} of {cid} missing in involved_annot")
                        continue
                    if a[0] != sample:
                        errors.append(f"gene {gid} strain {a[0]} != pav column {sample}")
                    w.writerow([cid, sample, gid, a[1], a[2], a[3], a[4]])
                    n_occ += 1

    with open(os.path.join(tbl, "checks.txt"), "w") as fh:
        fh.write(f"set={set_name} tool=pgap2 orth_id={orth_id}\n")
        fh.write(f"n_samples={n} n_clusters={n_clust} n_occurrences={n_occ}\n")
        fh.write(f"partitions={dict(part_count)}\n")
        fh.write(f"assertions={'ALL_PASS' if not errors else 'FAIL'}\n")
        for e in errors[:100]:
            fh.write(f"ERROR {e}\n")
    summary = {"set": set_name, "tool": "pgap2", "orth_id": orth_id,
               "n_samples": n, "total_clusters": n_clust,
               "strict_core": part_count["Strict_core"], "core": part_count["Core"],
               "soft_core": part_count["Soft_core"], "shell": part_count["Shell"],
               "cloud": part_count["Cloud"], "singletons": n_singletons,
               "n_occurrences": n_occ,
               "mean_genes_per_genome": f"{sum(gene_totals) / n:.1f}"}
    return summary, errors


def main():
    roster = read_roster()
    summaries = {}
    bad = 0
    for set_name, orth in [("main152", "0.98"), ("main152_orth95", "0.95")]:
        s, errs = process_pgap2(set_name, orth, roster)
        summaries[set_name] = s
        print(f"pgap2 {set_name}: {s}")
        if errs:
            bad += len(errs)
            for e in errs[:10]:
                print("  ERROR", e)
    import json
    with open(os.path.join(OUT, "m2a_pav_summary.json"), "w") as fh:
        json.dump(summaries, fh, indent=1)
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
