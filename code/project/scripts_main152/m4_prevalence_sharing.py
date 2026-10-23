#!/usr/bin/env python3
"""M4a: per-species cluster prevalence, seven-group sharing patterns,
group-conserved / group-limited clusters, per-genome burden by species group,
partition summary and the contract per-genome burden table.

Dictionary: pgap2_main152_orth0.98 (152 genomes; subflava 31 / perflava 114 /
flavescens 7). All counts are cluster-level detections on the binary PAV.
"Group-limited" = detected in >=1 genome of that group and 0 in both other
groups; it does NOT mean every genome of the group carries the cluster.
"""
import csv
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import RESULTS, read_roster, read_tsv, write_tsv, GROUPS, RUN_ID

DICT_ID = "pgap2_main152_orth0.98"
PAV_DIR = os.path.join(RESULTS, "02_pangenome/pgap2_main152/pav_tables")
OUT = os.path.join(RESULTS, "04_species_pangenome")

PARTITION_ORDER = ["Strict_core", "Core", "Soft_core", "Shell", "Cloud"]


def load_binary_pav():
    roster = read_roster()
    with open(os.path.join(PAV_DIR, "pav_binary.tsv"), newline="", encoding="utf-8") as fh:
        r = csv.reader(fh, delimiter="\t")
        header = next(r)
        assert header[0] == "#Clust", header[0]
        cols = header[1:]
        assert cols == roster, "PAV column order != roster order"
        clusters, rows = [], []
        for row in r:
            clusters.append(row[0])
            rows.append([int(x) for x in row[1:]])
    M = np.asarray(rows, dtype=np.uint8)
    assert M.shape == (13215, 152), M.shape
    assert set(np.unique(M)) <= {0, 1}
    assert (M.sum(axis=1) > 0).all(), "all-zero cluster row in dictionary"
    return clusters, M, roster


def main():
    clusters, M, roster = load_binary_pav()
    group_map = {r["strain"]: r["group_for_testing"] for r in
                 read_tsv(os.path.join(RESULTS, "00_freeze/input_packet_snapshot/00_input/cohort/main152_metadata.tsv"))}
    assert set(group_map) == set(roster)
    col = {g: np.array([i for i, a in enumerate(roster) if group_map[a] == g], dtype=int)
           for g in GROUPS}
    n_tot = {g: len(col[g]) for g in GROUPS}
    assert n_tot == {"subflava": 31, "perflava": 114, "flavescens": 7}

    part = {r["cluster"]: r["partition_recomputed"] for r in
            read_tsv(os.path.join(PAV_DIR, "partition_table.tsv"))}
    assert all(r["check"] == "MATCH" for r in
               read_tsv(os.path.join(PAV_DIR, "partition_table.tsv")))
    assert set(part) == set(clusters)

    pres = {g: M[:, col[g]].sum(axis=1) for g in GROUPS}  # per-cluster present counts

    # --- 1. cluster_species_prevalence.tsv (F4_FREQUENCIES) ---
    rows = []
    for i, c in enumerate(clusters):
        rows.append({
            "cluster": c,
            "subflava_present": int(pres["subflava"][i]), "subflava_total": 31,
            "perflava_present": int(pres["perflava"][i]), "perflava_total": 114,
            "flavescens_present": int(pres["flavescens"][i]), "flavescens_total": 7,
            "partition": part[c],
        })
    write_tsv(os.path.join(OUT, "cluster_species_prevalence.tsv"),
              ["cluster", "subflava_present", "subflava_total", "perflava_present",
               "perflava_total", "flavescens_present", "flavescens_total", "partition"], rows)
    assert len(rows) == 13215

    # --- 2. sharing_patterns.tsv (F4_SHARING): seven non-empty group-detection combos ---
    pattern_defs = [
        ("subflava_only", lambda s, p, f: (s > 0) & (p == 0) & (f == 0)),
        ("perflava_only", lambda s, p, f: (s == 0) & (p > 0) & (f == 0)),
        ("flavescens_only", lambda s, p, f: (s == 0) & (p == 0) & (f > 0)),
        ("subflava_perflava", lambda s, p, f: (s > 0) & (p > 0) & (f == 0)),
        ("subflava_flavescens", lambda s, p, f: (s > 0) & (p == 0) & (f > 0)),
        ("perflava_flavescens", lambda s, p, f: (s == 0) & (p > 0) & (f > 0)),
        ("all_three", lambda s, p, f: (s > 0) & (p > 0) & (f > 0)),
    ]
    s_, p_, f_ = pres["subflava"], pres["perflava"], pres["flavescens"]
    rows, total = [], 0
    for name, fn in pattern_defs:
        n = int(fn(s_, p_, f_).sum())
        rows.append({"pattern": name, "n_clusters": n, "dictionary_id": DICT_ID})
        total += n
    assert total == 13215, total
    write_tsv(os.path.join(OUT, "sharing_patterns.tsv"),
              ["pattern", "n_clusters", "dictionary_id"], rows)

    # --- 3. group_conserved.tsv ---
    # all_present: detected in every group genome (present == n).
    # ge95pct: 100*present >= 95*n on the real denominator (for n=7 this is 7/7).
    # group_limited: present >= 1 in this group AND 0 in both other groups.
    rows = []
    for g in GROUPS:
        n = n_tot[g]
        others = [o for o in GROUPS if o != g]
        all_present = int((pres[g] == n).sum())
        min_ge95 = -(-95 * n // 100)  # smallest integer with 100*p >= 95*n
        assert 100 * min_ge95 >= 95 * n and 100 * (min_ge95 - 1) < 95 * n
        ge95 = int((100 * pres[g] >= 95 * n).sum())
        limited = int(((pres[g] > 0) & (pres[others[0]] == 0) & (pres[others[1]] == 0)).sum())
        rows.append({
            "group": g, "n_genomes": n,
            "n_clusters_all_present": all_present,
            "all_present_definition": f"present in all {n}/{n} group genomes",
            "n_clusters_ge95pct": ge95,
            "ge95pct_min_present": min_ge95,
            "ge95pct_definition": f"present in >=95% of group genomes (100*present >= 95*{n}; minimum {min_ge95}/{n})",
            "n_clusters_group_limited": limited,
            "group_limited_definition": "present in >=1 genome of this group and 0 genomes in both other groups; NOT per-genome-universal and NOT a species-diagnostic marker",
            "dictionary_id": DICT_ID,
        })
    write_tsv(os.path.join(OUT, "group_conserved.tsv"),
              ["group", "n_genomes", "n_clusters_all_present", "all_present_definition",
               "n_clusters_ge95pct", "ge95pct_min_present", "ge95pct_definition",
               "n_clusters_group_limited", "group_limited_definition", "dictionary_id"], rows)

    # --- 4. per_genome_burden_by_group.tsv: four distinct counts x three groups ---
    burden = {r["accession"]: r for r in
              read_tsv(os.path.join(PAV_DIR, "per_genome_gene_burden.tsv"))}
    assert set(burden) == set(roster)
    metrics = [
        ("n_CDS_copies", "total CDS copies assigned to dictionary clusters (multi-copy counted)"),
        ("n_clusters", "distinct clusters detected (binary presence)"),
        ("n_accessory_clusters", "distinct Shell+Cloud clusters detected"),
        ("n_singletons", "distinct singleton clusters (present in exactly 1 of 152 genomes) detected"),
    ]
    burden_key = {"n_CDS_copies": "n_CDS_copies", "n_clusters": "n_clusters_present",
                  "n_accessory_clusters": "n_accessory_clusters",
                  "n_singletons": "n_singleton_clusters"}
    rows = []
    for metric, desc in metrics:
        for g in GROUPS:
            vals = np.array([int(burden[roster[i]][burden_key[metric]]) for i in col[g]], dtype=float)
            q1, q3 = np.percentile(vals, [25, 75])
            rows.append({
                "metric": metric, "metric_definition": desc, "group": g,
                "n": len(vals), "mean": f"{vals.mean():.6f}", "sd": f"{vals.std(ddof=1):.6f}",
                "median": f"{np.median(vals):.6f}", "q1": f"{q1:.6f}", "q3": f"{q3:.6f}",
                "iqr": f"{q3 - q1:.6f}", "min": int(vals.min()), "max": int(vals.max()),
            })
    write_tsv(os.path.join(OUT, "per_genome_burden_by_group.tsv"),
              ["metric", "metric_definition", "group", "n", "mean", "sd", "median",
               "q1", "q3", "iqr", "min", "max"], rows)

    # --- 5. partition_summary.tsv (F3_PARTITIONS) ---
    n_present = M.sum(axis=1)
    rows = []
    total = 0
    for pname in PARTITION_ORDER:
        idx = [i for i, c in enumerate(clusters) if part[c] == pname]
        n = len(idx)
        rows.append({"partition": pname, "n_clusters": n,
                     "n_singletons": int((n_present[idx] == 1).sum()),
                     "dictionary_id": DICT_ID})
        total += n
    assert total == 13215
    write_tsv(os.path.join(OUT, "partition_summary.tsv"),
              ["partition", "n_clusters", "n_singletons", "dictionary_id"], rows)

    # --- 6. per_genome_gene_burden.tsv (F3_GENOMES contract re-export with group) ---
    rows = []
    for a in roster:
        b = burden[a]
        rows.append({"accession": a, "species_group": group_map[a],
                     "n_CDS_copies": b["n_CDS_copies"], "n_clusters": b["n_clusters_present"],
                     "n_accessory_clusters": b["n_accessory_clusters"],
                     "n_singletons": b["n_singleton_clusters"]})
    write_tsv(os.path.join(OUT, "per_genome_gene_burden.tsv"),
              ["accession", "species_group", "n_CDS_copies", "n_clusters",
               "n_accessory_clusters", "n_singletons"], rows)

    print("M4a done. Sharing total = 13215 closed; rows written to", OUT)
    for g in GROUPS:
        print(g, "all_present:", int((pres[g] == n_tot[g]).sum()),
              "ge95:", int((100 * pres[g] >= 95 * n_tot[g]).sum()))


if __name__ == "__main__":
    main()
