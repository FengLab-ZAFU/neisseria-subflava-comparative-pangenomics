#!/usr/bin/env python3
"""Module A: pangenome frequencies, sharing, accumulation for 144+8 grouping."""
import csv, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import P, R, S, read_tsv, write_tsv, write_json

DICT_ID = "pgap2_main152_orth0.98"
PAV_DIR = os.path.join(R, "02_pangenome/pgap2_main152/pav_tables")
OUT = os.path.join(S, "04_species_pangenome")
SEED_BASE = 20261015000000
N_ITER = 1000
GROUP_ORDER = ["pooled152", "subflava", "perflava", "flavescens"]

def load_binary_pav():
    with open(os.path.join(PAV_DIR, "pav_binary.tsv"), newline='', encoding='utf-8') as fh:
        r = csv.reader(fh, delimiter='\t')
        header = next(r)
        assert header[0] == "#Clust"
        cols = header[1:]
        clusters, rows = [], []
        for row in r:
            clusters.append(row[0])
            rows.append([int(x) for x in row[1:]])
    M = np.asarray(rows, dtype=np.uint8)
    assert M.shape == (13215, 152), M.shape
    return clusters, M, cols

def q(vals, p):
    return float(np.percentile(vals, p))

def summarize(mat):
    out = []
    for j in range(mat.shape[1]):
        v = mat[:, j].astype(float)
        out.append({"mean": v.mean(), "median": float(np.median(v)),
                    "q025": q(v, 2.5), "q25": q(v, 25), "q75": q(v, 75), "q975": q(v, 97.5)})
    return out

def cliff_delta(a, b):
    gt = (a[:, None] > b[None, :]).sum()
    lt = (a[:, None] < b[None, :]).sum()
    return float((gt - lt) / (len(a) * len(b)))

def main():
    os.makedirs(OUT, exist_ok=True)
    clusters, M, cols = load_binary_pav()
    meta_overall = {r["strain"]: r for r in read_tsv(os.path.join(P, "00_cohort/metadata_overall152.tsv"))}
    meta_comp = {r["strain"]: r for r in read_tsv(os.path.join(P, "00_cohort/metadata_comparison144.tsv"))}
    bound = set(r.strip() for r in open(os.path.join(P, "00_cohort/boundary8_accessions.txt")) if r.strip())
    comp = set(meta_comp.keys())
    assert len(bound) == 8 and len(comp) == 144 and len(bound & comp) == 0

    # group column indices
    col_primary = {}
    for g in ["subflava", "perflava", "flavescens"]:
        col_primary[g] = np.array([i for i, a in enumerate(cols) if meta_comp.get(a, {}).get("group_for_testing") == g], dtype=int)
    col_sens = {}
    for g in ["subflava", "perflava", "flavescens"]:
        col_sens[g] = np.array([i for i, a in enumerate(cols) if meta_overall[a].get("group_for_sensitivity152") == g], dtype=int)
    col_all = np.arange(152)

    n_tot_primary = {g: len(col_primary[g]) for g in col_primary}
    assert n_tot_primary == {"subflava": 23, "perflava": 114, "flavescens": 7}

    part = {r["cluster"]: r["partition_recomputed"] for r in read_tsv(os.path.join(PAV_DIR, "partition_table.tsv"))}
    assert set(part) == set(clusters)

    pres_primary = {g: M[:, col_primary[g]].sum(axis=1) for g in col_primary}
    pres_sens = {g: M[:, col_sens[g]].sum(axis=1) for g in col_sens}
    pres_boundary = M[:, [i for i, a in enumerate(cols) if a in bound]].sum(axis=1)

    # 1. cluster_species_prevalence
    rows = []
    for i, c in enumerate(clusters):
        rows.append({
            "cluster": c,
            "subflava_present": int(pres_primary["subflava"][i]), "subflava_total": 23,
            "perflava_present": int(pres_primary["perflava"][i]), "perflava_total": 114,
            "flavescens_present": int(pres_primary["flavescens"][i]), "flavescens_total": 7,
            "boundary8_present": int(pres_boundary[i]), "boundary8_total": 8,
            "partition": part[c],
            "dictionary_id": DICT_ID,
        })
    write_tsv(os.path.join(OUT, "cluster_species_prevalence.tsv"), rows,
              ["cluster", "subflava_present", "subflava_total", "perflava_present", "perflava_total",
               "flavescens_present", "flavescens_total", "boundary8_present", "boundary8_total", "partition", "dictionary_id"])

    # 2. sharing_patterns on 144
    s_, p_, f_ = pres_primary["subflava"], pres_primary["perflava"], pres_primary["flavescens"]
    pattern_defs = [
        ("subflava_only", lambda s, p, f: (s > 0) & (p == 0) & (f == 0)),
        ("perflava_only", lambda s, p, f: (s == 0) & (p > 0) & (f == 0)),
        ("flavescens_only", lambda s, p, f: (s == 0) & (p == 0) & (f > 0)),
        ("subflava_perflava", lambda s, p, f: (s > 0) & (p > 0) & (f == 0)),
        ("subflava_flavescens", lambda s, p, f: (s > 0) & (p == 0) & (f > 0)),
        ("perflava_flavescens", lambda s, p, f: (s == 0) & (p > 0) & (f > 0)),
        ("all_three", lambda s, p, f: (s > 0) & (p > 0) & (f > 0)),
    ]
    sp_rows, total = [], 0
    for name, fn in pattern_defs:
        n = int(fn(s_, p_, f_).sum())
        sp_rows.append({"pattern": name, "n_clusters": n, "dictionary_id": DICT_ID})
        total += n
    write_tsv(os.path.join(OUT, "sharing_patterns.tsv"), sp_rows, ["pattern", "n_clusters", "dictionary_id"])

    # boundary8_only
    in144 = ((s_ > 0) | (p_ > 0) | (f_ > 0))
    only8 = (~in144) & (pres_boundary > 0)
    write_tsv(os.path.join(OUT, "boundary8_only_clusters.tsv"),
              [{"n_clusters": int(only8.sum()), "n_total_clusters": 13215,
                "definition": "present in at least one boundary8 genome and absent from all 144 comparison genomes",
                "dictionary_id": DICT_ID}],
              ["n_clusters", "n_total_clusters", "definition", "dictionary_id"])
    assert total + int(only8.sum()) == 13215

    # 3. group_conserved
    gc_rows = []
    for g in ["subflava", "perflava", "flavescens"]:
        n = n_tot_primary[g]
        others = [o for o in ["subflava", "perflava", "flavescens"] if o != g]
        all_present = int((pres_primary[g] == n).sum())
        min_ge95 = -(-95 * n // 100)
        ge95 = int((100 * pres_primary[g] >= 95 * n).sum())
        limited = int(((pres_primary[g] > 0) & (pres_primary[others[0]] == 0) & (pres_primary[others[1]] == 0)).sum())
        gc_rows.append({
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
    write_tsv(os.path.join(OUT, "group_conserved.tsv"), gc_rows,
              ["group", "n_genomes", "n_clusters_all_present", "all_present_definition",
               "n_clusters_ge95pct", "ge95pct_min_present", "ge95pct_definition",
               "n_clusters_group_limited", "group_limited_definition", "dictionary_id"])

    # 4. per_genome_burden_by_group
    burden = {r["accession"]: r for r in read_tsv(os.path.join(PAV_DIR, "per_genome_gene_burden.tsv"))}
    metrics = [
        ("n_CDS_copies", "total CDS copies assigned to dictionary clusters (multi-copy counted)"),
        ("n_clusters", "distinct clusters detected (binary presence)"),
        ("n_accessory_clusters", "distinct Shell+Cloud clusters detected"),
        ("n_singleton_clusters", "distinct singleton clusters detected"),
    ]
    burden_key = {"n_CDS_copies": "n_CDS_copies", "n_clusters": "n_clusters_present",
                  "n_accessory_clusters": "n_accessory_clusters", "n_singleton_clusters": "n_singleton_clusters"}
    bg_rows = []
    scopes = [("primary144", col_primary), ("sensitivity152", col_sens), ("descriptive8", {"boundary8": np.array([i for i, a in enumerate(cols) if a in bound], dtype=int)})]
    for scope_name, colmap in scopes:
        if scope_name == "descriptive8":
            groups = [("boundary8", colmap["boundary8"])]
        else:
            groups = [(g, colmap[g]) for g in ["subflava", "perflava", "flavescens"]]
        for metric, desc in metrics:
            for g, idx in groups:
                vals = np.array([int(burden[cols[i]][burden_key[metric]]) for i in idx], dtype=float)
                q1, q3 = np.percentile(vals, [25, 75])
                bg_rows.append({
                    "analysis_scope": scope_name, "metric": metric, "metric_definition": desc, "group": g,
                    "n": len(vals), "mean": f"{vals.mean():.6f}", "sd": f"{vals.std(ddof=1):.6f}",
                    "median": f"{np.median(vals):.6f}", "q1": f"{q1:.6f}", "q3": f"{q3:.6f}",
                    "iqr": f"{q3-q1:.6f}", "min": int(vals.min()), "max": int(vals.max()),
                })
    write_tsv(os.path.join(OUT, "per_genome_burden_by_group.tsv"), bg_rows,
              ["analysis_scope", "metric", "metric_definition", "group", "n", "mean", "sd", "median",
               "q1", "q3", "iqr", "min", "max"])

    # 5. Accumulation raw + summary + verification
    col_acc = {
        "pooled152": col_all,
        "subflava": col_primary["subflava"],
        "perflava": col_primary["perflava"],
        "flavescens": col_primary["flavescens"],
    }
    pan_store, core_store = {}, {}
    raw_path = os.path.join(OUT, "accumulation_raw.tsv")
    with open(raw_path, "w", newline='', encoding='utf-8') as fh:
        w = csv.writer(fh, delimiter='\t', lineterminator='\n')
        w.writerow(["group", "iteration", "seed", "n", "pan_clusters", "core_clusters"])
        for gi, g in enumerate(GROUP_ORDER):
            idx = col_acc[g]
            n_g = len(idx)
            pan = np.empty((N_ITER, n_g), dtype=np.int32)
            core = np.empty((N_ITER, n_g), dtype=np.int32)
            for it in range(1, N_ITER + 1):
                seed = SEED_BASE + gi * 100000 + it
                rng = np.random.default_rng(seed)
                sub = M[:, idx[rng.permutation(n_g)]]
                pan[it - 1] = np.logical_or.accumulate(sub, axis=1).sum(axis=0)
                core[it - 1] = np.logical_and.accumulate(sub, axis=1).sum(axis=0)
                for n in range(1, n_g + 1):
                    w.writerow([g, it, seed, n, int(pan[it - 1, n - 1]), int(core[it - 1, n - 1])])
            pan_store[g], core_store[g] = pan, core

    # summary
    sum_rows = []
    for g in GROUP_ORDER:
        ps, cs = summarize(pan_store[g]), summarize(core_store[g])
        for j in range(len(ps)):
            sum_rows.append({"group": g, "n": j + 1,
                             "pan_mean": f"{ps[j]['mean']:.6f}", "pan_median": f"{ps[j]['median']:.1f}",
                             "pan_q025": f"{ps[j]['q025']:.3f}", "pan_q25": f"{ps[j]['q25']:.3f}",
                             "pan_q75": f"{ps[j]['q75']:.3f}", "pan_q975": f"{ps[j]['q975']:.3f}",
                             "core_mean": f"{cs[j]['mean']:.6f}", "core_median": f"{cs[j]['median']:.1f}",
                             "core_q025": f"{cs[j]['q025']:.3f}", "core_q25": f"{cs[j]['q25']:.3f}",
                             "core_q75": f"{cs[j]['q75']:.3f}", "core_q975": f"{cs[j]['q975']:.3f}"})
    write_tsv(os.path.join(OUT, "accumulation_summary.tsv"), sum_rows,
              ["group", "n", "pan_mean", "pan_median", "pan_q025", "pan_q25", "pan_q75", "pan_q975",
               "core_mean", "core_median", "core_q025", "core_q25", "core_q75", "core_q975"])

    per_genome_clusters = M.sum(axis=0)
    vrows = []
    for g in GROUP_ORDER:
        idx = col_acc[g]
        sub = M[:, idx]
        union = int(sub.any(axis=1).sum())
        inter = int(sub.all(axis=1).sum())
        n_g = len(idx)
        pan_end_ok = bool(set(pan_store[g][:, -1].tolist()) == {union})
        core_end_ok = bool(set(core_store[g][:, -1].tolist()) == {inter})
        pg = per_genome_clusters[idx]
        n1_mean = float(pan_store[g][:, 0].mean())
        n1_expected = float(pg.mean())
        n1_se = float(pg.std(ddof=1) / np.sqrt(N_ITER))
        n1_ok = abs(n1_mean - n1_expected) <= 5 * n1_se and set(pan_store[g][:, 0].tolist()) <= set(pg.tolist()) and set(core_store[g][:, 0].tolist()) <= set(pg.tolist())
        vrows.append({"group": g, "n_max": n_g,
                      "pan_nmax_matrix_union": union,
                      "pan_nmax_curve_unique_values": ",".join(map(str, sorted(set(pan_store[g][:, -1].tolist())))),
                      "pan_nmax_match": pan_end_ok,
                      "core_nmax_matrix_intersection": inter,
                      "core_nmax_curve_unique_values": ",".join(map(str, sorted(set(core_store[g][:, -1].tolist())))),
                      "core_nmax_match": core_end_ok,
                      "n1_pan_mean_curve": f"{n1_mean:.6f}",
                      "n1_pan_mean_expected_pergenome": f"{n1_expected:.6f}",
                      "n1_mc_se": f"{n1_se:.6f}",
                      "n1_within_5se_and_values_subset": n1_ok})
    write_tsv(os.path.join(OUT, "accumulation_verification.tsv"), vrows,
              ["group", "n_max", "pan_nmax_matrix_union", "pan_nmax_curve_unique_values", "pan_nmax_match",
               "core_nmax_matrix_intersection", "core_nmax_curve_unique_values", "core_nmax_match",
               "n1_pan_mean_curve", "n1_pan_mean_expected_pergenome", "n1_mc_se", "n1_within_5se_and_values_subset"])

    # matched sampling summary
    def matched_rows(comparison, groups, n_max, pan_store_loc, core_store_loc):
        rr = []
        for g in groups:
            for metric, store in (("pan_clusters", pan_store_loc[g]), ("core_clusters", core_store_loc[g])):
                for n in range(1, n_max + 1):
                    v = store[:, n - 1].astype(float)
                    rr.append({"comparison": comparison, "group": g, "n": n, "metric": metric,
                               "mean": f"{v.mean():.6f}", "median": f"{np.median(v):.1f}",
                               "q025": f"{q(v, 2.5):.3f}", "q975": f"{q(v, 97.5):.3f}"})
        return rr

    ms_rows = matched_rows("three_groups_n1_7", ["subflava", "perflava", "flavescens"], 7, pan_store, core_store)
    ms_rows += matched_rows("subflava_perflava_n1_23", ["subflava", "perflava"], 23, pan_store, core_store)
    # sensitivity 152 n=7/23 using sensitivity group columns
    pan_sens, core_sens = {}, {}
    for g in ["subflava", "perflava", "flavescens"]:
        idx = col_sens[g]
        n_g = len(idx)
        pan = np.empty((N_ITER, n_g), dtype=np.int32)
        core = np.empty((N_ITER, n_g), dtype=np.int32)
        for it in range(1, N_ITER + 1):
            seed = SEED_BASE + 1000000 + (list(col_sens.keys()).index(g) + 1) * 100000 + it
            rng = np.random.default_rng(seed)
            sub = M[:, idx[rng.permutation(n_g)]]
            pan[it - 1] = np.logical_or.accumulate(sub, axis=1).sum(axis=0)
            core[it - 1] = np.logical_and.accumulate(sub, axis=1).sum(axis=0)
        pan_sens[g], core_sens[g] = pan, core
    ms_rows += matched_rows("sensitivity152_three_groups_n1_7", ["subflava", "perflava", "flavescens"], 7, pan_sens, core_sens)
    ms_rows += matched_rows("sensitivity152_subflava_perflava_n1_23", ["subflava", "perflava"], 23, pan_sens, core_sens)
    write_tsv(os.path.join(OUT, "matched_sampling_summary.tsv"), ms_rows,
              ["comparison", "group", "n", "metric", "mean", "median", "q025", "q975"])

    # matched sampling effects
    erows = []
    def effect_rows(comparison, n, groups, pan_loc, core_loc):
        rows = []
        for metric, stores in (("pan_clusters", pan_loc), ("core_clusters", core_loc)):
            for i in range(len(groups)):
                for j in range(i + 1, len(groups)):
                    ga, gb = groups[i], groups[j]
                    va = stores[ga][:, n - 1].astype(float)
                    vb = stores[gb][:, n - 1].astype(float)
                    pooled_sd = float(np.sqrt((va.var(ddof=1) + vb.var(ddof=1)) / 2))
                    rows.append({
                        "comparison": comparison, "metric": metric, "group_a": ga, "group_b": gb, "n": n,
                        "mean_a": f"{va.mean():.6f}", "mean_b": f"{vb.mean():.6f}",
                        "median_a": f"{np.median(va):.1f}", "median_b": f"{np.median(vb):.1f}",
                        "mean_difference_a_minus_b": f"{va.mean() - vb.mean():.6f}",
                        "median_difference_a_minus_b": f"{np.median(va) - np.median(vb):.1f}",
                        "pooled_sd": f"{pooled_sd:.6f}",
                        "cohens_d": f"{(va.mean() - vb.mean()) / pooled_sd:.6f}" if pooled_sd > 0 else "NA",
                        "cliffs_delta": f"{cliff_delta(va, vb):.6f}",
                        "note": "descriptive effect size over 1000 random-order iterations; not an independent-sample CI"})
        return rows

    erows += effect_rows("three_groups_n7", 7, ["subflava", "perflava", "flavescens"], pan_store, core_store)
    erows += effect_rows("subflava_perflava_n23", 23, ["subflava", "perflava"], pan_store, core_store)
    erows += effect_rows("sensitivity152_three_groups_n7", 7, ["subflava", "perflava", "flavescens"], pan_sens, core_sens)
    erows += effect_rows("sensitivity152_subflava_perflava_n23", 23, ["subflava", "perflava"], pan_sens, core_sens)
    write_tsv(os.path.join(OUT, "matched_sampling_effects.tsv"), erows,
              ["comparison", "metric", "group_a", "group_b", "n", "mean_a", "mean_b",
               "median_a", "median_b", "mean_difference_a_minus_b", "median_difference_a_minus_b",
               "pooled_sd", "cohens_d", "cliffs_delta", "note"])

    write_json(os.path.join(OUT, "run_info.json"), {
        "status": "DONE", "dictionary_id": DICT_ID, "n_clusters": 13215,
        "n_iter": N_ITER, "seed_base": SEED_BASE,
        "primary_groups": n_tot_primary,
        "sensitivity_groups": {g: len(col_sens[g]) for g in col_sens},
        "sharing_total_144": total, "boundary8_only": int(only8.sum()),
        "notes": "All primary144 and sensitivity152 accumulation curves recomputed deterministically."
    })
    print("Module A done.")
    print("n_tot_primary", n_tot_primary)
    print("sharing_total_144", total, "boundary8_only", int(only8.sum()))

if __name__ == "__main__":
    main()
