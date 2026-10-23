#!/usr/bin/env python3
"""M4b: random-order accumulation curves (1000 iterations) on the binary PAV,
curve summaries, matched-sampling comparison summaries and n=7 effect sizes.

Groups: pooled152 (all 152), subflava (31), perflava (114), flavescens (7).
Per iteration a fresh random genome order is drawn from a deterministic seed:
seed = 20261015000000 + group_index*100000 + iteration (1..1000), so the random
stream never depends on execution/concurrency order. pan = union cluster count,
core = intersection cluster count of the first n genomes.

Endpoint verification (independent of the curve machinery): n=max pan/core
must equal the whole-group union/intersection recomputed from the matrix;
n=1 pan mean is Monte-Carlo (tolerance check vs per-genome cluster counts).
"""
import csv
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import RESULTS, read_roster, read_tsv, GROUPS

PAV_DIR = os.path.join(RESULTS, "02_pangenome/pgap2_main152/pav_tables")
OUT = os.path.join(RESULTS, "04_species_pangenome")
SEED_BASE = 20261015000000
N_ITER = 1000
ACC_GROUPS = ["pooled152", "subflava", "perflava", "flavescens"]  # group_index 0..3


def load_binary_pav():
    roster = read_roster()
    with open(os.path.join(PAV_DIR, "pav_binary.tsv"), newline="", encoding="utf-8") as fh:
        r = csv.reader(fh, delimiter="\t")
        header = next(r)
        cols = header[1:]
        assert cols == roster
        clusters, rows = [], []
        for row in r:
            clusters.append(row[0])
            rows.append([int(x) for x in row[1:]])
    return clusters, np.asarray(rows, dtype=bool), roster


def q(vals, p):
    return float(np.percentile(vals, p))


def summarize(mat):
    """mat: (n_iter, max_n) -> per-n stat dicts for pan and core."""
    out = []
    for j in range(mat.shape[1]):
        v = mat[:, j].astype(float)
        out.append({"mean": v.mean(), "median": float(np.median(v)),
                    "q025": q(v, 2.5), "q25": q(v, 25), "q75": q(v, 75),
                    "q975": q(v, 97.5)})
    return out


def main():
    clusters, M, roster = load_binary_pav()
    group_map = {r["strain"]: r["group_for_testing"] for r in
                 read_tsv(os.path.join(RESULTS, "00_freeze/input_packet_snapshot/00_input/cohort/main152_metadata.tsv"))}
    col = {"pooled152": np.arange(152)}
    for g in GROUPS:
        col[g] = np.array([i for i, a in enumerate(roster) if group_map[a] == g], dtype=int)

    os.makedirs(OUT, exist_ok=True)
    raw_path = os.path.join(OUT, "accumulation_raw.tsv")
    pan_store, core_store = {}, {}
    with open(raw_path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["group", "iteration", "seed", "n", "pan_clusters", "core_clusters"])
        for gi, g in enumerate(ACC_GROUPS):
            idx = col[g]
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

    # --- accumulation_summary.tsv (recomputed from raw stores) ---
    rows = []
    for g in ACC_GROUPS:
        ps, cs = summarize(pan_store[g]), summarize(core_store[g])
        for j in range(len(ps)):
            rows.append({"group": g, "n": j + 1,
                         "pan_mean": f"{ps[j]['mean']:.6f}", "pan_median": f"{ps[j]['median']:.1f}",
                         "pan_q025": f"{ps[j]['q025']:.3f}", "pan_q25": f"{ps[j]['q25']:.3f}",
                         "pan_q75": f"{ps[j]['q75']:.3f}", "pan_q975": f"{ps[j]['q975']:.3f}",
                         "core_mean": f"{cs[j]['mean']:.6f}", "core_median": f"{cs[j]['median']:.1f}",
                         "core_q025": f"{cs[j]['q025']:.3f}", "core_q25": f"{cs[j]['q25']:.3f}",
                         "core_q75": f"{cs[j]['q75']:.3f}", "core_q975": f"{cs[j]['q975']:.3f}"})
    with open(os.path.join(OUT, "accumulation_summary.tsv"), "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    # --- endpoint verification (independent recomputation from the matrix) ---
    per_genome_clusters = M.sum(axis=0)  # clusters per genome
    vrows = []
    for g in ACC_GROUPS:
        idx = col[g]
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
        n1_ok = abs(n1_mean - n1_expected) <= 5 * n1_se and \
            set(pan_store[g][:, 0].tolist()) <= set(pg.tolist()) and \
            set(core_store[g][:, 0].tolist()) <= set(pg.tolist())
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
    with open(os.path.join(OUT, "accumulation_verification.tsv"), "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(vrows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(vrows)
    assert all(r["pan_nmax_match"] and r["core_nmax_match"] and r["n1_within_5se_and_values_subset"]
               for r in vrows), vrows

    # --- matched_sampling_summary.tsv (F3_MATCHED) ---
    def matched_rows(comparison, groups, n_max):
        rr = []
        for g in groups:
            for metric, store in (("pan_clusters", pan_store[g]), ("core_clusters", core_store[g])):
                for n in range(1, n_max + 1):
                    v = store[:, n - 1].astype(float)
                    rr.append({"comparison": comparison, "group": g, "n": n, "metric": metric,
                               "mean": f"{v.mean():.6f}", "median": f"{np.median(v):.1f}",
                               "q025": f"{q(v, 2.5):.3f}", "q975": f"{q(v, 97.5):.3f}"})
        return rr

    rows = matched_rows("three_groups_n1_7", GROUPS, 7)
    rows += matched_rows("subflava_perflava_n1_31", ["subflava", "perflava"], 31)
    with open(os.path.join(OUT, "matched_sampling_summary.tsv"), "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    # --- matched_sampling_effects.tsv: descriptive pairwise effect sizes ---
    # n=7 across the three groups (the shared sampling range), plus
    # subflava-vs-perflava at n=31 (their shared range endpoint).
    def cliff_delta(a, b):
        gt = (a[:, None] > b[None, :]).sum()
        lt = (a[:, None] < b[None, :]).sum()
        return float((gt - lt) / (len(a) * len(b)))

    erows = []
    pairs = [("three_groups_n7", 7, GROUPS),
             ("subflava_perflava_n31", 31, ["subflava", "perflava"])]
    for comparison, n, groups in pairs:
        for metric, stores in (("pan_clusters", pan_store), ("core_clusters", core_store)):
            for i in range(len(groups)):
                for j in range(i + 1, len(groups)):
                    ga, gb = groups[i], groups[j]
                    va = stores[ga][:, n - 1].astype(float)
                    vb = stores[gb][:, n - 1].astype(float)
                    pooled_sd = float(np.sqrt((va.var(ddof=1) + vb.var(ddof=1)) / 2))
                    erows.append({
                        "comparison": comparison, "metric": metric,
                        "group_a": ga, "group_b": gb, "n": n,
                        "mean_a": f"{va.mean():.6f}", "mean_b": f"{vb.mean():.6f}",
                        "median_a": f"{np.median(va):.1f}", "median_b": f"{np.median(vb):.1f}",
                        "mean_difference_a_minus_b": f"{va.mean() - vb.mean():.6f}",
                        "median_difference_a_minus_b": f"{np.median(va) - np.median(vb):.1f}",
                        "pooled_sd": f"{pooled_sd:.6f}",
                        "cohens_d": f"{(va.mean() - vb.mean()) / pooled_sd:.6f}" if pooled_sd > 0 else "NA",
                        "cliffs_delta": f"{cliff_delta(va, vb):.6f}",
                        "note": "descriptive effect size over 1000 random-order iterations; not an independent-sample CI"})
    with open(os.path.join(OUT, "matched_sampling_effects.tsv"), "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(erows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(erows)

    print("M4b accumulation done:", raw_path)
    for r in vrows:
        print(r["group"], "pan_end", r["pan_nmax_matrix_union"], r["pan_nmax_match"],
              "core_end", r["core_nmax_matrix_intersection"], r["core_nmax_match"],
              "n1", r["n1_pan_mean_curve"], "vs", r["n1_pan_mean_expected_pergenome"])


if __name__ == "__main__":
    main()
