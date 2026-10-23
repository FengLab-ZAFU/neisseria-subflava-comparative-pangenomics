#!/usr/bin/env python3
"""M4c (part 1): accessory feature sets, Jaccard distance matrices, PCoA,
and permutation-design audit / label permutations for the R statistics step.

Feature sets (dictionary pgap2_main152_orth0.98, binary PAV):
  shell_cloud   = clusters in partitions Shell or Cloud   (main)
  all_variable  = clusters present in 1..151 of 152 genomes (sensitivity)

Jaccard distance = 1 - |A&B|/|A|B||; both-zero pair -> 0, one-zero pair -> 1
(none expected; rule recorded either way). PCoA = classical MDS on D:
uncorrected eigenvalues are stored first; scale = max(abs(eigenvalues));
Lingoes correction only when min(eigenvalue) < -1e-8*scale and scale > 0,
otherwise tolerance-level negatives are treated as 0. Fraction denominators
are the sums of the reported positive eigenvalues; coordinates are written
for every positive axis, with a deterministic sign convention (the largest
absolute element of each axis is made positive).

Permutation designs (for m4_distance_stats.R):
  free                              = all 152 labels exchangeable
  restricted_within_study_namespace = labels shuffled only within blocks
Distinct label configurations count same-label swaps as identical. If the
label space were < 999 non-identity configurations we would enumerate it
exactly; here both spaces are astronomically large, so 999 unique non-identity
configurations are sampled per design (deterministic seeds).
"""
import csv
import math
import os
import sys
from collections import Counter, defaultdict

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import RESULTS, read_roster, read_tsv, write_tsv, GROUPS

PAV_DIR = os.path.join(RESULTS, "02_pangenome/pgap2_main152/pav_tables")
META = os.path.join(RESULTS, "00_freeze/input_packet_snapshot/00_input/cohort/main152_metadata.tsv")
OUT = os.path.join(RESULTS, "05_structure")
DICT_ID = "pgap2_main152_orth0.98"
TOL = 1e-8
N_PERM = 999
SEED_FREE = 20261015990001
SEED_RESTRICTED = 20261015990002


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
    return clusters, np.asarray(rows, dtype=np.uint8), roster


def write_matrix(path, clusters, X, roster):
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["#Clust"] + roster)
        for c, row in zip(clusters, X):
            w.writerow([c] + row.tolist())


def jaccard_distance(X):
    """X: (n_features, n_samples) binary -> (n, n) Jaccard distance."""
    Xf = X.astype(np.float64)
    inter = Xf.T @ Xf
    a = Xf.sum(axis=0)
    union = a[:, None] + a[None, :] - inter
    with np.errstate(invalid="ignore", divide="ignore"):
        D = 1.0 - inter / union
    both_zero = union == 0.0
    D[both_zero] = 0.0  # two all-zero profiles: distance 0 by convention
    np.fill_diagonal(D, 0.0)
    return D, both_zero


def pcoa(D):
    """Classical MDS with the plan's Lingoes rule. Returns dict of results."""
    n = D.shape[0]
    J = np.eye(n) - np.ones((n, n)) / n

    def eig_of(Dm):
        B = -0.5 * J @ (Dm ** 2) @ J
        B = (B + B.T) / 2.0
        w, V = np.linalg.eigh(B)
        return w[::-1], V[:, ::-1]

    w_unc, _ = eig_of(D)
    scale = float(np.max(np.abs(w_unc)))
    use_lingoes = scale > 0 and float(w_unc.min()) < -TOL * scale
    if use_lingoes:
        c = -float(w_unc.min())
        D2 = D ** 2 + c
        np.fill_diagonal(D2, 0.0)
        Dm = np.sqrt(D2)
        w_rep, V = eig_of(Dm)
        constant = c
    else:
        w_rep, V = w_unc, eig_of(D)[1]
        constant = 0.0
    pos = w_rep > 0.0
    coords = V[:, pos] * np.sqrt(w_rep[pos])
    # deterministic sign: largest-|element| of each axis positive
    for j in range(coords.shape[1]):
        k = int(np.argmax(np.abs(coords[:, j])))
        if coords[k, j] < 0:
            coords[:, j] *= -1.0
    denom = float(w_rep[pos].sum())
    return {"eigenvalues_uncorrected": w_unc, "eigenvalues_reported": w_rep,
            "correction": "lingoes" if use_lingoes else "none",
            "lingoes_constant": constant, "scale": scale,
            "positive_mask": pos, "coords": coords, "positive_sum": denom}


def write_distance_matrix(path, D, roster):
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["accession"] + roster)
        for a, row in zip(roster, D):
            w.writerow([a] + [f"{v:.12g}" for v in row])


def multinomial_count(counts):
    """Exact number of distinct label configurations (big int)."""
    n = sum(counts)
    num = math.factorial(n)
    for v in counts:
        num //= math.factorial(v)
    return num


def sample_configs(labels, blocks, n_want, seed):
    """Sample unique non-identity label configurations.

    labels: list of group labels (roster order); blocks: list of block ids.
    Returns (list_of_tuples, complete_space_size_as_big_int).
    """
    identity = tuple(labels)
    by_block = defaultdict(list)
    for i, b in enumerate(blocks):
        by_block[b].append(i)
    block_labels = {b: tuple(labels[i] for i in idx) for b, idx in by_block.items()}
    total = 1
    for b, idx in by_block.items():
        total *= multinomial_count(list(Counter(block_labels[b]).values()))
    rng = np.random.default_rng(seed)
    seen = {identity}
    configs = []
    if total - 1 <= n_want:  # exact enumeration of all non-identity configs
        from itertools import permutations as iperm
        def enum_block(bl):
            uniq = sorted(set(iperm(bl)))
            return uniq
        import itertools
        per_block = {b: enum_block(block_labels[b]) for b in by_block}
        blocks_order = list(by_block)
        for combo in itertools.product(*(per_block[b] for b in blocks_order)):
            cfg = [None] * len(labels)
            for b, lab in zip(blocks_order, combo):
                for pos, i in enumerate(by_block[b]):
                    cfg[i] = lab[pos]
            t = tuple(cfg)
            if t not in seen:
                seen.add(t)
                configs.append(t)
        return configs, total
    attempts = 0
    while len(configs) < n_want and attempts < n_want * 50:
        attempts += 1
        cfg = list(labels)
        for b, idx in by_block.items():
            perm = rng.permutation(len(idx))
            vals = [cfg[idx[p]] for p in perm]
            for pos, i in enumerate(idx):
                cfg[i] = vals[pos]
        t = tuple(cfg)
        if t not in seen:
            seen.add(t)
            configs.append(t)
    return configs, total


def main():
    os.makedirs(OUT, exist_ok=True)
    clusters, M, roster = load_binary_pav()
    meta = {r["strain"]: r for r in read_tsv(META)}
    group_of = {a: meta[a]["group_for_testing"] for a in roster}
    ns_of = {a: meta[a]["study_namespace"] for a in roster}
    part = {r["cluster"]: r["partition_recomputed"] for r in
            read_tsv(os.path.join(PAV_DIR, "partition_table.tsv"))}
    n_present = M.sum(axis=1)

    feat_defs = {
        "shell_cloud": [i for i, c in enumerate(clusters) if part[c] in ("Shell", "Cloud")],
        "all_variable": [i for i in range(len(clusters)) if 1 <= n_present[i] <= 151],
    }

    # --- permutation design audit + label configurations ---
    labels = [group_of[a] for a in roster]
    blocks_free = ["__ALL__"] * 152
    blocks_res = [ns_of[a] for a in roster]
    designs = {
        "free": (blocks_free, SEED_FREE),
        "restricted_within_study_namespace": (blocks_res, SEED_RESTRICTED),
    }
    audit_rows = []
    spec_rows = []
    for design, (blocks, seed) in designs.items():
        by_block = defaultdict(list)
        for i, b in enumerate(blocks):
            by_block[b].append(i)
        for b in sorted(by_block, key=lambda x: (x != "__ALL__", x)):
            idx = by_block[b]
            cnt = Counter(labels[i] for i in idx)
            ncfg = multinomial_count(list(cnt.values()))
            audit_rows.append({
                "permutation_design": design, "block_id": b, "block_size": len(idx),
                "n_subflava": cnt.get("subflava", 0), "n_perflava": cnt.get("perflava", 0),
                "n_flavescens": cnt.get("flavescens", 0),
                "distinct_label_configurations_in_block": str(ncfg),
                "total_label_space": "NA", "estimable": "NA", "enumeration_used": "NA",
                "B_used": "NA", "seed": "NA",
                "notes": "same-label swaps do not count as new configurations"})
        total = 1
        for b, idx in by_block.items():
            total *= multinomial_count(list(Counter(labels[i] for i in idx).values()))
        configs, total_check = sample_configs(labels, blocks, N_PERM, seed)
        assert total_check == total
        estimable = (total - 1) >= 1
        enumeration = (total - 1) <= N_PERM
        if not estimable:
            audit_rows.append({
                "permutation_design": design, "block_id": "__TOTAL__", "block_size": 152,
                "n_subflava": 31, "n_perflava": 114, "n_flavescens": 7,
                "distinct_label_configurations_in_block": "NA",
                "total_label_space": str(total), "estimable": "no",
                "enumeration_used": "no", "B_used": 0, "seed": seed,
                "notes": "no non-identity label rearrangement exists -> NOT_ESTIMABLE"})
            spec_rows.append({"permutation_design": design, "estimable": "no",
                              "n_configs": 0, "complete_enumeration": "no", "seed": seed})
            continue
        assert len(configs) == min(N_PERM, total - 1)
        lab_path = os.path.join(OUT, f"permutation_labels_{'free' if design == 'free' else 'restricted'}.tsv")
        with open(lab_path, "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh, delimiter="\t", lineterminator="\n")
            w.writerow(["perm_id"] + roster)
            for k, cfg in enumerate(configs, 1):
                w.writerow([k] + list(cfg))
        audit_rows.append({
            "permutation_design": design, "block_id": "__TOTAL__", "block_size": 152,
            "n_subflava": 31, "n_perflava": 114, "n_flavescens": 7,
            "distinct_label_configurations_in_block": "NA",
            "total_label_space": str(total), "estimable": "yes",
            "enumeration_used": "yes" if enumeration else "no",
            "B_used": len(configs), "seed": seed,
            "notes": ("label space smaller than 999 -> exact enumeration of all non-identity "
                      "configurations; exact tail probability uses the full label space as denominator")
            if enumeration else
            ("Monte Carlo: 999 unique non-identity label configurations sampled; "
             "upper-tail P=(1+b)/(B+1); identity configuration excluded from B")})
        spec_rows.append({"permutation_design": design, "estimable": "yes",
                          "n_configs": len(configs),
                          "complete_enumeration": "yes" if enumeration else "no",
                          "seed": seed,
                          "labels_file": os.path.basename(lab_path),
                          "total_label_space": str(total)})
    write_tsv(os.path.join(OUT, "permutation_design_audit.tsv"),
              ["permutation_design", "block_id", "block_size", "n_subflava", "n_perflava",
               "n_flavescens", "distinct_label_configurations_in_block", "total_label_space",
               "estimable", "enumeration_used", "B_used", "seed", "notes"], audit_rows)
    write_tsv(os.path.join(OUT, "permutation_spec.tsv"),
              ["permutation_design", "estimable", "n_configs", "complete_enumeration",
               "seed", "labels_file", "total_label_space"],
              [{k: r.get(k, "") for k in
                ["permutation_design", "estimable", "n_configs", "complete_enumeration",
                 "seed", "labels_file", "total_label_space"]} for r in spec_rows])

    # --- feature sets, matrices, Jaccard, PCoA ---
    import hashlib
    manifest_rows = []
    eig_rows = []
    coord_rows = []
    axes_rows = []
    pcoa_spec = []
    for fs, idx in feat_defs.items():
        fc = [clusters[i] for i in idx]
        X = M[idx, :]
        feat_dir = os.path.join(OUT, "feature_sets")
        os.makedirs(feat_dir, exist_ok=True)
        with open(os.path.join(feat_dir, f"feature_clusters_{fs}.tsv"), "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh, delimiter="\t", lineterminator="\n")
            w.writerow(["feature_set", "cluster", "partition", "n_genomes_present"])
            for i in idx:
                w.writerow([fs, clusters[i], part[clusters[i]], int(n_present[i])])
        mat_path = os.path.join(feat_dir, f"feature_matrix_{fs}.tsv")
        write_matrix(mat_path, fc, X, roster)
        mat_sha = hashlib.sha256(open(mat_path, "rb").read()).hexdigest()

        D, both_zero = jaccard_distance(X)
        zero_samples = [roster[i] for i in range(152) if X[:, i].sum() == 0]
        assert D.shape == (152, 152)
        assert np.allclose(D, D.T, atol=1e-12) and np.all(np.diag(D) == 0.0)
        dm_path = os.path.join(OUT, f"jaccard_distance_{fs}.tsv")
        write_distance_matrix(dm_path, D, roster)
        dm_sha = hashlib.sha256(open(dm_path, "rb").read()).hexdigest()

        pc = pcoa(D)
        n_pos = int(pc["positive_mask"].sum())
        w_unc, w_rep = pc["eigenvalues_uncorrected"], pc["eigenvalues_reported"]
        for ax in range(152):
            eig_rows.append({
                "feature_set": fs, "axis": ax + 1,
                "eigenvalue_uncorrected": f"{w_unc[ax]:.12g}",
                "eigenvalue_reported": f"{w_rep[ax]:.12g}",
                "positive_in_reported": bool(pc["positive_mask"][ax]),
                "correction": pc["correction"],
                "included_in_fraction_denominator": bool(pc["positive_mask"][ax]),
                "fraction_of_positive_sum": f"{w_rep[ax] / pc['positive_sum']:.12g}" if pc["positive_mask"][ax] else "0",
            })
        fr1 = w_rep[0] / pc["positive_sum"] if n_pos >= 1 else float("nan")
        fr2 = w_rep[1] / pc["positive_sum"] if n_pos >= 2 else float("nan")
        status = "OK" if n_pos >= 2 else ("DEGENERATE_ONE_POSITIVE_AXIS" if n_pos == 1 else "DEGENERATE_NO_POSITIVE_AXIS")
        for i, a in enumerate(roster):
            coord_rows.append({
                "feature_set": fs, "accession": a, "species_group": group_of[a],
                "PC1": f"{pc['coords'][i, 0]:.10g}" if n_pos >= 1 else "NA",
                "PC2": f"{pc['coords'][i, 1]:.10g}" if n_pos >= 2 else "NA",
                "axis1_fraction": f"{fr1:.8f}" if n_pos >= 1 else "NA",
                "axis2_fraction": f"{fr2:.8f}" if n_pos >= 2 else "NA",
                "correction": pc["correction"],
                "status": status,
            })
            for j in range(n_pos):
                axes_rows.append({"feature_set": fs, "accession": a, "axis": j + 1,
                                  "coordinate": f"{pc['coords'][i, j]:.10g}"})
        pcoa_spec.append({"feature_set": fs, "correction": pc["correction"],
                          "lingoes_constant": f"{pc['lingoes_constant']:.12g}",
                          "scale": f"{pc['scale']:.12g}",
                          "n_positive_axes": n_pos,
                          "axis1_fraction": f"{fr1:.8f}" if n_pos >= 1 else "NA",
                          "axis2_fraction": f"{fr2:.8f}" if n_pos >= 2 else "NA",
                          "fraction_denominator": "sum of all reported positive eigenvalues",
                          "fraction_basis": "corrected" if pc["correction"] == "lingoes" else "uncorrected",
                          "min_eigenvalue_uncorrected": f"{w_unc.min():.12g}",
                          "status": status})
        manifest_rows.append({
            "feature_set": fs, "dictionary_id": DICT_ID,
            "definition": ("clusters in partitions Shell or Cloud" if fs == "shell_cloud"
                           else "clusters present in 1..151 of 152 genomes (all non-strict-core variable)"),
            "n_features": len(idx),
            "feature_list_path": f"feature_sets/feature_clusters_{fs}.tsv",
            "matrix_path": f"feature_sets/feature_matrix_{fs}.tsv",
            "matrix_sha256": mat_sha,
            "distance_matrix_path": f"jaccard_distance_{fs}.tsv",
            "distance_matrix_sha256": dm_sha,
            "n_all_zero_samples": len(zero_samples),
            "all_zero_samples": ",".join(zero_samples) if zero_samples else "none",
            "jaccard_zero_rule": "intersection/union; both-zero pair -> distance 0, one-zero pair -> distance 1 (no all-zero sample present)" ,
            "symmetric": True, "diagonal_zero": True,
        })

    write_tsv(os.path.join(OUT, "eigenvalues.tsv"),
              ["feature_set", "axis", "eigenvalue_uncorrected", "eigenvalue_reported",
               "positive_in_reported", "correction", "included_in_fraction_denominator",
               "fraction_of_positive_sum"], eig_rows)
    write_tsv(os.path.join(OUT, "pcoa_coordinates.tsv"),
              ["feature_set", "accession", "species_group", "PC1", "PC2",
               "axis1_fraction", "axis2_fraction", "correction", "status"], coord_rows)
    write_tsv(os.path.join(OUT, "pcoa_positive_axes.tsv"),
              ["feature_set", "accession", "axis", "coordinate"], axes_rows)
    write_tsv(os.path.join(OUT, "pcoa_summary.tsv"),
              ["feature_set", "correction", "lingoes_constant", "scale", "n_positive_axes",
               "axis1_fraction", "axis2_fraction", "fraction_denominator", "fraction_basis",
               "min_eigenvalue_uncorrected", "status"], pcoa_spec)
    write_tsv(os.path.join(OUT, "feature_set_manifest.tsv"),
              ["feature_set", "dictionary_id", "definition", "n_features", "feature_list_path",
               "matrix_path", "matrix_sha256", "distance_matrix_path", "distance_matrix_sha256",
               "n_all_zero_samples", "all_zero_samples", "jaccard_zero_rule",
               "symmetric", "diagonal_zero"], manifest_rows)

    for r in pcoa_spec:
        print(r["feature_set"], "correction:", r["correction"], "pos_axes:", r["n_positive_axes"],
              "axis1:", r["axis1_fraction"], "axis2:", r["axis2_fraction"],
              "min_eig:", r["min_eigenvalue_uncorrected"])
    print("M4c part1 done ->", OUT)


if __name__ == "__main__":
    main()
