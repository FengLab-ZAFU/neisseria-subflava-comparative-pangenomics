#!/usr/bin/env python3
"""Generate PERMANOVA label configurations for Module B.

Seeds and rules match analysis_contract.json; uses numpy.random.default_rng.
Writes label-configuration matrices and a partial permutation-design audit.
"""
from pathlib import Path
import csv
import math
import os
import sys
from collections import Counter, defaultdict

import numpy as np

P = str(Path(os.environ["NEISSERIA_PROJECT_ROOT"]) / 'boundary8_split_handoff_20261019')
S = str(Path(os.environ["NEISSERIA_PROJECT_ROOT"]) / 'results_revision_20261019_boundary8_analysis')
OUT = os.path.join(S, "05_structure")
os.makedirs(OUT, exist_ok=True)

N_PERM_TARGET = 999
SEED_FREE = 20261015990001
SEED_RESTRICTED = 20261015990002


def read_meta(path):
    with open(path, newline="", encoding="utf-8") as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        return list(rdr)


def multinomial_count(counts):
    n = sum(counts)
    num = math.factorial(n)
    for v in counts:
        num //= math.factorial(v)
    return num


def total_label_configs(labels, blocks):
    by_block = defaultdict(list)
    for i, b in enumerate(blocks):
        by_block[b].append(i)
    total = 1
    for b, idx in by_block.items():
        total *= multinomial_count(list(Counter(labels[i] for i in idx).values()))
    return total


def generate_label_configs(labels, blocks, n_want, seed):
    """Sample unique non-identity label configurations within blocks."""
    labels = list(labels)
    blocks = list(blocks)
    identity = tuple(labels)
    by_block = defaultdict(list)
    for i, b in enumerate(blocks):
        by_block[b].append(i)
    block_labels = {b: tuple(labels[i] for i in idx) for b, idx in by_block.items()}
    rng = np.random.default_rng(seed)
    seen = {identity}
    configs = []
    if total_label_configs(labels, blocks) - 1 <= n_want:
        # Exact enumeration (not expected for our scopes)
        import itertools
        per_block = {b: sorted(set(itertools.permutations(block_labels[b]))) for b in by_block}
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
        return configs
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
        if t in seen:
            continue
        seen.add(t)
        configs.append(t)
    return configs


def write_label_matrix(path, configs, accessions):
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["perm_id"] + accessions)
        for k, cfg in enumerate(configs, 1):
            w.writerow([k] + list(cfg))


def main():
    mo = read_meta(os.path.join(P, "00_cohort", "metadata_overall152.tsv"))
    mp = read_meta(os.path.join(P, "00_cohort", "metadata_comparison144.tsv"))
    acc_overall = [r["strain"] for r in mo]
    acc_primary = [r["strain"] for r in mp]

    mo_by_acc = {r["strain"]: r for r in mo}
    mp_by_acc = {r["strain"]: r for r in mp}

    labels_overall_sens = [mo_by_acc[a]["group_for_sensitivity152"] for a in acc_overall]
    blocks_overall = [mo_by_acc[a]["study_namespace"] for a in acc_overall]
    labels_primary = [mp_by_acc[a]["display_group"] for a in acc_primary]
    blocks_primary = [mp_by_acc[a]["study_namespace"] for a in acc_primary]

    scopes = {
        "primary144": {
            "acc": acc_primary,
            "labels": labels_primary,
            "blocks_free": ["__ALL__"] * 144,
            "blocks_restricted": blocks_primary,
            "n_sub": sum(1 for x in labels_primary if x == "subflava"),
            "n_per": sum(1 for x in labels_primary if x == "perflava"),
            "n_flav": sum(1 for x in labels_primary if x == "flavescens"),
        },
        "sensitivity152": {
            "acc": acc_overall,
            "labels": labels_overall_sens,
            "blocks_free": ["__ALL__"] * 152,
            "blocks_restricted": blocks_overall,
            "n_sub": sum(1 for x in labels_overall_sens if x == "subflava"),
            "n_per": sum(1 for x in labels_overall_sens if x == "perflava"),
            "n_flav": sum(1 for x in labels_overall_sens if x == "flavescens"),
        },
    }

    audit_rows = []
    for scope_name, sc in scopes.items():
        n = len(sc["acc"])
        # Free design
        total_free = multinomial_count([sc["n_sub"], sc["n_per"], sc["n_flav"]])
        configs_free = generate_label_configs(sc["labels"], sc["blocks_free"], N_PERM_TARGET, SEED_FREE)
        write_label_matrix(os.path.join(OUT, f"permutation_labels_{scope_name}_free.tsv"),
                           configs_free, sc["acc"])
        audit_rows.append({
            "permutation_design": "free", "block_id": "__ALL__", "block_size": n,
            "n_subflava": sc["n_sub"], "n_perflava": sc["n_per"], "n_flavescens": sc["n_flav"],
            "distinct_label_configurations_in_block": str(total_free),
            "total_label_space": "NA", "estimable": "NA", "enumeration_used": "NA",
            "B_used": "NA", "seed": "NA",
            "notes": "same-label swaps do not count as new configurations",
            "scope": scope_name,
        })
        audit_rows.append({
            "permutation_design": "free", "block_id": "__TOTAL__", "block_size": n,
            "n_subflava": sc["n_sub"], "n_perflava": sc["n_per"], "n_flavescens": sc["n_flav"],
            "distinct_label_configurations_in_block": "NA",
            "total_label_space": str(total_free), "estimable": "yes", "enumeration_used": "no",
            "B_used": len(configs_free), "seed": str(SEED_FREE),
            "notes": f"Monte Carlo: {len(configs_free)} unique non-identity label configurations sampled; upper-tail P=(1+b)/(B+1); identity configuration excluded from B",
            "scope": scope_name,
        })

        # Restricted design
        by_block = defaultdict(list)
        for i, b in enumerate(sc["blocks_restricted"]):
            by_block[b].append(i)
        for b in sorted(by_block):
            idx = by_block[b]
            cnt = Counter(sc["labels"][i] for i in idx)
            ncfg = multinomial_count(list(cnt.values()))
            audit_rows.append({
                "permutation_design": "restricted_within_study_namespace", "block_id": b,
                "block_size": len(idx),
                "n_subflava": cnt.get("subflava", 0), "n_perflava": cnt.get("perflava", 0),
                "n_flavescens": cnt.get("flavescens", 0),
                "distinct_label_configurations_in_block": str(ncfg),
                "total_label_space": "NA", "estimable": "NA", "enumeration_used": "NA",
                "B_used": "NA", "seed": "NA",
                "notes": "same-label swaps do not count as new configurations",
                "scope": scope_name,
            })
        total_restricted = total_label_configs(sc["labels"], sc["blocks_restricted"])
        configs_restricted = generate_label_configs(sc["labels"], sc["blocks_restricted"], N_PERM_TARGET, SEED_RESTRICTED)
        write_label_matrix(os.path.join(OUT, f"permutation_labels_{scope_name}_restricted.tsv"),
                           configs_restricted, sc["acc"])
        audit_rows.append({
            "permutation_design": "restricted_within_study_namespace", "block_id": "__TOTAL__",
            "block_size": n,
            "n_subflava": sc["n_sub"], "n_perflava": sc["n_per"], "n_flavescens": sc["n_flav"],
            "distinct_label_configurations_in_block": "NA",
            "total_label_space": str(total_restricted), "estimable": "yes", "enumeration_used": "no",
            "B_used": len(configs_restricted), "seed": str(SEED_RESTRICTED),
            "notes": f"Monte Carlo: {len(configs_restricted)} unique non-identity label configurations sampled; upper-tail P=(1+b)/(B+1); identity configuration excluded from B",
            "scope": scope_name,
        })

    with open(os.path.join(OUT, "permutation_design_audit_partial.tsv"), "w", newline="", encoding="utf-8") as fh:
        cols = ["permutation_design", "block_id", "block_size", "n_subflava", "n_perflava",
                "n_flavescens", "distinct_label_configurations_in_block", "total_label_space",
                "estimable", "enumeration_used", "B_used", "seed", "notes", "scope"]
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(audit_rows)

    print("PERMANOVA permutations generated.")


if __name__ == "__main__":
    main()
