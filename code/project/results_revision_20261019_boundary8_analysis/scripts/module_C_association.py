#!/usr/bin/env python3
"""Module C: gene association + COG enrichment for boundary8_split_handoff_20261019.

Inputs:
  P/00_cohort/analysis_contract.json
  P/00_cohort/metadata_overall152.tsv
  P/00_cohort/metadata_comparison144.tsv
  R/02_pangenome/pgap2_main152/pav_tables/pav_binary.tsv
  R/02_pangenome/pgap2_main152/pav_tables/partition_table.tsv
  R/06_function/cluster_cog_map.tsv

Outputs (S/06_function/):
  cluster_filter_ledger.tsv
  association_all_tests.tsv
  associated_candidates_ranked.tsv
  COG_enrichment_all.tsv
  cog_foreground_background.tsv
  run_info.json
"""

import csv
import json
import math
import os
import sys
from collections import defaultdict
from pathlib import Path

import scipy.stats as st
from scipy.optimize import brentq, minimize_scalar

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
ROOT = Path(os.environ["NEISSERIA_PROJECT_ROOT"])
P = ROOT / "boundary8_split_handoff_20261019"
R = ROOT / "results_revision_20261015_main152_species"
S = ROOT / "results_revision_20261019_boundary8_analysis"

OUT_DIR = S / "06_function"
OUT_DIR.mkdir(parents=True, exist_ok=True)

CONTRACT = P / "00_cohort/analysis_contract.json"
META_OVERALL = P / "00_cohort/metadata_overall152.tsv"
META_COMPARISON = P / "00_cohort/metadata_comparison144.tsv"
PAV_BINARY = R / "02_pangenome/pgap2_main152/pav_tables/pav_binary.tsv"
PARTITION_TABLE = R / "02_pangenome/pgap2_main152/pav_tables/partition_table.tsv"
CLUSTER_COG_MAP = R / "06_function/cluster_cog_map.tsv"

GROUPS = ["subflava", "perflava", "flavescens"]
TARGET_GROUP_SIZES_PRIMARY = {"subflava": 23, "perflava": 114, "flavescens": 7}

# ---------------------------------------------------------------------------
# Formatting helpers (mimic R format(x, digits=15, scientific=FALSE))
# ---------------------------------------------------------------------------
def fmt_r(x):
    if x is None or (isinstance(x, float) and math.isnan(x)):
        return "NA"
    if isinstance(x, float) and math.isinf(x):
        return "Inf" if x > 0 else "-Inf"
    if x == 0:
        return "0"
    if not isinstance(x, (int, float)):
        x = float(x)
    exp = math.floor(math.log10(abs(x)))
    decimals = max(0, 15 - exp - 1)
    return ("{:." + str(decimals) + "f}").format(x)

# ---------------------------------------------------------------------------
# Fisher exact test: two-sided p via scipy; conditional MLE OR + exact CI via
# non-central hypergeometric distribution (scipy.stats.nchypergeom_fisher).
# ---------------------------------------------------------------------------
def fisher_test(a, b, c, d):
    """Return dict with OR, CI_low, CI_high, p, status, reason.

    Matches R fisher.test(two.sided) conditional MLE OR and exact CI.
    """
    total = a + b + c + d
    k = a + c          # population successes
    n = a + b          # sample size

    # No variation within evaluable genomes
    if b == 0 and d == 0:
        return {"OR": "NA", "CI_low": "NA", "CI_high": "NA",
                "p": "NA", "status": "NOT_ESTIMABLE",
                "reason": "no_variation_all_evaluable_carriers"}

    # p-value (scipy two-sided equals R two-sided for 2x2)
    pval = st.fisher_exact([[a, b], [c, d]], alternative="two-sided").pvalue

    # Conditional MLE OR
    if a == 0 or d == 0:
        or_est = 0.0
    elif b == 0 or c == 0:
        or_est = float("inf")
    else:
        def nll(logpsi):
            return -st.nchypergeom_fisher.logpmf(a, M=total, n=n,
                                                  odds=math.exp(logpsi), N=k)
        res = minimize_scalar(nll, bounds=(-25, 35), method="bounded")
        or_est = math.exp(res.x)

    # Exact conditional 95% CI by tail inversion
    alpha2 = 0.025
    if a == 0 or d == 0:
        ci_low = 0.0
    else:
        try:
            ci_low = brentq(
                lambda psi: st.nchypergeom_fisher.sf(a - 1, M=total, n=n,
                                                      odds=psi, N=k) - alpha2,
                1e-12, 1e12)
        except ValueError:
            ci_low = 0.0

    if b == 0 or c == 0:
        ci_high = float("inf")
    else:
        try:
            ci_high = brentq(
                lambda psi: st.nchypergeom_fisher.cdf(a, M=total, n=n,
                                                       odds=psi, N=k) - alpha2,
                1e-12, 1e12)
        except ValueError:
            ci_high = float("inf")

    return {
        "OR": fmt_r(or_est) if math.isfinite(or_est) else ("Inf" if or_est > 0 else "0"),
        "CI_low": fmt_r(ci_low) if math.isfinite(ci_low) else "Inf",
        "CI_high": fmt_r(ci_high) if math.isfinite(ci_high) else "Inf",
        "p": fmt_r(pval),
        "status": "TESTED",
        "reason": "",
    }

# ---------------------------------------------------------------------------
# BH FDR
# ---------------------------------------------------------------------------
def bh(pvals):
    m = len(pvals)
    order = sorted(range(m), key=lambda i: pvals[i])
    q = [0.0] * m
    prev = 1.0
    for rank in range(m, 0, -1):
        i = order[rank - 1]
        prev = min(prev, pvals[i] * m / rank)
        q[i] = prev
    return q

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    contract = json.loads(CONTRACT.read_text())
    assoc_cfg = contract["association"]

    # ---- metadata ----------------------------------------------------------------
    overall_strains = []
    comparison_strains = []
    group_overall = {}
    group_comparison = {}
    boundary_strains = []

    with open(META_OVERALL) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            overall_strains.append(r["strain"])
            group_overall[r["strain"]] = r["group_for_sensitivity152"] or "unassigned"
            if r["included_species_comparison144"] == "1":
                comparison_strains.append(r["strain"])
                group_comparison[r["strain"]] = r["group_for_testing"] or "unassigned"
            else:
                boundary_strains.append(r["strain"])

    assert len(overall_strains) == 152, len(overall_strains)
    assert len(comparison_strains) == 144, len(comparison_strains)
    assert len(boundary_strains) == 8, len(boundary_strains)
    assert set(comparison_strains).isdisjoint(set(boundary_strains))

    comparison_groups = {g: [s for s in comparison_strains
                             if group_comparison[s] == g] for g in GROUPS}
    sizes144 = {g: len(comparison_groups[g]) for g in GROUPS}
    assert sizes144 == TARGET_GROUP_SIZES_PRIMARY, sizes144

    # ---- partition table ---------------------------------------------------------
    partition = {}
    with open(PARTITION_TABLE) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            partition[r["cluster"]] = r["partition_recomputed"]

    # ---- cluster COG map ---------------------------------------------------------
    cmap = {}
    with open(CLUSTER_COG_MAP) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            cmap[r["cluster"]] = r
    annotated = {cl for cl, r in cmap.items() if r["annotation_state"] == "annotated"}
    cog_letters = {cl: set(cmap[cl]["cog_letters"]) - {"-"} for cl in annotated}

    # ---- PAV binary --------------------------------------------------------------
    ledger = []
    tests = []
    cluster_group_counts = {}

    with open(PAV_BINARY) as fh:
        rdr = csv.reader(fh, delimiter="\t")
        header = next(rdr)
        acc_idx = {acc: i for i, acc in enumerate(header[1:])}

        for row in rdr:
            cluster = row[0]
            vals = [int(x) for x in row[1:]]

            # carrier counts across scopes
            carrier152 = sum(vals)
            carrier144 = sum(vals[acc_idx[s]] for s in comparison_strains)
            carrier8 = sum(vals[acc_idx[s]] for s in boundary_strains)

            counts144 = {g: sum(vals[acc_idx[s]] for s in comparison_groups[g])
                         for g in GROUPS}
            cluster_group_counts[cluster] = counts144

            part = partition.get(cluster, "NA")

            if part not in ("Shell", "Cloud"):
                ledger.append({
                    "cluster": cluster, "partition": part,
                    "carrier144": carrier144, "carrier152": carrier152,
                    "carrier8": carrier8,
                    "evaluable_for_association": False,
                    "status": "FILTER_partition_not_shell_cloud",
                    "reason": f"partition={part}; only Shell/Cloud are candidates",
                })
                continue

            if carrier144 < assoc_cfg["carrier_minimum"]:
                ledger.append({
                    "cluster": cluster, "partition": part,
                    "carrier144": carrier144, "carrier152": carrier152,
                    "carrier8": carrier8,
                    "evaluable_for_association": False,
                    "status": "FILTER_carrier_lt5",
                    "reason": f"carrier144={carrier144} < {assoc_cfg['carrier_minimum']}",
                })
                continue

            ledger.append({
                "cluster": cluster, "partition": part,
                "carrier144": carrier144, "carrier152": carrier152,
                "carrier8": carrier8,
                "evaluable_for_association": True,
                "status": "ELIGIBLE",
                "reason": "Shell/Cloud and carrier>=5 in primary144 evaluable genomes",
            })

            for tg in GROUPS:
                rest = [g for g in GROUPS if g != tg]
                a = counts144[tg]
                b = sizes144[tg] - a
                c = sum(counts144[g] for g in rest)
                d = sum(sizes144[g] for g in rest) - c
                tests.append({
                    "feature": cluster,
                    "target_group": tg,
                    "a": a, "b": b, "c": c, "d": d,
                    "n_target_evaluable": sizes144[tg],
                    "n_rest_evaluable": sum(sizes144[g] for g in rest),
                    "target_frequency": a / sizes144[tg],
                    "rest_frequency": c / sum(sizes144[g] for g in rest),
                    "partition": part,
                    "carrier144": carrier144,
                })

    # ---- write ledger ------------------------------------------------------------
    with open(OUT_DIR / "cluster_filter_ledger.tsv", "w", newline="") as fh:
        cols = ["cluster", "partition", "carrier144", "carrier152", "carrier8",
                "evaluable_for_association", "status", "reason"]
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=cols, lineterminator="\n")
        w.writeheader()
        w.writerows(ledger)

    n_eligible = sum(1 for r in ledger if r["evaluable_for_association"])

    # ---- association Fisher tests ------------------------------------------------
    fisher_results = []
    for t in tests:
        res = fisher_test(t["a"], t["b"], t["c"], t["d"])
        fisher_results.append(res)

    # BH over family gene_species_main (all valid tests)
    valid_ps = [float(r["p"]) if r["status"] == "TESTED" else None
                for r in fisher_results]
    valid_indices = [i for i, p in enumerate(valid_ps) if p is not None]
    qs = bh([valid_ps[i] for i in valid_indices])
    qmap = {valid_indices[i]: qs[i] for i in range(len(valid_indices))}

    association_rows = []
    for i, t in enumerate(tests):
        f = fisher_results[i]
        row = {
            "feature": t["feature"],
            "target_group": t["target_group"],
            "a": t["a"], "b": t["b"], "c": t["c"], "d": t["d"],
            "n_target_evaluable": t["n_target_evaluable"],
            "n_rest_evaluable": t["n_rest_evaluable"],
            "target_frequency": fmt_r(t["target_frequency"]),
            "rest_frequency": fmt_r(t["rest_frequency"]),
            "OR": f["OR"],
            "CI_low": f["CI_low"],
            "CI_high": f["CI_high"],
            "p": f["p"],
            "q": fmt_r(qmap[i]) if i in qmap else "NA",
            "status": f["status"],
            "reason": f["reason"],
        }
        association_rows.append(row)

    with open(OUT_DIR / "association_all_tests.tsv", "w", newline="") as fh:
        cols = ["feature", "target_group", "a", "b", "c", "d",
                "n_target_evaluable", "n_rest_evaluable",
                "target_frequency", "rest_frequency",
                "OR", "CI_low", "CI_high", "p", "q", "status", "reason"]
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=cols, lineterminator="\n")
        w.writeheader()
        w.writerows(association_rows)

    # ---- associated candidates ranked -------------------------------------------
    state_rank = {"annotated": 0, "matched_no_cog": 1,
                  "matched_no_emapper_row": 2, "no_sequence_match": 3}
    candidates = []
    for i, t in enumerate(tests):
        if fisher_results[i]["status"] != "TESTED":
            continue
        tf = t["target_frequency"]
        rf = t["rest_frequency"]
        q = qmap[i]
        if tf > assoc_cfg["candidate_target_frequency_gt"] and \
           rf < assoc_cfg["candidate_rest_frequency_lt"] and \
           q < assoc_cfg["candidate_q_lt"]:
            cm = cmap.get(t["feature"], {})
            candidates.append({
                "target_group": t["target_group"],
                "cluster": t["feature"],
                "prevalence_difference": tf - rf,
                "q": q,
                "product": cm.get("product", "NA"),
                "annotation_status": cm.get("annotation_state", "NA"),
                "cog_letters": cm.get("cog_letters", "NA"),
            })

    candidates.sort(key=lambda x: (-x["prevalence_difference"], x["q"],
                                   state_rank.get(x["annotation_status"], 9),
                                   x["cluster"]))

    with open(OUT_DIR / "associated_candidates_ranked.tsv", "w", newline="") as fh:
        cols = ["target_group", "cluster", "rank", "prevalence_difference",
                "q", "product", "annotation_status", "cog_letters"]
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=cols, lineterminator="\n")
        w.writeheader()
        for rank, c in enumerate(candidates, 1):
            w.writerow({
                "target_group": c["target_group"],
                "cluster": c["cluster"],
                "rank": rank,
                "prevalence_difference": fmt_r(c["prevalence_difference"]),
                "q": fmt_r(c["q"]),
                "product": c["product"],
                "annotation_status": c["annotation_status"],
                "cog_letters": c["cog_letters"],
            })

    cand_by_group = defaultdict(set)
    for c in candidates:
        cand_by_group[c["target_group"]].add(c["cluster"])

    # ---- COG enrichment ----------------------------------------------------------
    eligible_clusters = {r["cluster"] for r in ledger
                         if r["evaluable_for_association"]}
    bg_all = eligible_clusters & annotated

    cog_tests = []
    cog_test_meta = []

    for g in GROUPS:
        fg = cand_by_group.get(g, set()) & annotated
        set_id = f"assoc_{g}"
        if not fg:
            # still emit rows for bookkeeping, but mark not estimable
            for cat in sorted({l for cl in bg_all for l in cog_letters[cl]}):
                cog_tests.append({
                    "family_id": "COG_species_associated",
                    "set_id": set_id, "category": cat,
                    "a": 0, "b": 0, "c": 0, "d": 0,
                    "fg_clusters": sorted(fg), "bg_clusters": sorted(bg_all),
                    "status": "NOT_ESTIMABLE",
                    "reason": "empty foreground (no annotated associated candidates)",
                })
            continue
        bg = bg_all - fg
        cats = sorted({l for cl in fg | bg for l in cog_letters[cl]})
        for cat in cats:
            a = sum(1 for cl in fg if cat in cog_letters[cl])
            c = sum(1 for cl in bg if cat in cog_letters[cl])
            cog_tests.append({
                "family_id": "COG_species_associated",
                "set_id": set_id, "category": cat,
                "a": a, "b": len(fg) - a, "c": c, "d": len(bg) - c,
                "fg_clusters": sorted(fg), "bg_clusters": sorted(bg),
                "status": "TESTED", "reason": "",
            })

    cog_fisher_results = []
    for t in cog_tests:
        if t["status"] == "TESTED":
            res = fisher_test(t["a"], t["b"], t["c"], t["d"])
        else:
            res = {"OR": "NA", "CI_low": "NA", "CI_high": "NA",
                   "p": "NA", "status": t["status"], "reason": t["reason"]}
        cog_fisher_results.append(res)

    # BH over COG_species_associated family
    valid_cog_idx = [i for i, r in enumerate(cog_fisher_results)
                     if r["status"] == "TESTED"]
    cog_qs = bh([float(cog_fisher_results[i]["p"]) for i in valid_cog_idx])
    cog_qmap = {valid_cog_idx[i]: cog_qs[i] for i in range(len(valid_cog_idx))}

    with open(OUT_DIR / "COG_enrichment_all.tsv", "w", newline="") as fh:
        cols = ["family_id", "set_id", "category",
                "foreground_count", "foreground_total",
                "background_count", "background_total",
                "OR", "CI_low", "CI_high", "p", "q", "status", "reason"]
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=cols, lineterminator="\n")
        w.writeheader()
        for i, t in enumerate(cog_tests):
            f = cog_fisher_results[i]
            if t["status"] == "TESTED":
                w.writerow({
                    "family_id": t["family_id"],
                    "set_id": t["set_id"],
                    "category": t["category"],
                    "foreground_count": t["a"],
                    "foreground_total": t["a"] + t["b"],
                    "background_count": t["c"],
                    "background_total": t["c"] + t["d"],
                    "OR": f["OR"], "CI_low": f["CI_low"], "CI_high": f["CI_high"],
                    "p": f["p"],
                    "q": fmt_r(cog_qmap[i]),
                    "status": f["status"], "reason": f["reason"],
                })
            else:
                w.writerow({
                    "family_id": t["family_id"],
                    "set_id": t["set_id"],
                    "category": t["category"],
                    "foreground_count": 0, "foreground_total": 0,
                    "background_count": 0, "background_total": 0,
                    "OR": "NA", "CI_low": "NA", "CI_high": "NA",
                    "p": "NA", "q": "NA",
                    "status": f["status"], "reason": f["reason"],
                })

    # ---- foreground/background membership ---------------------------------------
    with open(OUT_DIR / "cog_foreground_background.tsv", "w", newline="") as fh:
        cols = ["cluster", "COG_category", "is_foreground", "is_background", "target_group"]
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=cols, lineterminator="\n")
        w.writeheader()
        for t in cog_tests:
            g = t["set_id"].replace("assoc_", "")
            for cl in t["fg_clusters"]:
                for cat in sorted(cog_letters.get(cl, set())):
                    w.writerow({
                        "cluster": cl, "COG_category": cat,
                        "is_foreground": True, "is_background": False,
                        "target_group": g,
                    })
            for cl in t["bg_clusters"]:
                for cat in sorted(cog_letters.get(cl, set())):
                    w.writerow({
                        "cluster": cl, "COG_category": cat,
                        "is_foreground": False, "is_background": True,
                        "target_group": g,
                    })

    # ---- run_info.json ------------------------------------------------------------
    n_cog_tested = sum(1 for r in cog_fisher_results if r["status"] == "TESTED")
    sig_cog = []
    for i in valid_cog_idx:
        if cog_qmap[i] < 0.05:
            or_str = cog_fisher_results[i]["OR"]
            if or_str == "Inf":
                direction = "enriched_inf"
            elif or_str == "0":
                direction = "depleted_zero"
            else:
                direction = "enriched" if float(or_str) > 1 else "depleted"
            sig_cog.append((cog_tests[i]["set_id"], cog_tests[i]["category"],
                            or_str, direction, cog_qmap[i]))
    enriched_cats = sorted({f"{sid}:{cat}(OR={or_str})"
                            for sid, cat, or_str, direction, _ in sig_cog
                            if direction == "enriched"})
    depleted_cats = sorted({f"{sid}:{cat}(OR={or_str})"
                            for sid, cat, or_str, direction, _ in sig_cog
                            if direction in ("depleted", "depleted_zero")})
    enriched_inf_cats = sorted({f"{sid}:{cat}(OR={or_str})"
                                for sid, cat, or_str, direction, _ in sig_cog
                                if direction == "enriched_inf"})

    run_info = {
        "module": "C",
        "date": contract.get("date"),
        "dictionary_id": contract.get("dictionary_id"),
        "parameters": {
            "carrier_minimum": assoc_cfg["carrier_minimum"],
            "carrier_denominator": assoc_cfg["carrier_denominator"],
            "fisher_alternative": assoc_cfg["fisher_alternative"],
            "OR": assoc_cfg["OR"],
            "CI": assoc_cfg["CI"],
            "candidate_target_frequency_gt": assoc_cfg["candidate_target_frequency_gt"],
            "candidate_rest_frequency_lt": assoc_cfg["candidate_rest_frequency_lt"],
            "candidate_q_lt": assoc_cfg["candidate_q_lt"],
            "use_unrounded_thresholds": assoc_cfg["use_unrounded_thresholds"],
        },
        "bh_families": {
            "gene_species_main": len(valid_indices),
            "COG_species_associated": len(valid_cog_idx),
        },
        "counts": {
            "clusters_total": len(ledger),
            "clusters_eligible_for_association": n_eligible,
            "association_tests_total": len(tests),
            "association_tests_tested": len(valid_indices),
            "association_tests_not_estimable": len(tests) - len(valid_indices),
            "associated_candidates_total": len(candidates),
            "associated_candidates_by_group": {g: len(cand_by_group.get(g, set()))
                                               for g in GROUPS},
            "cog_enrichment_tests_total": len(cog_tests),
            "cog_enrichment_tests_tested": n_cog_tested,
            "cog_enrichment_q_lt_0.05": len(sig_cog),
            "cog_enriched_categories": enriched_cats,
            "cog_enriched_inf_categories": enriched_inf_cats,
            "cog_depleted_categories": depleted_cats,
        },
        "status": "PASS",
    }
    with open(OUT_DIR / "run_info.json", "w") as fh:
        json.dump(run_info, fh, indent=2)

    # ---- stderr summary ----------------------------------------------------------
    sys.stderr.write(json.dumps({
        "clusters_eligible": n_eligible,
        "association_tests": len(tests),
        "candidates": {g: len(cand_by_group.get(g, set())) for g in GROUPS},
        "cog_enriched_categories": enriched_cats,
    }, indent=1) + "\n")

    return 0

if __name__ == "__main__":
    sys.exit(main())
