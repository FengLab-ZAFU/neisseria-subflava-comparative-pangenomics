#!/usr/bin/env python3
"""Module D1 (AMR / tet(M)) for boundary8_split_handoff_20261019.

Produces the AMR-related outputs required in S/07_amr_tetm per
boundary8_split_handoff_20261019/00_cohort/required_outputs.tsv module D.
"""

import csv
import json
import math
import os
import sys

from scipy import stats
from scipy.optimize import brentq
from scipy.special import gammaln
from statsmodels.stats.multitest import multipletests

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import P, R, S, read_tsv, write_tsv, write_json

OUT = os.path.join(S, "07_amr_tetm")
PRIMARY_GROUPS = ["subflava", "perflava", "flavescens"]


def fmt6(x):
    if x is None or (isinstance(x, float) and math.isnan(x)):
        return "NA"
    return f"{x:.6f}"


def fmt_or_na(x):
    if x is None or (isinstance(x, float) and math.isnan(x)):
        return "NA"
    if hasattr(x, "item"):
        x = x.item()
    return repr(x)


def wilson_ci(k, n, alpha=0.05):
    """Wilson score interval clipped to [0, 1]."""
    if n == 0:
        return (float("nan"), float("nan"))
    p = k / n
    z = stats.norm.ppf(1 - alpha / 2)
    denom = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / denom
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / denom
    return (max(0.0, centre - half), min(1.0, centre + half))


# -----------------------------------------------------------------------------
# Conditional MLE OR + exact CI for 2x2 tables (noncentral hypergeometric)
# -----------------------------------------------------------------------------

def _nch_support(n1, n2, m1):
    return max(0, m1 - n2), min(n1, m1)


def _nch_terms(n1, n2, m1, logpsi):
    k_min, k_max = _nch_support(n1, n2, m1)
    ks = list(range(k_min, k_max + 1))
    logw = []
    for k in ks:
        lw = (
            gammaln(n1 + 1) - gammaln(k + 1) - gammaln(n1 - k + 1)
            + gammaln(n2 + 1) - gammaln(m1 - k + 1) - gammaln(n2 - m1 + k + 1)
            + k * logpsi
        )
        logw.append(lw)
    maxlog = max(logw)
    w = [math.exp(lw - maxlog) for lw in logw]
    s = sum(w)
    p = [wi / s for wi in w]
    return ks, p


def conditional_mle_or(a, b, c, d):
    n1 = a + b
    n2 = c + d
    m1 = a + c
    k_min, k_max = _nch_support(n1, n2, m1)
    if a <= k_min:
        return 0.0
    if a >= k_max:
        return float("inf")

    def f(logpsi):
        ks, p = _nch_terms(n1, n2, m1, logpsi)
        return sum(k * pk for k, pk in zip(ks, p)) - a

    lo, hi = -20.0, 20.0
    while f(lo) > 0 and lo > -1000:
        lo *= 2
    while f(hi) < 0 and hi < 1000:
        hi *= 2
    return math.exp(brentq(f, lo, hi))


def conditional_or_ci(a, b, c, d, alpha=0.05):
    n1 = a + b
    n2 = c + d
    m1 = a + c
    k_min, k_max = _nch_support(n1, n2, m1)

    if a == k_min:
        lower = 0.0
    else:
        def fl(logpsi):
            ks, p = _nch_terms(n1, n2, m1, logpsi)
            return sum(pk for k, pk in zip(ks, p) if k >= a) - alpha / 2
        lo, hi = -20.0, 20.0
        while fl(lo) > 0 and lo > -1000:
            lo *= 2
        while fl(hi) < 0 and hi < 1000:
            hi *= 2
        lower = math.exp(brentq(fl, lo, hi))

    if a == k_max:
        upper = float("inf")
    else:
        def fu(logpsi):
            ks, p = _nch_terms(n1, n2, m1, logpsi)
            return sum(pk for k, pk in zip(ks, p) if k <= a) - alpha / 2
        lo, hi = -20.0, 20.0
        while fu(lo) > 0 and lo > -1000:
            lo *= 2
        while fu(hi) < 0 and hi < 1000:
            hi *= 2
        upper = math.exp(brentq(fu, lo, hi))

    return lower, upper


def fisher_exact_2x2(a, b, c, d):
    """Return (OR_MLE, CI_low, CI_high, p_two_sided)."""
    _, p = stats.fisher_exact([[a, b], [c, d]], alternative="two-sided")
    or_mle = conditional_mle_or(a, b, c, d)
    ci_low, ci_high = conditional_or_ci(a, b, c, d)
    return or_mle, ci_low, ci_high, p


def main():
    os.makedirs(OUT, exist_ok=True)

    contract = json.load(open(os.path.join(P, "00_cohort", "analysis_contract.json")))
    meta_overall = {r["strain"]: r for r in read_tsv(os.path.join(P, "00_cohort", "metadata_overall152.tsv"))}
    matrix = read_tsv(os.path.join(P, "02_evidence", "AMR_feature_evaluable_matrix152.tsv"))
    manifest = read_tsv(os.path.join(S, "08_defense_mge_vfdb", "participant_subset_manifest.tsv"))
    per_genome = {r["accession"]: r for r in read_tsv(os.path.join(R, "07_amr_tetm", "amr_tetm_per_genome.tsv"))}

    features = []
    for r in matrix:
        if r["feature"] not in features:
            features.append(r["feature"])
    mat_lookup = {(r["accession"], r["feature"]): r for r in matrix}

    overall_accs = list(meta_overall.keys())
    primary_accs = [acc for acc in overall_accs if meta_overall[acc]["included_species_comparison144"] == "1"]
    desc_accs = [acc for acc in overall_accs if meta_overall[acc]["included_species_comparison144"] == "0"]

    # -------------------------------------------------------------------------
    # 1. carriage_by_group.tsv
    # -------------------------------------------------------------------------
    carriage_rows = []
    for feat in features:
        # overall152
        eval_n = pos_n = 0
        for acc in overall_accs:
            row = mat_lookup[(acc, feat)]
            eval_n += int(row["evaluable"])
            pos_n += int(row["positive"])
        p = pos_n / eval_n if eval_n else float("nan")
        lo, hi = wilson_ci(pos_n, eval_n)
        carriage_rows.append({
            "analysis_scope": "overall152",
            "group": "overall152",
            "feature": feat,
            "n_total": 152,
            "n_evaluable": eval_n,
            "n_positive": pos_n,
            "proportion": fmt6(p),
            "wilson_low": fmt6(lo),
            "wilson_high": fmt6(hi),
            "status": "OK",
            "reason": "",
        })
        # primary144 groups
        for grp in PRIMARY_GROUPS:
            accs = [acc for acc in primary_accs if meta_overall[acc]["group_for_testing"] == grp]
            n_total = len(accs)
            eval_n = pos_n = 0
            for acc in accs:
                row = mat_lookup[(acc, feat)]
                eval_n += int(row["evaluable"])
                pos_n += int(row["positive"])
            p = pos_n / eval_n if eval_n else float("nan")
            lo, hi = wilson_ci(pos_n, eval_n)
            carriage_rows.append({
                "analysis_scope": "primary144",
                "group": grp,
                "feature": feat,
                "n_total": n_total,
                "n_evaluable": eval_n,
                "n_positive": pos_n,
                "proportion": fmt6(p),
                "wilson_low": fmt6(lo),
                "wilson_high": fmt6(hi),
                "status": "OK",
                "reason": "",
            })
        # descriptive8
        n_total = len(desc_accs)
        eval_n = pos_n = 0
        for acc in desc_accs:
            row = mat_lookup[(acc, feat)]
            eval_n += int(row["evaluable"])
            pos_n += int(row["positive"])
        p = pos_n / eval_n if eval_n else float("nan")
        lo, hi = wilson_ci(pos_n, eval_n)
        carriage_rows.append({
            "analysis_scope": "descriptive8",
            "group": "descriptive8",
            "feature": feat,
            "n_total": n_total,
            "n_evaluable": eval_n,
            "n_positive": pos_n,
            "proportion": fmt6(p),
            "wilson_low": fmt6(lo),
            "wilson_high": fmt6(hi),
            "status": "OK",
            "reason": "",
        })

    write_tsv(os.path.join(OUT, "carriage_by_group.tsv"), carriage_rows,
              ["analysis_scope", "group", "feature", "n_total", "n_evaluable",
               "n_positive", "proportion", "wilson_low", "wilson_high",
               "status", "reason"])

    # -------------------------------------------------------------------------
    # 2. amr_association_all.tsv
    # -------------------------------------------------------------------------
    assoc_rows = []
    valid_p = {}
    for feat in features:
        carriers_144 = 0
        eval_144 = 0
        for acc in primary_accs:
            row = mat_lookup[(acc, feat)]
            eval_144 += int(row["evaluable"])
            carriers_144 += int(row["positive"])
        filtered = carriers_144 < 5
        for grp in PRIMARY_GROUPS:
            accs_target = [acc for acc in primary_accs if meta_overall[acc]["group_for_testing"] == grp]
            a = b = 0
            for acc in accs_target:
                row = mat_lookup[(acc, feat)]
                if int(row["evaluable"]):
                    b += 1
                    if int(row["positive"]):
                        a += 1
            b = b - a
            c = carriers_144 - a
            d = (eval_144 - (a + b)) - c
            target_n = a + b
            rest_n = c + d
            target_freq = a / target_n if target_n else float("nan")
            rest_freq = c / rest_n if rest_n else float("nan")

            if filtered:
                status = "FILTERED"
                reason = "carrier<5_in_primary144"
                or_val = ci_low = ci_high = p = q = "NA"
            elif target_n == 0 or rest_n == 0 or (a == 0 and c == 0) or (a == target_n and c == rest_n):
                status = "NOT_ESTIMABLE"
                reason = "no_variation_or_empty_margin"
                or_val = ci_low = ci_high = p = q = "NA"
            else:
                status = "VALID"
                reason = ""
                or_val, ci_low, ci_high, p = fisher_exact_2x2(a, b, c, d)
                valid_p[(feat, grp)] = p
                or_val = fmt_or_na(or_val)
                ci_low = fmt_or_na(ci_low)
                ci_high = fmt_or_na(ci_high)
                p = fmt_or_na(p)
                q = "NA"

            assoc_rows.append({
                "family_id": "AMR_prevalence",
                "feature": feat,
                "target_group": grp,
                "a": a,
                "b": b,
                "c": c,
                "d": d,
                "target_n": target_n,
                "rest_n": rest_n,
                "target_frequency": fmt_or_na(target_freq),
                "rest_frequency": fmt_or_na(rest_freq),
                "OR": or_val,
                "CI_low": ci_low,
                "CI_high": ci_high,
                "p": p,
                "q": q,
                "status": status,
                "reason": reason,
            })

    if valid_p:
        keys = list(valid_p.keys())
        ps = [valid_p[k] for k in keys]
        qs = multipletests(ps, method="fdr_bh")[1]
        q_map = dict(zip(keys, qs))
        for r in assoc_rows:
            if r["status"] == "VALID":
                r["q"] = fmt_or_na(q_map[(r["feature"], r["target_group"])])

    write_tsv(os.path.join(OUT, "amr_association_all.tsv"), assoc_rows,
              ["family_id", "feature", "target_group", "a", "b", "c", "d",
               "target_n", "rest_n", "target_frequency", "rest_frequency",
               "OR", "CI_low", "CI_high", "p", "q", "status", "reason"])

    # -------------------------------------------------------------------------
    # 3. feature_evaluable_matrix.tsv
    # -------------------------------------------------------------------------
    fem_rows = []
    for r in matrix:
        acc = r["accession"]
        m = meta_overall.get(acc, {})
        fem_rows.append({
            "accession": acc,
            "feature": r["feature"],
            "terminal_state": r["terminal_state"],
            "positive": r["positive"],
            "evaluable": r["evaluable"],
            "loose_only": r["loose_only"],
            "RGI_run_status": r["RGI_run_status"],
            "definition": r["definition"],
            "group_for_testing": m.get("group_for_testing", ""),
            "display_group": m.get("display_group", ""),
        })
    write_tsv(os.path.join(OUT, "feature_evaluable_matrix.tsv"), fem_rows,
              ["accession", "feature", "terminal_state", "positive", "evaluable",
               "loose_only", "RGI_run_status", "definition", "group_for_testing",
               "display_group"])

    # -------------------------------------------------------------------------
    # 4. sensitivity_amr_tetm.tsv (participant subset, n=95)
    # -------------------------------------------------------------------------
    selected = {r["accession"]: r for r in manifest if r["selected"] == "1"}
    sens_rows = []
    valid_p_subset = {}
    for feat in features:
        carriers_subset = 0
        eval_subset = 0
        for acc in selected:
            row = mat_lookup[(acc, feat)]
            eval_subset += int(row["evaluable"])
            carriers_subset += int(row["positive"])
        filtered_subset = carriers_subset < 5
        for grp in PRIMARY_GROUPS:
            accs = [acc for acc, r in selected.items() if r["species_group"] == grp]
            n_total = len(accs)
            eval_n = pos_n = 0
            for acc in accs:
                row = mat_lookup[(acc, feat)]
                eval_n += int(row["evaluable"])
                pos_n += int(row["positive"])
            p = pos_n / eval_n if eval_n else float("nan")
            lo, hi = wilson_ci(pos_n, eval_n)

            a = pos_n
            b = eval_n - pos_n
            c = carriers_subset - a
            d = (eval_subset - eval_n) - c
            target_n = a + b
            rest_n = c + d
            target_freq = a / target_n if target_n else float("nan")
            rest_freq = c / rest_n if rest_n else float("nan")

            if filtered_subset:
                status = "FILTERED"
                reason = "carrier<5_in_participant_subset"
                or_val = ci_low = ci_high = pval = q = "NA"
            elif target_n == 0 or rest_n == 0 or (a == 0 and c == 0) or (a == target_n and c == rest_n):
                status = "NOT_ESTIMABLE"
                reason = "no_variation_or_empty_margin"
                or_val = ci_low = ci_high = pval = q = "NA"
            else:
                status = "VALID"
                reason = ""
                or_val, ci_low, ci_high, pval = fisher_exact_2x2(a, b, c, d)
                valid_p_subset[(feat, grp)] = pval
                or_val = fmt_or_na(or_val)
                ci_low = fmt_or_na(ci_low)
                ci_high = fmt_or_na(ci_high)
                pval = fmt_or_na(pval)
                q = "NA"

            sens_rows.append({
                "feature": feat,
                "target_group": grp,
                "n_total": n_total,
                "n_evaluable": eval_n,
                "n_positive": pos_n,
                "proportion": fmt6(p),
                "wilson_low": fmt6(lo),
                "wilson_high": fmt6(hi),
                "target_frequency": fmt_or_na(target_freq),
                "rest_frequency": fmt_or_na(rest_freq),
                "OR": or_val,
                "CI_low": ci_low,
                "CI_high": ci_high,
                "p": pval,
                "q": q,
                "status": status,
                "reason": reason,
            })

    if valid_p_subset:
        keys = list(valid_p_subset.keys())
        ps = [valid_p_subset[k] for k in keys]
        qs = multipletests(ps, method="fdr_bh")[1]
        q_map = dict(zip(keys, qs))
        for r in sens_rows:
            if r["status"] == "VALID":
                r["q"] = fmt_or_na(q_map[(r["feature"], r["target_group"])])

    write_tsv(os.path.join(OUT, "sensitivity_amr_tetm.tsv"), sens_rows,
              ["feature", "target_group", "n_total", "n_evaluable", "n_positive",
               "proportion", "wilson_low", "wilson_high", "target_frequency",
               "rest_frequency", "OR", "CI_low", "CI_high", "p", "q", "status",
               "reason"])

    # -------------------------------------------------------------------------
    # 5. run_info.json + benchmark check
    # -------------------------------------------------------------------------
    def tetm_counts(accs):
        eval_n = pos_n = 0
        for acc in accs:
            row = mat_lookup[(acc, "gene:tet(M)")]
            eval_n += int(row["evaluable"])
            pos_n += int(row["positive"])
        return {"positive": pos_n, "evaluable": eval_n}

    obs_genome = {
        "overall152": tetm_counts(overall_accs),
        "primary144": tetm_counts(primary_accs),
        "descriptive8": tetm_counts(desc_accs),
    }
    for grp in PRIMARY_GROUPS:
        accs = [acc for acc in primary_accs if meta_overall[acc]["group_for_testing"] == grp]
        obs_genome[grp] = tetm_counts(accs)

    def locus_counts(accs):
        return sum(int(per_genome[acc]["n_tetM_loci"]) for acc in accs if acc in per_genome)

    obs_locus = {
        "overall152": locus_counts(overall_accs),
        "primary144": locus_counts(primary_accs),
        "descriptive8": locus_counts(desc_accs),
    }
    for grp in PRIMARY_GROUPS:
        accs = [acc for acc in primary_accs if meta_overall[acc]["group_for_testing"] == grp]
        obs_locus[grp] = locus_counts(accs)

    expected = {
        "overall152": {"positive_genomes": 38, "loci": 57},
        "primary144": {"positive_genomes": 35, "loci": 54},
        "descriptive8": {"positive_genomes": 3, "loci": 3},
        "subflava": {"positive_genomes": 0, "evaluable_genomes": 23},
        "perflava": {"positive_genomes": 35, "evaluable_genomes": 114},
        "flavescens": {"positive_genomes": 0, "evaluable_genomes": 7},
    }

    benchmark_pass = True
    notes = []
    for scope in ["overall152", "primary144", "descriptive8"]:
        if obs_genome[scope]["positive"] != expected[scope]["positive_genomes"]:
            benchmark_pass = False
            notes.append(f"{scope} positive genome count {obs_genome[scope]['positive']} != expected {expected[scope]['positive_genomes']}")
        if obs_locus[scope] != expected[scope]["loci"]:
            benchmark_pass = False
            notes.append(f"{scope} tet(M) locus count {obs_locus[scope]} != expected {expected[scope]['loci']}")
    for grp in PRIMARY_GROUPS:
        if obs_genome[grp]["positive"] != expected[grp]["positive_genomes"]:
            benchmark_pass = False
            notes.append(f"{grp} positive genome count {obs_genome[grp]['positive']} != expected {expected[grp]['positive_genomes']}")
        if obs_genome[grp]["evaluable"] != expected[grp]["evaluable_genomes"]:
            benchmark_pass = False
            notes.append(f"{grp} evaluable genome count {obs_genome[grp]['evaluable']} != expected {expected[grp]['evaluable_genomes']}")

    run_info = {
        "module": "D1",
        "status": "COMPLETE",
        "date": contract.get("date", "2026-10-19"),
        "expected_tetM_counts": expected,
        "observed_tetM_genome_counts": obs_genome,
        "observed_tetM_locus_counts": obs_locus,
        "benchmark_pass": benchmark_pass,
        "notes": "; ".join(notes) if notes else "Observed tet(M) counts match expected benchmarks.",
        "output_summary": {
            "carriage_by_group_rows": len(carriage_rows),
            "amr_association_all_rows": len(assoc_rows),
            "feature_evaluable_matrix_rows": len(fem_rows),
            "sensitivity_amr_tetm_rows": len(sens_rows),
            "n_features": len(features),
            "n_primary144_association_tests": len(assoc_rows),
            "n_valid_primary144_association_tests": len(valid_p),
            "n_participant_subset": len(selected),
        },
    }
    write_json(os.path.join(OUT, "run_info.json"), run_info)

    print("Module D1 outputs written to", OUT)
    print("Benchmark pass:", benchmark_pass)
    print("Observed tet(M) genome counts:", obs_genome)
    print("Observed tet(M) locus counts:", obs_locus)
    return 0 if benchmark_pass else 1


if __name__ == "__main__":
    sys.exit(main())
