#!/usr/bin/env python3
"""M5a finalize: join Fisher results, BH within family gene_species_main,
candidate filter + ranking. Also prepares the COG enrichment Fisher input
(partition family; species_associated family added after candidates known).
"""
from pathlib import Path
import os
import csv, json, sys
from collections import defaultdict

ROOT = str(Path(os.environ["NEISSERIA_PROJECT_ROOT"]))
RES = f"{ROOT}/results_revision_20261015_main152_species"
OUT = f"{RES}/06_function"
PGAP2 = f"{RES}/02_pangenome/pgap2_main152"
GROUPS = ["subflava", "perflava", "flavescens"]

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

def main():
    prep = list(csv.DictReader(open(f"{OUT}/association_rows_prep.tsv"), delimiter="\t"))
    fisher = {r["job_id"]: r for r in csv.DictReader(open(f"{OUT}/fisher_output_association.tsv"), delimiter="\t")}
    assert len(prep) == len(fisher)

    cmap = {}
    for r in csv.DictReader(open(f"{OUT}/cluster_cog_map.tsv"), delimiter="\t"):
        cmap[r["cluster"]] = r

    pvals = [float(fisher[f"assoc_{i}"]["p"]) for i in range(len(prep))]
    qs = bh(pvals)

    rows = []
    for i, r in enumerate(prep):
        f = fisher[f"assoc_{i}"]
        rows.append({
            "cluster": r["cluster"], "target_group": r["target_group"],
            "a": r["a"], "b": r["b"], "c": r["c"], "d": r["d"],
            "target_n": r["target_n"], "rest_n": r["rest_n"],
            "target_frequency": r["target_frequency"], "rest_frequency": r["rest_frequency"],
            "prevalence_difference": r["prevalence_difference"],
            "freq_subflava": r["freq_subflava"], "freq_perflava": r["freq_perflava"],
            "freq_flavescens": r["freq_flavescens"],
            "OR": f["OR"], "CI_low": f["CI_low"], "CI_high": f["CI_high"],
            "p": f["p"], "q": f"{qs[i]:.10g}",
            "family_id": "gene_species_main", "status": "TESTED",
        })

    cols = ["cluster", "target_group", "a", "b", "c", "d", "target_n", "rest_n",
            "target_frequency", "rest_frequency", "prevalence_difference",
            "freq_subflava", "freq_perflava", "freq_flavescens",
            "OR", "CI_low", "CI_high", "p", "q", "family_id", "status"]
    with open(f"{OUT}/association_all_tests.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=cols, lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    # candidates: pre-rounding thresholds
    state_rank = {"annotated": 0, "matched_no_cog": 1, "matched_no_emapper_row": 2,
                  "no_sequence_match": 3}
    cand = []
    for i, r in enumerate(prep):
        tf = float(r["target_frequency"]); rf = float(r["rest_frequency"])
        if tf > 0.80 and rf < 0.20 and qs[i] < 0.05:
            cm = cmap.get(r["cluster"], {})
            cand.append({
                "target_group": r["target_group"], "cluster": r["cluster"],
                "prevalence_difference": float(r["prevalence_difference"]),
                "q": qs[i],
                "product": cm.get("product", "NA"),
                "annotation_status": cm.get("annotation_state", "NA"),
            })
    cand.sort(key=lambda x: (-x["prevalence_difference"], x["q"],
                             state_rank.get(x["annotation_status"], 9), x["cluster"]))
    with open(f"{OUT}/associated_candidates_ranked.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["target_group", "cluster", "rank", "prevalence_difference", "q",
                    "product", "annotation_status"])
        for rank, c in enumerate(cand, 1):
            w.writerow([c["target_group"], c["cluster"], rank,
                        f"{c['prevalence_difference']:.10g}", f"{c['q']:.10g}",
                        c["product"], c["annotation_status"]])

    n_by_group = defaultdict(int)
    for c in cand:
        n_by_group[c["target_group"]] += 1
    summary = {"n_tests": len(rows), "n_candidates": len(cand),
               "n_candidates_by_group": dict(n_by_group),
               "candidate_clusters_by_group": {g: sorted(c["cluster"] for c in cand if c["target_group"] == g) for g in GROUPS}}
    json.dump(summary, open(f"{OUT}/association_summary.json", "w"), indent=1)
    sys.stderr.write(json.dumps({k: v for k, v in summary.items() if k != "candidate_clusters_by_group"}) + "\n")

if __name__ == "__main__":
    main()
