#!/usr/bin/env python3
"""Compare the 144-genome primary and 152-genome sensitivity scopes.

Pair results by feature, statistic and permutation design, retaining denominators,
effect sizes, P/q values, testability and multiplicity-family sizes. Analyses
without a comparable result are listed explicitly.
"""
from pathlib import Path
import csv, os, sys

ROOT = str(Path(os.environ["NEISSERIA_PROJECT_ROOT"]))
P = f"{ROOT}/boundary8_split_handoff_20261019"
S = f"{ROOT}/results_revision_20261019_boundary8_analysis"
F = f"{ROOT}/results_revision_20261019_boundary8_permdisp_fix"

def read_tsv(path):
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))

def write_tsv(path, fields, rows):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

FIELDS = ["module", "feature_or_item", "statistic", "design", "paired",
          "primary144_statistic", "sensitivity152_statistic",
          "primary144_effect_size", "sensitivity152_effect_size",
          "primary144_p", "sensitivity152_p", "primary144_q", "sensitivity152_q",
          "primary144_denominator", "sensitivity152_denominator",
          "primary144_status", "sensitivity152_status",
          "bh_family_id", "bh_family_size", "method_change", "grouping_change", "note"]

GRP_CHANGE = "subflava 23 (primary144) vs 31 (sensitivity152 reference assignment: 8 boundary genomes reassigned to their reference-ANI group); perflava 114 and flavescens 7 unchanged"
rows = []

def gap(module, item, stat, design, why):
    rows.append({"module": module, "feature_or_item": item, "statistic": stat, "design": design,
                 "paired": "NO_GAP", "primary144_statistic": "", "sensitivity152_statistic": "",
                 "primary144_effect_size": "", "sensitivity152_effect_size": "",
                 "primary144_p": "", "sensitivity152_p": "", "primary144_q": "", "sensitivity152_q": "",
                 "primary144_denominator": "", "sensitivity152_denominator": "",
                 "primary144_status": "", "sensitivity152_status": "",
                 "bh_family_id": "", "bh_family_size": "",
                 "method_change": "", "grouping_change": "", "note": why})

# ---------------- Module B: PERMANOVA + PERMDISP ----------------
ds_new = {r["analysis_id"]: r for r in read_tsv(f"{F}/05_structure/distance_statistics.tsv")}
for test in ("PERMANOVA", "PERMDISP"):
    for fs in ("shell_cloud", "all_variable"):
        for design in ("free", "restricted_within_study_namespace"):
            a = ds_new[f"{test.lower()}_{fs}_primary144_{design}"]
            b = ds_new[f"{test.lower()}_{fs}_sensitivity152_{design}"]
            if test == "PERMANOVA":
                mchg = "none: identical manual pseudo-F and P=(1+b)/(B+1) over delivered label configurations in both scopes"
                eff_a, eff_b = f"R2={a['R2']}", f"R2={b['R2']}"
                fam, famsz = "not_applicable_no_BH", "structure distance tests are not BH-adjusted by design"
            else:
                mchg = ("v1 (PARTIAL): group[P] label permutation over raw distances; "
                        "fix: vegan::permutest(bd, permutations=P_delivered) Freedman-Lane residual permutation — "
                        "identical corrected method in BOTH scopes, same delivered P matrices")
                eff_a, eff_b = "", ""
                fam, famsz = "not_applicable_no_BH", "structure distance tests are not BH-adjusted by design"
            rows.append({"module": "B_structure", "feature_or_item": fs, "statistic": test, "design": design,
                         "paired": "YES",
                         "primary144_statistic": f"F={a['statistic']}", "sensitivity152_statistic": f"F={b['statistic']}",
                         "primary144_effect_size": eff_a, "sensitivity152_effect_size": eff_b,
                         "primary144_p": a["p"], "sensitivity152_p": b["p"],
                         "primary144_q": "NA", "sensitivity152_q": "NA",
                         "primary144_denominator": f"n={a['n']}, df={a['df1']},{a['df2']}, B={a['B_used']}",
                         "sensitivity152_denominator": f"n={b['n']}, df={b['df1']},{b['df2']}, B={b['B_used']}",
                         "primary144_status": a["status"], "sensitivity152_status": b["status"],
                         "bh_family_id": fam, "bh_family_size": famsz,
                         "method_change": mchg, "grouping_change": GRP_CHANGE, "note": ""})

# ---------------- Module B: Mantel (GAP) ----------------
gap("B_structure", "patristic_vs_jaccard", "Mantel_Pearson", "free_whole_genome_labels",
    "Mantel is defined on the overall152 patristic tree only; there is no primary144 same-method counterpart (tree covers all 152 taxa). Reused unchanged r=0.7966 (shell_cloud) / 0.7372 (all_variable), P=0.001.")

# ---------------- Module A: matched sampling ----------------
mse = read_tsv(f"{S}/04_species_pangenome/matched_sampling_effects.tsv")
by = {}
for r in mse:
    by.setdefault((r["comparison"], r["metric"], r["group_a"], r["group_b"]), r)
    key = (r["comparison"], r["metric"], r["group_a"], r["group_b"])
for (comp, metric, ga, gb), r in sorted(by.items()):
    if comp.startswith("sensitivity152_"):
        base = comp[len("sensitivity152_"):]
        r144 = by.get((base, metric, ga, gb))
        rows.append({"module": "A_pangenome", "feature_or_item": f"matched_sampling_{base}", "statistic": metric,
                     "design": f"{ga}_vs_{gb}_n{r['n']}", "paired": "YES" if r144 else "NO_GAP",
                     "primary144_statistic": f"mean_diff={r144['mean_difference_a_minus_b']}" if r144 else "",
                     "sensitivity152_statistic": f"mean_diff={r['mean_difference_a_minus_b']}",
                     "primary144_effect_size": f"cohens_d={r144['cohens_d']}" if r144 else "",
                     "sensitivity152_effect_size": f"cohens_d={r['cohens_d']}",
                     "primary144_p": "NA", "sensitivity152_p": "NA", "primary144_q": "NA", "sensitivity152_q": "NA",
                     "primary144_denominator": f"n={r144['n']} per group, 1000 iterations" if r144 else "",
                     "sensitivity152_denominator": f"n={r['n']} per group, 1000 iterations",
                     "primary144_status": "DESCRIPTIVE" if r144 else "", "sensitivity152_status": "DESCRIPTIVE",
                     "bh_family_id": "not_applicable_no_test", "bh_family_size": "descriptive effect sizes, no p-values",
                     "method_change": "none: identical matched-size random-order accumulation", "grouping_change": GRP_CHANGE,
                     "note": "descriptive effect size over 1000 random-order iterations; not an independent-sample CI"})

# ---------------- Module D1: AMR carriage ----------------
fem = read_tsv(f"{S}/07_amr_tetm/feature_evaluable_matrix.tsv")
meta = {r["strain"]: r for r in read_tsv(f"{P}/00_cohort/metadata_overall152.tsv")}
meta144 = {r["strain"] for r in read_tsv(f"{P}/00_cohort/metadata_comparison144.tsv")}
feats = sorted({r["feature"] for r in fem})
groups = ("subflava", "perflava", "flavescens")
for feat in feats:
    for g in groups:
        rows144 = [r for r in fem if r["feature"] == feat and r["accession"] in meta144
                   and meta[r["accession"]]["group_for_testing"] == g]
        rows152 = [r for r in fem if r["feature"] == feat
                   and meta[r["accession"]]["group_for_sensitivity152"] == g]
        def cnt(rr):
            ev = [x for x in rr if x["evaluable"] == "1"]
            pos = sum(int(x["positive"]) for x in ev)
            return len(rr), len(ev), pos
        t144, e144, p144 = cnt(rows144)
        t152, e152, p152 = cnt(rows152)
        rows.append({"module": "D1_AMR", "feature_or_item": feat, "statistic": "carriage_prevalence", "design": f"group={g}",
                     "paired": "YES",
                     "primary144_statistic": f"{p144}/{e144}", "sensitivity152_statistic": f"{p152}/{e152}",
                     "primary144_effect_size": f"proportion={p144/e144:.6f}" if e144 else "",
                     "sensitivity152_effect_size": f"proportion={p152/e152:.6f}" if e152 else "",
                     "primary144_p": "NA", "sensitivity152_p": "NA", "primary144_q": "NA", "sensitivity152_q": "NA",
                     "primary144_denominator": f"n_total={t144}, n_evaluable={e144}",
                     "sensitivity152_denominator": f"n_total={t152}, n_evaluable={e152}",
                     "primary144_status": "DESCRIPTIVE", "sensitivity152_status": "DESCRIPTIVE",
                     "bh_family_id": "not_applicable_descriptive_prevalence", "bh_family_size": "prevalence rows are not hypothesis tests",
                     "method_change": "none: same evaluable-matrix definitions (Perfect/Strict; tet(M) four-state) in both scopes",
                     "grouping_change": GRP_CHANGE,
                     "note": "boundary8 contributes only to the sensitivity152 reference-assignment column"})

# AMR association Fisher: primary144 only
gap("D1_AMR", "AMR_feature_x_group", "Fisher association", "primary144 target/rest",
    "AMR Fisher association was run on primary144 only (family AMR_prevalence, 36 tests, 9 valid). A sensitivity152-reference-assignment Fisher run does not exist; boundary8 are excluded from target/rest by design, so no same-method pair can be formed.")

# ---------------- Module D2: defense burden KW ----------------
dstats = read_tsv(f"{S}/08_defense_mge_vfdb/defense_stats.tsv")
kw = {(r["analysis_id"], r["burden"]): r for r in dstats}
for burden in ("n_complete_systems", "n_complete_families", "n_partial_systems"):
    a = kw[("primary144", burden)]
    b = kw[("sensitivity152", burden)]
    rows.append({"module": "D2_defense", "feature_or_item": burden, "statistic": "Kruskal-Wallis", "design": "3 groups",
                 "paired": "YES",
                 "primary144_statistic": f"H={a['H']}", "sensitivity152_statistic": f"H={b['H']}",
                 "primary144_effect_size": f"epsilon2={a['epsilon2']}", "sensitivity152_effect_size": f"epsilon2={b['epsilon2']}",
                 "primary144_p": a["p"], "sensitivity152_p": b["p"],
                 "primary144_q": a["q_BH_within_3_burdens"], "sensitivity152_q": b["q_BH_within_3_burdens"],
                 "primary144_denominator": f"n={a['n']}, k={a['k_groups']}", "sensitivity152_denominator": f"n={b['n']}, k={b['k_groups']}",
                 "primary144_status": a["status"], "sensitivity152_status": b["status"],
                 "bh_family_id": "defense_burden_KW_omnibus", "bh_family_size": "3 burdens within scope",
                 "method_change": "none: identical R stats::kruskal.test in both scopes", "grouping_change": GRP_CHANGE, "note": ""})

gap("D2_defense", "defense_family_x_group", "Fisher association", "primary144 target/rest",
    "defense_family Fisher association was run on primary144 only (family defense_prevalence, 243 rows, 111 tested). No sensitivity152-reference-assignment run exists; boundary8 excluded from target/rest by design.")
gap("D2_defense", "vfdb_feature_x_group", "Fisher association", "primary144 target/rest",
    "VFDB Fisher association was run on primary144 only (family VFDB_prevalence, 57 rows, 42 tested). No sensitivity152-reference-assignment run exists.")
gap("D2_defense", "assembly_covariate_model", "rank regression", "primary144",
    "The assembly-covariate rank model was fitted on 144 complete cases only; no sensitivity152 fit exists (participant subset n=95 is a separate sensitivity, not the reference-assignment scope).")

# ---------------- Module C: association (GAP) ----------------
gap("C_function", "gene_species_association", "Fisher + BH (gene_species_main)", "primary144 target/rest",
    "Gene association was run on primary144 only (9,495 tests over 3,165 eligible clusters x 3 targets; candidates 36/26/119). A sensitivity152-reference run was not part of the validated scope: the eligible-cluster filter uses carrier144, so no same-method pair exists.")

# ---------------- Module E: interval null (GAP) ----------------
gap("E_intervals", "mge_positional_null", "permutation null (P=(1+b)/1001, joint BH)", "pooled152 + three primary groups",
    "The positional null was computed for pooled152 and the three primary144 groups (subflava recomputed at n=23, same RNG namespace main152_m8_20261015). A subflava-31 sensitivity152-reference null was not run, so no 144-vs-152 same-method pair exists; pooled152 already contains all boundary8 genomes.")

# ---------------- PCoA (GAP) ----------------
gap("B_structure", "pcoa_overall152", "wcmdscale Lingoes coordinates", "overall152 only",
    "PCoA is an overall152 ordination interface (formal merged files pcoa_coordinates_overall152.tsv / pcoa_eigenvalues_overall152.tsv); there is no primary144 plotting interface to pair with.")

write_tsv(f"{F}/12_delivery/sensitivity152_vs_primary144.tsv", FIELDS, rows)
print(f"contrast rows: {len(rows)} "
      f"(paired={sum(1 for r in rows if r['paired']=='YES')}, gap={sum(1 for r in rows if r['paired']=='NO_GAP')})")
