#!/usr/bin/env python3
"""M7c summarize: parse per-genome ABRicate/VFDB raw TSVs for the main152 cohort.

Outputs (under results_revision_20261015_main152_species/08_defense_mge_vfdb/):
  VFDB_presence_long.tsv   accession x candidate-feature long table
  vfdb/vfdb_run_status.tsv one row per genome (run status, n hits)
  vfdb_prevalence.tsv      feature x group prevalence with Wilson 95% CI
  vfdb_candidate_catalog.tsv per-feature catalog (VFDB id/product, carriers, group freqs)
  vfdb/vfdb_input_binding.tsv accession/path/sha256 ledger
  vfdb/vfdb_test_contingency.tsv input for m7c_vfdb_tests.R (a,b,c,d per feature x group)
"""
from pathlib import Path
import csv, math, os, sys

PROJECT_ROOT = str(Path(os.environ["NEISSERIA_PROJECT_ROOT"]))
RESULTS = os.path.join(PROJECT_ROOT, "results_revision_20261015_main152_species")
OUT = os.path.join(RESULTS, "08_defense_mge_vfdb")
RAW = os.path.join(OUT, "vfdb", "raw")
COHORT = os.path.join(RESULTS, "00_freeze/input_packet_snapshot/00_input/cohort")
BINDING = os.path.join(RESULTS, "00_preflight/input_check/input_binding.tsv")

DB_SNAPSHOT = ("vfdb_4592nucl_2026-Apr-03_"
               "sha256:e827dddb1841f7d31cd1bcf03dbfd45a0959648786ee70ae21a647475de15dce")
GROUPS = ["subflava", "perflava", "flavescens"]
MIN_CARRIER = 5

accessions = [l.strip() for l in open(os.path.join(COHORT, "main152_accessions.txt")) if l.strip()]
assert len(accessions) == 152 and len(set(accessions)) == 152

group_of = {}
with open(os.path.join(COHORT, "main152_metadata.tsv")) as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        group_of[row["strain"]] = row["group_for_testing"]
assert set(group_of) == set(accessions)

# input binding ledger (reuse verified hashes from M0)
binding_rows = []
with open(BINDING) as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        if row["accession"] in group_of:
            binding_rows.append(row)
assert len(binding_rows) == 152
with open(os.path.join(OUT, "vfdb", "vfdb_input_binding.tsv"), "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(["accession", "fasta_path", "expected_sha256", "actual_sha256", "m0_binding_status"])
    for r in sorted(binding_rows, key=lambda r: accessions.index(r["accession"])):
        w.writerow([r["accession"], r["path"], r["expected_sha256"], r["actual_sha256"], r["status"]])

# exit codes from the run script
exit_code = {}
with open(os.path.join(RAW, "run_exit_codes.tsv")) as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        exit_code[row["accession"]] = row["exit_code"]

# parse raw abricate output
# columns: #FILE SEQUENCE START END STRAND GENE COVERAGE COVERAGE_MAP GAPS %COVERAGE %IDENTITY DATABASE ACCESSION PRODUCT RESISTANCE
hits = {}          # acc -> {feature: {"ident": max, "cov": max, "acc": vfdb id, "prod": product}}
run_status = {}    # acc -> dict
for acc in accessions:
    ec = exit_code.get(acc, "NOT_RUN")
    raw_path = os.path.join(RAW, f"{acc}.tsv")
    acc_hits = {}
    if ec == "0" and os.path.isfile(raw_path):
        with open(raw_path) as fh:
            rd = csv.reader(fh, delimiter="\t")
            header = next(rd, None)
            for rec in rd:
                if len(rec) < 14:
                    continue
                gene, cov, ident, vfdb_acc, prod = rec[5], rec[9], rec[10], rec[12], rec[13]
                cur = acc_hits.setdefault(gene, {"ident": 0.0, "cov": 0.0, "acc": vfdb_acc, "prod": prod})
                try:
                    cur["ident"] = max(cur["ident"], float(ident))
                    cur["cov"] = max(cur["cov"], float(cov))
                except ValueError:
                    pass
    hits[acc] = acc_hits
    if ec == "0":
        st = "SUCCESS_HITS" if acc_hits else "SUCCESS_NO_HITS"
    elif ec == "INPUT_MISSING":
        st = "FAILED_INPUT_MISSING"
    elif ec == "NOT_RUN":
        st = "NOT_RUN"
    else:
        st = "FAILED_EXIT_" + ec
    run_status[acc] = {"accession": acc, "group": group_of[acc], "exit_code": ec,
                       "run_status": st, "n_hits": len(acc_hits)}

with open(os.path.join(OUT, "vfdb", "vfdb_run_status.tsv"), "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(["accession", "group", "exit_code", "run_status", "n_vfdb_hits", "raw_tsv"])
    for acc in accessions:
        r = run_status[acc]
        w.writerow([acc, r["group"], r["exit_code"], r["run_status"], r["n_hits"],
                    os.path.join(RAW, acc + ".tsv")])

evaluable = [a for a in accessions if run_status[a]["exit_code"] == "0"]
n_eval = len(evaluable)

# candidate feature universe = any feature detected in >=1 evaluable genome
feature_meta = {}  # feature -> {acc, prod, carriers:set}
for acc in evaluable:
    for feat, h in hits[acc].items():
        m = feature_meta.setdefault(feat, {"acc": h["acc"], "prod": h["prod"], "carriers": set()})
        m["carriers"].add(acc)
features = sorted(feature_meta)

# presence long table: 152 accessions x all candidate features
with open(os.path.join(OUT, "VFDB_presence_long.tsv"), "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(["accession", "feature", "detected", "run_status",
                "identity_percent", "coverage_percent", "DB_snapshot"])
    for acc in accessions:
        st = run_status[acc]["run_status"]
        ok = run_status[acc]["exit_code"] == "0"
        for feat in features:
            h = hits[acc].get(feat)
            if ok and h:
                w.writerow([acc, feat, 1, st, f'{h["ident"]:.3f}', f'{h["cov"]:.3f}', DB_SNAPSHOT])
            else:
                w.writerow([acc, feat, 0, st, "NA", "NA", DB_SNAPSHOT])

def wilson(k, n, z=1.959963984540054):
    if n == 0:
        return ("NA", "NA", "NA")
    p = k / n
    denom = 1 + z * z / n
    centre = p + z * z / (2 * n)
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    return (f"{p:.6f}", f"{max(0.0, (centre - half) / denom):.6f}",
            f"{min(1.0, (centre + half) / denom):.6f}")

# prevalence: feature x group
group_eval = {g: [a for a in evaluable if group_of[a] == g] for g in GROUPS}
group_total = {g: sum(1 for a in accessions if group_of[a] == g) for g in GROUPS}
with open(os.path.join(OUT, "vfdb_prevalence.tsv"), "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(["feature", "group", "n_total", "n_evaluable", "n_positive",
                "prevalence", "wilson95_low", "wilson95_high", "status"])
    for feat in features:
        carriers = feature_meta[feat]["carriers"]
        for g in GROUPS:
            ne = len(group_eval[g])
            k = sum(1 for a in group_eval[g] if a in carriers)
            p, lo, hi = wilson(k, ne)
            status = "OK" if ne > 0 else "NO_EVALUABLE_RUNS"
            w.writerow([feat, g, group_total[g], ne, k, p, lo, hi, status])

# candidate catalog: one row per feature
with open(os.path.join(OUT, "vfdb_candidate_catalog.tsv"), "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(["feature", "vfdb_accession", "vfdb_product", "n_genomes_positive",
                "n_evaluable_total", "overall_prevalence",
                "n_positive_subflava", "freq_subflava",
                "n_positive_perflava", "freq_perflava",
                "n_positive_flavescens", "freq_flavescens",
                "min_identity_observed", "min_coverage_observed"])
    for feat in sorted(features, key=lambda f: (-len(feature_meta[f]["carriers"]), f)):
        m = feature_meta[feat]
        carriers = m["carriers"]
        row = [feat, m["acc"], m["prod"], len(carriers), n_eval,
               f"{len(carriers)/n_eval:.6f}" if n_eval else "NA"]
        for g in GROUPS:
            ne = len(group_eval[g])
            k = sum(1 for a in group_eval[g] if a in carriers)
            row += [k, f"{k/ne:.6f}" if ne else "NA"]
        minid = min(hits[a][feat]["ident"] for a in carriers)
        mincov = min(hits[a][feat]["cov"] for a in carriers)
        row += [f"{minid:.3f}", f"{mincov:.3f}"]
        w.writerow(row)

# contingency table for Fisher tests (R step)
with open(os.path.join(OUT, "vfdb", "vfdb_test_contingency.tsv"), "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(["family_id", "feature", "target_group", "a", "b", "c", "d",
                "carrier_total", "n_evaluable_total", "status", "reason"])
    for feat in features:
        carriers = feature_meta[feat]["carriers"]
        npos = len(carriers)
        if npos < MIN_CARRIER:
            base_status, base_reason = "FILTERED_BELOW_CARRIER_MIN", f"carrier_total<{MIN_CARRIER}"
        elif npos == n_eval:
            base_status, base_reason = "NOT_ESTIMABLE", "no_variation_all_positive"
        else:
            base_status, base_reason = "TESTED", "NA"
        for g in GROUPS:
            rest = [x for x in GROUPS if x != g]
            a = sum(1 for x in group_eval[g] if x in carriers)
            b = len(group_eval[g]) - a
            c = sum(1 for rg in rest for x in group_eval[rg] if x in carriers)
            d = sum(len(group_eval[rg]) for rg in rest) - c
            w.writerow(["VFDB_prevalence", feat, g, a, b, c, d, npos, n_eval,
                        base_status, base_reason])

print(f"evaluable={n_eval}/152  features={len(features)}")
print("status counts:")
from collections import Counter
for k, v in sorted(Counter(run_status[a]["run_status"] for a in accessions).items()):
    print(f"  {k}: {v}")
