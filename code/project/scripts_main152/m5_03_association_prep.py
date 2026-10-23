#!/usr/bin/env python3
"""M5a prep: one-versus-rest 2x2 tables for all eligible Shell+Cloud clusters
(carrier>=5) x 3 target groups. Writes fisher input and the filter ledger.
"""
from pathlib import Path
import os
import csv, json, sys

ROOT = str(Path(os.environ["NEISSERIA_PROJECT_ROOT"]))
RES = f"{ROOT}/results_revision_20261015_main152_species"
PGAP2 = f"{RES}/02_pangenome/pgap2_main152"
OUT = f"{RES}/06_function"
META = f"{RES}/00_freeze/input_packet_snapshot/00_input/cohort/main152_metadata.tsv"

GROUPS = ["subflava", "perflava", "flavescens"]

def main():
    group_of = {}
    with open(META) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            group_of[r["strain"]] = r["group_for_testing"]
    sizes = {g: sum(1 for v in group_of.values() if v == g) for g in GROUPS}
    assert sizes == {"subflava": 31, "perflava": 114, "flavescens": 7}, sizes

    partition = {}
    with open(f"{PGAP2}/pav_tables/partition_table.tsv") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            partition[r["cluster"]] = r["partition_recomputed"]

    rows = []          # eligible test rows
    ledger = []        # all clusters w/ filter disposition
    with open(f"{PGAP2}/pav_tables/pav_binary.tsv") as fh:
        rdr = csv.reader(fh, delimiter="\t")
        header = next(rdr)
        accessions = header[1:]
        gidx = {g: [i for i, a in enumerate(accessions) if group_of[a] == g] for g in GROUPS}
        for f in rdr:
            cl = f[0]
            vals = [int(x) for x in f[1:]]
            carrier = sum(vals)
            part = partition[cl]
            counts = {g: sum(vals[i] for i in gidx[g]) for g in GROUPS}
            if part not in ("Shell", "Cloud"):
                ledger.append((cl, part, carrier, "not_tested", "partition_not_shell_cloud"))
                continue
            if carrier < 5:
                ledger.append((cl, part, carrier, "not_tested", "carrier_lt_5"))
                continue
            ledger.append((cl, part, carrier, "tested", "eligible_shell_cloud_carrier_ge5"))
            for tg in GROUPS:
                rest = [g for g in GROUPS if g != tg]
                a = counts[tg]
                b = sizes[tg] - a
                c = sum(counts[g] for g in rest)
                d = sum(sizes[g] for g in rest) - c
                rows.append({
                    "cluster": cl, "target_group": tg, "a": a, "b": b, "c": c, "d": d,
                    "target_n": sizes[tg], "rest_n": sum(sizes[g] for g in rest),
                    "target_frequency": a / sizes[tg],
                    "rest_frequency": c / sum(sizes[g] for g in rest),
                    "freq_subflava": counts["subflava"] / sizes["subflava"],
                    "freq_perflava": counts["perflava"] / sizes["perflava"],
                    "freq_flavescens": counts["flavescens"] / sizes["flavescens"],
                    "partition": part, "carrier": carrier,
                })

    with open(f"{OUT}/fisher_input_association.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["job_id", "a", "b", "c", "d"])
        for i, r in enumerate(rows):
            w.writerow([f"assoc_{i}", r["a"], r["b"], r["c"], r["d"]])
    with open(f"{OUT}/association_rows_prep.tsv", "w", newline="") as fh:
        cols = ["cluster", "target_group", "a", "b", "c", "d", "target_n", "rest_n",
                "target_frequency", "rest_frequency", "prevalence_difference",
                "freq_subflava", "freq_perflava", "freq_flavescens", "partition", "carrier"]
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=cols, lineterminator="\n")
        w.writeheader()
        for r in rows:
            r["prevalence_difference"] = r["target_frequency"] - r["rest_frequency"]
            w.writerow(r)
    with open(f"{OUT}/cluster_filter_ledger.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["cluster", "partition", "carrier", "disposition", "reason"])
        w.writerows(ledger)

    n_clusters = len({r["cluster"] for r in rows})
    summary = {"group_sizes": sizes, "n_eligible_clusters": n_clusters,
               "n_tests": len(rows),
               "ledger_counts": {}}
    for _cl, _p, _c, disp, reason in ledger:
        summary["ledger_counts"][reason] = summary["ledger_counts"].get(reason, 0) + 1
    with open(f"{OUT}/association_prep_summary.json", "w") as fh:
        json.dump(summary, fh, indent=1)
    sys.stderr.write(json.dumps(summary) + "\n")

if __name__ == "__main__":
    main()
