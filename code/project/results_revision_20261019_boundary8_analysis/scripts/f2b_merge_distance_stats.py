#!/usr/bin/env python3
"""Replace PERMDISP rows with residual-permutation results; retain other tests.

The input rows are archived before recalculation for comparison.
"""
from pathlib import Path
import csv, os

ROOT = str(Path(os.environ["NEISSERIA_PROJECT_ROOT"]))
S = f"{ROOT}/results_revision_20261019_boundary8_analysis"
F = f"{ROOT}/results_revision_20261019_boundary8_permdisp_fix"
OUT = f"{F}/05_structure"

def read_tsv(path):
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))

def write_tsv(path, fields, rows):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

old = read_tsv(f"{S}/05_structure/distance_statistics.tsv")
fields = list(old[0].keys())
new_pd = {r["analysis_id"]: r for r in read_tsv(f"{OUT}/permdisp_new_rows.tsv")}
legacy = [r for r in old if r["test"] == "PERMDISP"]
merged = []
for r in old:
    if r["test"] == "PERMDISP":
        n = new_pd[r["analysis_id"]]
        merged.append({k: n.get(k, "") for k in fields})
    else:
        merged.append(r)
write_tsv(f"{OUT}/distance_statistics.tsv", fields, merged)
write_tsv(f"{OUT}/legacy/permdisp_before_recalculation.tsv", fields, legacy)
print(f"merged {len(merged)} rows ({len(legacy)} PERMDISP replaced); input rows archived")
