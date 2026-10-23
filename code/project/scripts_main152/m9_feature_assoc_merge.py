#!/usr/bin/env python3
"""Merge AMR / defense / VFDB feature association tests into
08_defense_mge_vfdb/feature_association_all.tsv (contract F6_FEATURE_TESTS).

Required columns: family_id;feature;target_group;a;b;c;d;OR;CI_low;CI_high;p;q;status;reason
"""
from pathlib import Path
import csv, sys, pathlib

ROOT = pathlib.Path("results_revision_20261015_main152_species")
SOURCES = [
    ROOT / "07_amr_tetm" / "amr_association_all.tsv",
    ROOT / "08_defense_mge_vfdb" / "defense_family_tests.tsv",
    ROOT / "08_defense_mge_vfdb" / "vfdb_feature_tests.tsv",
]
COLS = ["family_id", "feature", "target_group", "a", "b", "c", "d",
        "OR", "CI_low", "CI_high", "p", "q", "status", "reason"]
OUT = ROOT / "08_defense_mge_vfdb" / "feature_association_all.tsv"

rows = []
for src in SOURCES:
    with open(src, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            missing = [c for c in COLS if c not in r]
            if missing:
                sys.exit(f"{src}: missing columns {missing}")
            rows.append({c: (r[c] if r[c] not in (None, "") else "NA") for c in COLS})

fams = sorted({r["family_id"] for r in rows})
print(f"families: {fams}")
assert len(fams) == 3, "expected exactly 3 family_ids"
with open(OUT, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=COLS, delimiter="\t", lineterminator="\n")
    w.writeheader()
    w.writerows(rows)
print(f"wrote {OUT} rows={len(rows)}")
