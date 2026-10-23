#!/usr/bin/env python3
"""A3: build candidate duplicate groups for local coordinate alignment.

Group types:
  provenance   - same study_namespace + same non-unknown participant_id
  culture_synonym - identical/normalized strain name across records (ATCC/LMG/CCUG/NCTC)
  sg0010       - SG0010 baselines explicitly paired (repeat sequencing)

Output: 02_independence/candidate_groups.tsv (group_id, type, members, reason)
"""
import csv
import re
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent.parent
OUT = ROOT / "results_revision_20261001" / "02_independence"


def norm_strain(s):
    s = (s or "").strip()
    if not s:
        return ""
    s = s.upper().replace("STRAIN", "").replace(":", " ").strip()
    s = re.sub(r"\s+", " ", s)
    return s


rows = list(csv.DictReader(open(OUT / "provenance.tsv"), delimiter="\t"))
groups = []
gid = 0

# provenance groups
by_pp = defaultdict(list)
for r in rows:
    if r["participant_id"] not in ("unknown", "") and r["participant_confidence"] == "high" \
       and not r["participant_id"].startswith("not_applicable"):
        by_pp[(r["study_namespace"], r["participant_id"])].append(r)
for (study, part), rs in sorted(by_pp.items()):
    if len(rs) > 1:
        gid += 1
        groups.append((f"G{gid:04d}", "provenance", ";".join(x["assembly_accession"] for x in rs),
                       f"{study}:{part}"))

# culture synonyms
by_strain = defaultdict(list)
for r in rows:
    ns = norm_strain(r["strain"])
    if ns and any(tok in ns for tok in ("ATCC", "LMG", "CCUG", "NCTC", "DSM")):
        by_strain[ns].append(r)
for ns, rs in sorted(by_strain.items()):
    accs = {x["assembly_accession"] for x in rs}
    if len(accs) > 1:
        gid += 1
        groups.append((f"G{gid:04d}", "culture_synonym", ";".join(sorted(accs)), f"strain:{ns}"))

# SG0010 baselines
sg = [r["assembly_accession"] for r in rows if r["derivative_flag"] == "natural_baseline_repeat"]
if len(sg) > 1:
    gid += 1
    groups.append((f"G{gid:04d}", "sg0010_repeat", ";".join(sorted(sg)), "SG0010 baseline resequencing"))

with open(OUT / "candidate_groups.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["group_id", "group_type", "members", "reason"])
    w.writerows(groups)
print(f"wrote {len(groups)} candidate groups -> {OUT/'candidate_groups.tsv'}")
for g in groups:
    print(" ", g[0], g[1], len(g[2].split(";")), g[3][:60])
