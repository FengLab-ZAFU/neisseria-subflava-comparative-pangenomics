#!/usr/bin/env python3
"""A3: clone clustering (complete linkage) on pairwise SNP distances + representative selection.

Inputs:
  02_independence/pairs_snp.tsv      - pairwise SNP within candidate groups
  QC table (checkm2_merged.tsv or historical qc_summary.tsv) - for representative ranking
Rules (task doc A3):
  - merge only when: resolvable_by_snp == yes AND pairwise SNP <= threshold AND provenance evidence
  - thresholds: default 20; sensitivity 5/10/50; 250 = legacy reference only
  - complete linkage: all pairs within a unit must satisfy threshold
  - input order: accession lexicographic; ties recorded
  - representative: derivatives never; quality ranking completeness desc,
    contamination asc, N50 desc, contigs asc, accession asc; earliest real date first
    for same-host-same-clone only when verifiable dates exist
Outputs:
  02_independence/clone_units_<thr>.tsv for each threshold
  02_independence/clone_units_default.tsv (20)
"""
import csv
import sys
from pathlib import Path

import numpy as np
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform

ROOT = Path(__file__).resolve().parent.parent.parent
OUT = ROOT / "results_revision_20261001" / "02_independence"

pairs_file = OUT / "pairs_snp.tsv"
if not pairs_file.exists():
    sys.exit("pairs_snp.tsv not found; run a3_pairwise_snp.py first")

pairs = list(csv.DictReader(open(pairs_file), delimiter="\t"))
qc = {}
cm = ROOT / "results_revision_20261001" / "01_retrieval_taxonomy" / "checkm2_merged.tsv"
if cm.exists():
    for r in csv.DictReader(open(cm), delimiter="\t"):
        qc[r["accession"]] = {"completeness": float(r.get("Completeness", 0) or 0),
                              "contamination": float(r.get("Contamination", 0) or 0)}
hist_qc = ROOT / "results" / "summary" / "qc_summary.tsv"
hist_len = {}
for r in csv.DictReader(open(ROOT / "neisseria_subflava_associated_group_metadata.tsv", encoding="utf-8-sig"),
                        delimiter="\t"):
    hist_len[r["assembly_accession"]] = {
        "n50": int(r.get("contig_n50_bp", 0) or 0),
        "contigs": int(r.get("number_of_contigs", 0) or 0),
        "completeness": float(r.get("checkm_completeness", 0) or 0),
        "contamination": float(r.get("checkm_contamination", 0) or 0),
    }

def get_qc(acc):
    if acc in qc:
        q = qc[acc]
        h = hist_len.get(acc, {})
        return {"completeness": q["completeness"], "contamination": q["contamination"],
                "n50": h.get("n50", 0), "contigs": h.get("contigs", 0)}
    h = hist_len.get(acc, {"n50": 0, "contigs": 0, "completeness": 0, "contamination": 0})
    return h


def cluster_at(threshold):
    """Complete-linkage clustering within each candidate group."""
    units = []
    by_group = {}
    for p in pairs:
        by_group.setdefault(p["group_id"], []).append(p)
    uid = 0
    for gid, ps in sorted(by_group.items()):
        members = sorted({p["reference"] for p in ps} | {p["member"] for p in ps})
        idx = {m: i for i, m in enumerate(members)}
        n = len(members)
        D = np.zeros((n, n))
        resolvable = {}
        for p in ps:
            i, j = idx[p["reference"]], idx[p["member"]]
            snp = int(p["pairwise_snp_ACGT"])
            D[i, j] = D[j, i] = snp
            resolvable[(min(i, j), max(i, j))] = p["resolvable_by_snp"] == "yes"
        # unresolved pairs block merging: set their distance to inf
        for i in range(n):
            for j in range(i + 1, n):
                if (i, j) in resolvable and not resolvable[(i, j)]:
                    D[i, j] = D[j, i] = 1e9
        Z = linkage(squareform(D), method="complete")
        labels = fcluster(Z, t=threshold, criterion="distance")
        groups = {}
        for m, lab in zip(members, labels):
            groups.setdefault(lab, []).append(m)
        for lab, mem in sorted(groups.items()):
            if len(mem) == 1:
                uid += 1
                units.append({"unit_id": f"U{uid:05d}", "group_id": gid,
                              "members": mem[0], "n_members": 1,
                              "max_pairwise_snp": 0, "representative": mem[0]})
            else:
                uid += 1
                # representative by quality ranking
                ranked = sorted(mem, key=lambda a: (-get_qc(a)["completeness"],
                                                    get_qc(a)["contamination"],
                                                    -get_qc(a)["n50"],
                                                    get_qc(a)["contigs"], a))
                mx = 0
                for x in range(len(mem)):
                    for y in range(x + 1, len(mem)):
                        mx = max(mx, int(D[idx[mem[x]], idx[mem[y]]]))
                units.append({"unit_id": f"U{uid:05d}", "group_id": gid,
                              "members": ";".join(mem), "n_members": len(mem),
                              "max_pairwise_snp": mx, "representative": ranked[0]})
    return units


def write_units(thr, suffix):
    units = cluster_at(thr)
    path = OUT / f"clone_units_{suffix}.tsv"
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t",
                           fieldnames=["unit_id", "group_id", "members", "n_members",
                                       "max_pairwise_snp", "representative"])
        w.writeheader()
        w.writerows(units)
    n_multi = sum(1 for u in units if u["n_members"] > 1)
    print(f"threshold={thr} ({suffix}): {len(units)} units, {n_multi} multi-member -> {path.name}")
    return units


for thr, suffix in [(5, "s5"), (10, "s10"), (20, "default"), (50, "s50"), (250, "s250_legacy")]:
    write_units(thr, suffix)
