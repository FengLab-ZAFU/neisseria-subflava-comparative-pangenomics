#!/usr/bin/env python3
"""Build canonical units using source identity and participant-level complete linkage.

Difference from fix3 h4: source-identity groups are atomic NODES; the
participant SNP rule then clusters NODES (identity groups and/or singleton
assemblies) within each verified (namespace, participant). Node-node distance
= max over cross-member assembly pairs (complete linkage; any unresolvable
pair = INF, blocks the merge). This way a Rochester participant with several
isolates, each isolate submitted twice, still collapses to one unit when the
isolates are <=threshold SNP apart, while the duplicate-submission pairs keep
their source-identity basis (SNP reported, never claimed <= threshold).

Representative: identity-alone units -> QC rank; participant_clone units ->
earliest TRUE collection date then QC rank (dates are year-precision for
Table-S1-verified Rochester records).

Outputs: results_revision_20261010_sourcefix_v7/03_units_v7/
"""
import os
import csv
import re
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform

ROOT = Path(os.environ["NEISSERIA_PROJECT_ROOT"])
OLD = ROOT / "results_revision_20261001"
F2 = ROOT / "results_revision_20261008_fix2"
F3 = ROOT / "results_revision_20261009_fix3"
SUP = ROOT / "results_revision_20261010_sourcefix_v7"
OUT = SUP / "03_units_v7"
OUT.mkdir(parents=True, exist_ok=True)
INF = 1e12
THRESHOLDS = [5, 10, 20, 50, 250]

# ---------- inputs ----------
prov = {r["assembly_accession"]: r for r in csv.DictReader(
    open(SUP / "inputs_v7/provenance_corrected_v7.tsv"), delimiter="\t")}
pairs = list(csv.DictReader(open(F2 / "01_canonical_units/pairs_snp_completed.tsv"), delimiter="\t"))
dist = {}
for p in pairs:
    if p["snp_status"] in ("OK", "EMPTY_OK") and p["resolvable_by_snp"] == "yes":
        d = int(p["pairwise_snp_ACGT"])
        k = frozenset((p["acc_a"], p["acc_b"]))
        if k not in dist or d < dist[k]:
            dist[k] = d

universe = set()
for p in pairs:
    universe.update((p["acc_a"], p["acc_b"]))
member_v4, member_exp = {}, {}
for r in csv.DictReader(open(F3 / "01_reference_v4/rescreen_v4.tsv"), delimiter="\t"):
    member_v4[r["accession"]] = r["associated_group_member_v4"] == "True"
    member_exp[r["accession"]] = r["associated_group_member_expanded"] == "True"
    if member_exp[r["accession"]]:
        universe.add(r["accession"])
for f in (F2 / "03_cohorts").glob("cohort_*.tsv"):
    for r in csv.DictReader(open(f), delimiter="\t"):
        universe.add(r["accession"])

igroups = {}
acc2ig = {}
for g in csv.DictReader(open(SUP / "inputs_v7/same_strain_identity_groups_v7.tsv"), delimiter="\t"):
    if g["decision"] == "merged":
        mem = sorted(g["members"].split(";"))
        igroups[g["group_id"]] = mem
        for a in mem:
            acc2ig[a] = g["group_id"]
            universe.add(a)

# ---------- QC / dates ----------
qc = {}
for r in csv.DictReader(open(OLD / "01_retrieval_taxonomy/checkm2_merged.tsv"), delimiter="\t"):
    qc[r["accession"]] = (float(r.get("Completeness", 0) or 0), float(r.get("Contamination", 0) or 0))
hist = {r["assembly_accession"]: r for r in csv.DictReader(
    open(ROOT / "neisseria_subflava_associated_group_metadata.tsv", encoding="utf-8-sig"), delimiter="\t")}
def qcm(acc):
    if acc in qc:
        return qc[acc]
    h = hist.get(acc, {})
    return (float(h.get("checkm_completeness", 0) or 0), float(h.get("checkm_contamination", 0) or 0))
n50, nctg = {}, {}
for r in csv.DictReader(open(F2 / "02_assembly_covariates/assembly_covariates_fasta.tsv"), delimiter="\t"):
    if r["n_contigs"] != "NA":
        n50[r["accession"]] = int(r["n50_bp"])
        nctg[r["accession"]] = int(r["n_contigs"])

def parse_date(s):
    if not s:
        return None
    s = s.strip().lower()
    if s in ("", "unknown", "na", "n/a", "none", "missing", "not applicable",
             "not provided", "not collected", "not available", "not determined"):
        return None
    m = re.search(r"(\d{4})(?:[-/](\d{1,2}))?(?:[-/](\d{1,2}))?", s)
    if not m:
        return None
    return (int(m.group(1)), int(m.group(2) or 0), int(m.group(3) or 0))

def fmt_date(t):
    if not t:
        return "NA"
    y, m, d = t
    if m == 0:
        return f"{y}"
    if d == 0:
        return f"{y}-{m:02d}"
    return f"{y}-{m:02d}-{d:02d}"

def col_date(acc):
    return parse_date(hist.get(acc, {}).get("collection_date", "")) or \
           parse_date(prov.get(acc, {}).get("collection_date", ""))

ELIGIBLE_FLAGS = {"natural", "natural_baseline_repeat"}
def source_eligible(acc):
    return prov.get(acc, {}).get("derivative_flag") in ELIGIBLE_FLAGS

def qc_pass(acc):
    c, x = qcm(acc)
    if c == 0 and acc not in qc and not hist.get(acc, {}).get("checkm_completeness"):
        return False  # QC unevaluable
    return c >= 90.0 and x <= 5.0 and n50.get(acc, 0) >= 10000

def membership_pass(acc):
    return member_v4.get(acc, False)

def eligible_for_cohort(acc):
    return source_eligible(acc) and qc_pass(acc) and membership_pass(acc)

def qc_rank(acc):
    c, x = qcm(acc)
    return (-c, x, -n50.get(acc, 0), nctg.get(acc, 99999), acc)

def rep_pick(mems, rule):
    elig = [a for a in mems if source_eligible(a)]
    pool = elig if elig else sorted(mems)
    if rule == "identity" or len(mems) == 1:
        ranked = sorted(pool, key=qc_rank)
    else:
        ranked = sorted(pool, key=lambda a: (col_date(a) or (9999, 0, 0), qc_rank(a)))
    return ranked[0], ranked

# ---------- nodes ----------
# node id: identity group id, or accession for singletons
def node_of(acc):
    return acc2ig.get(acc, acc)

nodes = {}  # nid -> members
for gid, mem in igroups.items():
    nodes[gid] = sorted(mem)

def node_participant(nid):
    """(ns, pid) if every member is identifiable with the SAME participant; else None."""
    mem = nodes.get(nid, [nid])
    keys = set()
    for a in mem:
        p = prov.get(a, {})
        if p.get("participant_status") != "identifiable":
            return None
        keys.add((p.get("study_namespace", ""), p.get("participant_id", "")))
    return keys.pop() if len(keys) == 1 else None

def node_dist(n1, n2):
    """max over cross-member pairs; INF if any unresolvable."""
    m1, m2 = nodes.get(n1, [n1]), nodes.get(n2, [n2])
    worst = 0
    for a in m1:
        for b in m2:
            d = dist.get(frozenset((a, b)))
            if d is None:
                return INF
            worst = max(worst, d)
    return worst

# participant groups over nodes
pnodes = defaultdict(set)
for acc, p in prov.items():
    if p.get("participant_status") == "identifiable":
        nid = node_of(acc)
        key = node_participant(nid)
        if key:
            pnodes[key].add(nid)

def cluster_nodes(nids, threshold):
    nids = sorted(nids)
    if len(nids) == 1:
        return [nids]
    idx = {n: i for i, n in enumerate(nids)}
    D = np.zeros((len(nids), len(nids)))
    for i in range(len(nids)):
        for j in range(i + 1, len(nids)):
            D[i, j] = D[j, i] = node_dist(nids[i], nids[j])
    Z = linkage(squareform(D, checks=False), method="complete")
    lab = fcluster(Z, t=threshold, criterion="distance")
    cl = defaultdict(list)
    for n, l in zip(nids, lab):
        cl[l].append(n)
    return [sorted(v) for v in cl.values()]

# ---------- build units ----------
units_by_t = {t: [] for t in THRESHOLDS}
clustered_nodes = set()
for key, nids in sorted(pnodes.items()):
    per_t = {t: cluster_nodes(nids, t) for t in THRESHOLDS}
    for t in THRESHOLDS:
        for c in per_t[t]:
            mems = sorted({a for n in c for a in nodes.get(n, [n])})
            if len(c) == 1 and c[0] in igroups:
                basis, group = "identity", c[0]
            elif len(c) == 1:
                basis, group = "singleton", ""
            else:
                basis = "participant_clone"
                group = f"{key[0]}:{key[1]}"
                if any(n in igroups for n in c):
                    group += " (with source-identity duplicate subgroups)"
            units_by_t[t].append({"members": mems, "basis": basis, "group": group})
    clustered_nodes.update(nids)
# identity groups not in any participant group (members not identifiable)
for gid, mem in sorted(igroups.items()):
    if gid in clustered_nodes:
        continue
    for t in THRESHOLDS:
        units_by_t[t].append({"members": sorted(mem), "basis": "identity", "group": gid})
# everything else in universe
assigned = {a for t in (20,) for u in units_by_t[t] for a in u["members"]}
for a in sorted(universe - assigned):
    for t in THRESHOLDS:
        units_by_t[t].append({"members": [a], "basis": "singleton", "group": ""})
# bookkeeping: every screened provenance record keeps a singleton unit even if it
# is outside the analysis universe (non-member, no pairs) - never silently dropped
assigned = {a for u in units_by_t[20] for a in u["members"]}
for a in sorted(set(prov) - assigned):
    for t in THRESHOLDS:
        units_by_t[t].append({"members": [a], "basis": "singleton", "group": ""})

# ---------- write unit tables + representative ledger (s20) ----------
rep_ledger = []
unit_tables = {}
a2u20 = {}
for t in THRESHOLDS:
    us = sorted(units_by_t[t], key=lambda u: min(u["members"]))
    rows = []
    for i, u in enumerate(us, 1):
        uid = f"CU{i:05d}"
        mem = u["members"]
        if t == 20:
            rule = "identity" if u["basis"] == "identity" else "participant"
            rep, ranked = rep_pick(mem, rule)
            rrule = {"identity": "qc_rank_same_isolate",
                     "participant": "earliest_date_then_qc"}[rule] if len(mem) > 1 else "singleton"
            for a in mem:
                rep_ledger.append({
                    "unit_id": uid, "accession": a, "n_members": len(mem),
                    "unit_basis": u["basis"], "unit_group": u["group"],
                    "source_eligible": source_eligible(a),
                    "qc_pass": qc_pass(a),
                    "membership_pass": membership_pass(a),
                    "eligible_for_cohort": eligible_for_cohort(a),
                    "completeness": qcm(a)[0], "contamination": qcm(a)[1],
                    "n50_bp": n50.get(a, "NA"), "n_contigs": nctg.get(a, "NA"),
                    "collection_date": fmt_date(col_date(a)),
                    "chosen_representative": rep, "is_representative": a == rep,
                    "rep_rule": rrule,
                    "rank_position": (ranked.index(a) + 1) if a in ranked else "NA"})
            a2u20.update({a: uid for a in mem})
        maxd = "NA"
        ds = [dist.get(frozenset((mem[x], mem[y]))) for x in range(len(mem)) for y in range(x + 1, len(mem))]
        have = [d for d in ds if d is not None]
        if have:
            maxd = max(have)
        rows.append({"unit_id": uid, "members": ";".join(mem), "n_members": len(mem),
                     "unit_basis": u["basis"], "unit_group": u["group"],
                     "max_pairwise_snp": maxd,
                     "representative": rep if t == 20 else "",
                     "rep_rule": rrule if t == 20 else ""})
    unit_tables[t] = rows
    with open(OUT / f"canonical_units_s{t}_v7.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

with open(OUT / "accession_to_unit_v7.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["accession", "unit_id"])
    for a in sorted(a2u20):
        w.writerow([a, a2u20[a]])
with open(OUT / "representative_ledger_v7.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, delimiter="\t", fieldnames=list(rep_ledger[0].keys()))
    w.writeheader()
    w.writerows(rep_ledger)

# ---------- intra-unit pair validation (s20) ----------
pv = []
for r in unit_tables[20]:
    mem = r["members"].split(";")
    if len(mem) < 2:
        continue
    for x in range(len(mem)):
        for y in range(x + 1, len(mem)):
            a, b = mem[x], mem[y]
            d = dist.get(frozenset((a, b)))
            if acc2ig.get(a) and acc2ig.get(a) == acc2ig.get(b):
                mb = "source_identity"
                note = "source-identity pair (same isolate/library): SNP reported, merge justified by provenance not SNP"
            else:
                ok = (d is not None and d <= 20)
                mb = "snp_clone_threshold"
                note = "PASS" if ok else "FAIL: cross-isolate pair exceeds 20 SNP or unresolvable within unit"
            pv.append({"unit_id": r["unit_id"], "acc_a": a, "acc_b": b,
                       "unit_basis": r["unit_basis"], "unit_group": r["unit_group"],
                       "identity_group_pair": acc2ig.get(a) if acc2ig.get(a) == acc2ig.get(b) else "",
                       "pairwise_snp": d if d is not None else "UNRESOLVED",
                       "merge_basis": mb, "validation": note})
with open(OUT / "intraunit_pair_validation_v7.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, delimiter="\t", fieldnames=list(pv[0].keys()))
    w.writeheader()
    w.writerows(pv)
print("intra-unit pair validation FAILs:", sum(1 for r in pv if r["validation"].startswith("FAIL")))

# ---------- threshold sensitivity ----------
with open(OUT / "threshold_sensitivity_v7.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["threshold_snp", "n_units", "n_multi_member_units", "largest_unit",
                "member_occurrences", "unique_members", "note"])
    for t in THRESHOLDS:
        rows = unit_tables[t]
        occ = sum(r["n_members"] for r in rows)
        uniq = len({a for r in rows for a in r["members"].split(";")})
        note = ("historical audit reference only; no pre-merge at any threshold" if t == 250
                else "threshold applied over isolate nodes within verified participant groups")
        w.writerow([t, len(rows), sum(1 for r in rows if r["n_members"] > 1),
                    max(r["n_members"] for r in rows), occ, uniq, note])

# ---------- near-neighbour hints ----------
hints = []
for k, d in sorted(dist.items(), key=lambda kv: kv[1]):
    if d > 20:
        continue
    a, b = sorted(k)
    if a2u20.get(a) and a2u20.get(a) == a2u20.get(b):
        continue
    hints.append({"acc_a": a, "acc_b": b, "pairwise_snp": d,
                  "unit_a": a2u20.get(a, ""), "unit_b": a2u20.get(b, ""),
                  "participant_a": f"{prov.get(a,{}).get('study_namespace','')}:{prov.get(a,{}).get('participant_id','')}",
                  "participant_b": f"{prov.get(b,{}).get('study_namespace','')}:{prov.get(b,{}).get('participant_id','')}",
                  "note": "near neighbour: no same-strain/same-participant source evidence -> NOT folded"})
with open(OUT / "near_neighbor_hints_v7.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, delimiter="\t", fieldnames=(list(hints[0].keys()) if hints else
                       ["acc_a", "acc_b", "pairwise_snp", "unit_a", "unit_b",
                        "participant_a", "participant_b", "note"]))
    w.writeheader()
    w.writerows(hints)

# ---------- unit diff vs fix3 v4 ----------
v4u = {r["accession"]: r["unit_id"] for r in csv.DictReader(
    open(F3 / "03_units_v4/accession_to_unit_v4.tsv"), delimiter="\t")}
v4units = {u["unit_id"]: u for u in csv.DictReader(
    open(F3 / "03_units_v4/canonical_units_s20_v4.tsv"), delimiter="\t")}
v5mem_of = {r["unit_id"]: r["members"].split(";") for r in unit_tables[20]}
diff = []
for a in sorted(set(v4u) | set(a2u20)):
    old, new = v4u.get(a, "NOT_IN_V4"), a2u20.get(a, "NOT_IN_V6")
    if old == new:
        continue
    oldmem = sorted(v4units.get(old, {}).get("members", "").split(";")) if old in v4units else []
    newmem = v5mem_of.get(new, [])
    cat = "renamed_same_membership" if oldmem and sorted(oldmem) == sorted(newmem) else "membership_changed"
    diff.append({"accession": a, "unit_v4": old, "unit_v7": new, "change_category": cat})
with open(OUT / "unit_diff_v7_vs_fix3.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, delimiter="\t", fieldnames=["accession", "unit_v4", "unit_v7", "change_category"])
    w.writeheader()
    w.writerows(diff)

# ---------- unit diff vs ns-fix v6 (2026-10-10 sourcefix v7) ----------
V6 = ROOT / "results_revision_20261010_nsfix_v6"
v5u = {r["accession"]: r["unit_id"] for r in csv.DictReader(
    open(V6 / "03_units_v6/accession_to_unit_v6.tsv"), delimiter="\t")}
v5units = {u["unit_id"]: u for u in csv.DictReader(
    open(V6 / "03_units_v6/canonical_units_s20_v6.tsv"), delimiter="\t")}
v6mem_of = {r["unit_id"]: r["members"].split(";") for r in unit_tables[20]}
diff5 = []
for a in sorted(set(v5u) | set(a2u20)):
    old, new = v5u.get(a, "NOT_IN_V6"), a2u20.get(a, "NOT_IN_V7")
    oldmem = sorted(v5units.get(old, {}).get("members", "").split(";")) if old in v5units else []
    newmem = v6mem_of.get(new, [])
    if oldmem and sorted(oldmem) == sorted(newmem):
        continue  # same membership (possibly renamed)
    diff5.append({"accession": a, "unit_v6": old, "unit_v7": new,
                  "members_v6": ";".join(oldmem), "members_v7": ";".join(newmem),
                  "participant_v6": f"{prov.get(a,{}).get('study_namespace','')}:{prov.get(a,{}).get('participant_id','')}",
                  "change_category": "membership_changed" if old != "NOT_IN_V6" and new != "NOT_IN_V7"
                                     else "added_or_removed"})
with open(OUT / "unit_diff_v7_vs_v6.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, delimiter="\t",
                       fieldnames=["accession", "unit_v6", "unit_v7", "members_v6", "members_v7",
                                   "participant_v6", "change_category"])
    w.writeheader()
    w.writerows(diff5)

# ---------- summary ----------
for t in THRESHOLDS:
    rows = unit_tables[t]
    multi = [r for r in rows if r["n_members"] > 1]
    print(f"s{t}: units={len(rows)} multi={len(multi)} "
          f"identity={sum(1 for r in multi if r['unit_basis']=='identity')} "
          f"participant_clone={sum(1 for r in multi if r['unit_basis']=='participant_clone')}")
print("near-neighbour hints (<=20 SNP across units):", len(hints))
print("unit diff vs fix3:", dict(Counter(d["change_category"] for d in diff)),
      "| accessions changed:", len(diff))
print("unit diff vs v6 (membership changed):", len(diff5))
for d in diff5:
    print("   ", d["accession"], d["unit_v6"], "->", d["unit_v7"], d["participant_v6"])
