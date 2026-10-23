#!/usr/bin/env python3
"""fix3-H4 (Phase A2): canonical units v4 - no hidden 250-SNP pre-merge.

Merge domains (task 2026-10-09 sec A2):
  1. same-strain / same-culture IDENTITY GROUPS (evidence-based, h3): contracted
     into one unit at every threshold as SOURCE duplicates. Representative by
     QC rank. Intra-group SNP reported, never claimed to be <= threshold.
  2. same (study_namespace, participant_id) groups with participant_status =
     identifiable: complete linkage on the ACTUAL member pair matrix at the
     working threshold (5/10/20/50/250). Missing / unresolvable / low-coverage
     pairs = +inf: they block merges, never create them. No 250-SNP super-node
     pre-contraction anywhere; 250 is reported as a historical audit threshold.
  3. everything else: singleton units. Cross-group low-SNP pairs are recorded
     as near-neighbour hints only (no same-strain evidence -> no folding).

Representatives (task A2): eligibility first (natural / natural_baseline_repeat
only), then identity groups by QC rank; participant clone units by TRUE earliest
collection date then QC rank (sample ordinals like S001 are NOT dates);
ties broken deterministically by accession. tet(M)/PAM/significance never used.

Outputs (results_revision_20261009_fix3/03_units_v4/):
  canonical_units_s{5,10,20,50,250}_v4.tsv
  accession_to_unit_v4.tsv
  representative_ledger_v4.tsv
  intraunit_pair_validation_v4.tsv
  threshold_sensitivity_v4.tsv
  near_neighbor_hints_v4.tsv
  unit_diff_v4_vs_v3.tsv
"""
import os
import csv
import re
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform

ROOT = Path(os.environ["NEISSERIA_PROJECT_ROOT"])
OLD = ROOT / "results_revision_20261001"
F2 = ROOT / "results_revision_20261008_fix2"
F3 = ROOT / "results_revision_20261009_fix3"
OUT = F3 / "03_units_v4"
OUT.mkdir(parents=True, exist_ok=True)
INF = 1e12
THRESHOLDS = [5, 10, 20, 50, 250]

# ---------- inputs ----------
prov = {r["assembly_accession"]: r for r in csv.DictReader(
    open(F3 / "02_provenance_v4/provenance_corrected.tsv"), delimiter="\t")}
pairs = list(csv.DictReader(open(F2 / "01_canonical_units/pairs_snp_completed.tsv"), delimiter="\t"))
dist = {}
for p in pairs:
    if p["snp_status"] in ("OK", "EMPTY_OK") and p["resolvable_by_snp"] == "yes":
        d = int(p["pairwise_snp_ACGT"])
        k = frozenset((p["acc_a"], p["acc_b"]))
        if k not in dist or d < dist[k]:
            dist[k] = d

# universe: pair members + expanded member pool + v3 cohort members + identity groups
universe = set()
for p in pairs:
    universe.update((p["acc_a"], p["acc_b"]))
for r in csv.DictReader(open(F3 / "01_reference_v4/rescreen_v4.tsv"), delimiter="\t"):
    if r["associated_group_member_expanded"] == "True":
        universe.add(r["accession"])
for f in (F2 / "03_cohorts").glob("cohort_*.tsv"):
    for r in csv.DictReader(open(f), delimiter="\t"):
        universe.add(r["accession"])

igroups = {}   # gid -> members (merged only)
acc2ig = {}
for g in csv.DictReader(open(F3 / "02_provenance_v4/same_strain_identity_groups.tsv"), delimiter="\t"):
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
    """true collection date -> sortable tuple; None if not a real date."""
    if not s:
        return None
    s = s.strip().lower()
    if s in ("", "unknown", "na", "n/a", "none", "missing", "not applicable",
             "not provided", "not collected", "not available"):
        return None
    m = re.search(r"(\d{4})(?:[-/](\d{1,2}))?(?:[-/](\d{1,2}))?", s)
    if not m:
        return None
    return (int(m.group(1)), int(m.group(2) or 0), int(m.group(3) or 0))

def col_date(acc):
    return parse_date(hist.get(acc, {}).get("collection_date", "")) or \
           parse_date(prov.get(acc, {}).get("collection_date", ""))

ELIGIBLE_FLAGS = {"natural", "natural_baseline_repeat"}
def eligible(acc):
    return prov.get(acc, {}).get("derivative_flag") in ELIGIBLE_FLAGS

def qc_rank(acc):
    c, x = qcm(acc)
    return (-c, x, -n50.get(acc, 0), nctg.get(acc, 99999), acc)

def rep_pick(mems, rule):
    """rule: 'identity' -> QC rank; 'participant' -> earliest true date then QC."""
    elig = [a for a in mems if eligible(a)]
    pool = elig if elig else sorted(mems)  # fallback: no eligible member
    if rule == "identity" or len(mems) == 1:
        ranked = sorted(pool, key=qc_rank)
    else:
        ranked = sorted(pool, key=lambda a: (col_date(a) or (9999, 0, 0), qc_rank(a)))
    return ranked[0], ranked

# ---------- build units ----------
# super-nodes: identity groups; participant groups for identifiable participants
pgroups = defaultdict(list)
for acc, p in prov.items():
    if p["participant_status"] == "identifiable" and acc not in acc2ig:
        pgroups[(p["study_namespace"], p["participant_id"])].append(acc)
# participant groups may include accessions outside `universe`; intersect later
for gid, mem in igroups.items():
    for a in mem:
        assert prov.get(a, {}).get("participant_status") != "identifiable", \
            f"identity group member {a} unexpectedly identifiable"

def complete_linkage_clusters(mems, threshold):
    mems = sorted(mems)
    if len(mems) == 1:
        return [mems]
    idx = {a: i for i, a in enumerate(mems)}
    D = np.full((len(mems), len(mems)), INF)
    np.fill_diagonal(D, 0)
    for k, d in dist.items():
        a, b = tuple(k)
        if a in idx and b in idx:
            D[idx[a], idx[b]] = D[idx[b], idx[a]] = d
    Z = linkage(squareform(D, checks=False), method="complete")
    lab = fcluster(Z, t=threshold, criterion="distance")
    cl = defaultdict(list)
    for a, l in zip(mems, lab):
        cl[l].append(a)
    return [sorted(v) for v in cl.values()]

# unit objects: {"members": [...], "basis": "identity"|"participant_clone"|"singleton",
#                "group": gid or (ns,pid)}
units_by_t = {t: [] for t in THRESHOLDS}
assigned = set()
for gid, mem in sorted(igroups.items()):
    mem = [a for a in mem if a in universe or True]  # identity groups always enter universe
    for t in THRESHOLDS:
        units_by_t[t].append({"members": sorted(mem), "basis": "identity", "group": gid})
    assigned.update(mem)
for key, mem in sorted(pgroups.items()):
    mem = sorted(set(mem) - assigned)
    if not mem:
        continue
    if len(mem) == 1:
        for t in THRESHOLDS:
            units_by_t[t].append({"members": mem, "basis": "singleton", "group": ""})
        assigned.update(mem)
        continue
    per_t = {t: complete_linkage_clusters(mem, t) for t in THRESHOLDS}
    for t in THRESHOLDS:
        for c in per_t[t]:
            basis = "participant_clone" if len(c) > 1 else "singleton"
            units_by_t[t].append({"members": c, "basis": basis,
                                  "group": f"{key[0]}:{key[1]}" if len(c) > 1 else ""})
    assigned.update(mem)
for a in sorted(universe - assigned):
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
            rrule = {"identity": "qc_rank_same_culture",
                     "participant": "earliest_date_then_qc"}[rule] if len(mem) > 1 else "singleton"
            for a in mem:
                rep_ledger.append({
                    "unit_id": uid, "accession": a, "n_members": len(mem),
                    "unit_basis": u["basis"], "unit_group": u["group"],
                    "eligible_natural": eligible(a),
                    "completeness": qcm(a)[0], "contamination": qcm(a)[1],
                    "n50_bp": n50.get(a, "NA"), "n_contigs": nctg.get(a, "NA"),
                    "collection_date_parsed": "-".join(f"{x:02d}" for x in col_date(a)) if col_date(a) else "NA",
                    "chosen_representative": rep, "is_representative": a == rep,
                    "rep_rule": rrule,
                    "rank_position": (ranked.index(a) + 1) if a in ranked else "NA"})
            a2u = {a: uid for a in mem}
            a2u20.update(a2u)
        else:
            rep = ""
            rrule = ""
        maxd = "NA"
        ds = [dist.get(frozenset((mem[x], mem[y]))) for x in range(len(mem)) for y in range(x + 1, len(mem))]
        have = [d for d in ds if d is not None]
        if have:
            maxd = max(have)
        rows.append({"unit_id": uid, "members": ";".join(mem), "n_members": len(mem),
                     "unit_basis": u["basis"], "unit_group": u["group"],
                     "max_pairwise_snp": maxd, "representative": rep if t == 20 else "",
                     "rep_rule": rrule})
    unit_tables[t] = rows
    with open(OUT / f"canonical_units_s{t}_v4.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

with open(OUT / "accession_to_unit_v4.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["accession", "unit_id"])
    for a in sorted(a2u20):
        w.writerow([a, a2u20[a]])
with open(OUT / "representative_ledger_v4.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, delimiter="\t", fieldnames=list(rep_ledger[0].keys()))
    w.writeheader()
    w.writerows(rep_ledger)

# ---------- intra-unit pair validation (s20) ----------
pv = []
for r in unit_tables[20]:
    mem = r["members"].split(";")
    if len(mem) < 2:
        continue
    basis = r["unit_basis"]
    for x in range(len(mem)):
        for y in range(x + 1, len(mem)):
            a, b = mem[x], mem[y]
            d = dist.get(frozenset((a, b)))
            if basis == "participant_clone":
                ok = (d is not None and d <= 20)
                mb = "snp_clone_threshold"
                note = "PASS" if ok else "FAIL: participant clone pair exceeds 20 SNP or unresolvable"
            else:
                ok = True
                mb = "source_identity"
                note = "source-identity exception: merge justified by provenance, SNP reported not claimed <=20"
            pv.append({"unit_id": r["unit_id"], "acc_a": a, "acc_b": b,
                       "unit_basis": basis, "unit_group": r["unit_group"],
                       "pairwise_snp": d if d is not None else "UNRESOLVED",
                       "merge_basis": mb, "validation": note})
with open(OUT / "intraunit_pair_validation_v4.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, delimiter="\t", fieldnames=list(pv[0].keys()))
    w.writeheader()
    w.writerows(pv)
n_fail = sum(1 for r in pv if r["validation"].startswith("FAIL"))
print("intra-unit pair validation FAILs:", n_fail)

# ---------- threshold sensitivity ----------
with open(OUT / "threshold_sensitivity_v4.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["threshold_snp", "n_units", "n_multi_member_units", "largest_unit",
                "member_occurrences", "unique_members", "note"])
    for t in THRESHOLDS:
        rows = unit_tables[t]
        occ = sum(r["n_members"] for r in rows)
        uniq = len({a for r in rows for a in r["members"].split(";")})
        note = ("historical audit reference only; no pre-merge at any threshold" if t == 250
                else "threshold actually applied within participant groups")
        w.writerow([t, len(rows), sum(1 for r in rows if r["n_members"] > 1),
                    max(r["n_members"] for r in rows), occ, uniq, note])

# ---------- near-neighbour hints (resolvable <=20 SNP pairs across s20 units) ----------
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
                  "note": "near neighbour: no same-strain/same-participant source evidence -> NOT folded (task sec 2.6/A2)"})
with open(OUT / "near_neighbor_hints_v4.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, delimiter="\t", fieldnames=list(hints[0].keys()) if hints else
                       ["acc_a", "acc_b", "pairwise_snp", "unit_a", "unit_b",
                        "participant_a", "participant_b", "note"])
    w.writeheader()
    w.writerows(hints)

# ---------- unit diff vs v3 ----------
v3u = {}
for u in csv.DictReader(open(F2 / "01_canonical_units/canonical_units_s20.tsv"), delimiter="\t"):
    v3u[u["unit_id"]] = u
v4_of = a2u20
diff = []
for uid, u in sorted(v3u.items()):
    mem = u["members"].split(";")
    tgt = sorted({v4_of.get(m, "NOT_IN_V4") for m in mem})
    if len(tgt) == 1 and tgt[0] != "NOT_IN_V4":
        t0 = tgt[0]
        v4mem = next(r for r in unit_tables[20] if r["unit_id"] == t0)["members"].split(";")
        cat = "unchanged" if sorted(v4mem) == sorted(mem) else "grew_or_shrunk"
    else:
        cat = "split_or_mixed"
    diff.append({"v3_unit": uid, "v3_members": u["members"], "v3_representative": u["representative"],
                 "v4_units": ";".join(tgt), "change_category": cat})
# v4 units with no v3 counterpart member
v3mem_all = {m for u in v3u.values() for m in u["members"].split(";")}
for r in unit_tables[20]:
    mem = r["members"].split(";")
    if not any(m in v3mem_all for m in mem):
        diff.append({"v3_unit": "NEW_IN_V4", "v3_members": "", "v3_representative": "",
                     "v4_units": r["unit_id"], "change_category": "new_unit"})
with open(OUT / "unit_diff_v4_vs_v3.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, delimiter="\t", fieldnames=list(diff[0].keys()))
    w.writeheader()
    w.writerows(diff)

# ---------- summary ----------
for t in THRESHOLDS:
    rows = unit_tables[t]
    multi = [r for r in rows if r["n_members"] > 1]
    print(f"s{t}: units={len(rows)} multi={len(multi)} "
          f"identity={sum(1 for r in multi if r['unit_basis']=='identity')} "
          f"participant_clone={sum(1 for r in multi if r['unit_basis']=='participant_clone')}")
print("near-neighbour hints (<=20 SNP across units):", len(hints))
print("unit diff categories:", {c: sum(1 for d in diff if d['change_category'] == c)
                                for c in {d['change_category'] for d in diff}})
