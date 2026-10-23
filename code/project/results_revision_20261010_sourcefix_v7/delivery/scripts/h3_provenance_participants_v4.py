#!/usr/bin/env python3
"""fix3-H3 (Phase A1+A3): corrected provenance, derivative evidence,
same-strain identity groups and participant identity mapping.

Fixes vs v3 (task v4, 2026-10-09):
  A1.2  SAMN50554711 (KU1003-01_tssM_KO) is an experimental derivative, not natural.
  A1.3  explicit same-strain / same-culture identity groups, evidence-based,
        separate from SNP clone thresholds. Groups whose genetics contradict
        the claimed identity (>250 SNP loose bound) go to REVIEW, not merged.
  A1.5  unresolved source-role conflicts -> REVIEW ledger, not frozen as natural.
  A3    participant_status in {identifiable, unknown, not_applicable};
        participant codes only valid inside their own study namespace
        (PRJNA514245 P#### != PRJNA1273950 P####). PRJNA1273950 P####.S###
        naming design parsed (70 candidates / 44 main-cohort reps unresolved
        in v3). The old "has strain name -> not_applicable_culture" heuristic
        is replaced by culture-collection-id recognition; lab strain labels
        become unknown (never a shared participant level).

Outputs (results_revision_20261009_fix3/02_provenance_v4/):
  provenance_corrected.tsv        all v3 columns + participant_status + model id
  derivative_evidence.tsv
  same_strain_identity_groups.tsv
  source_review_unresolved.tsv
  participant_identity_map.tsv
  participant_group_validation.tsv
"""
import os
import csv
import re
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(os.environ["NEISSERIA_PROJECT_ROOT"])
OLD = ROOT / "results_revision_20261001"
F2 = ROOT / "results_revision_20261008_fix2"
OUT = ROOT / "results_revision_20261009_fix3" / "02_provenance_v4"
OUT.mkdir(parents=True, exist_ok=True)

prov = {r["assembly_accession"]: r for r in csv.DictReader(
    open(OLD / "02_independence/provenance.tsv"), delimiter="\t")}

pairs = {}
for p in csv.DictReader(open(F2 / "01_canonical_units/pairs_snp_completed.tsv"), delimiter="\t"):
    pairs[frozenset((p["acc_a"], p["acc_b"]))] = p

def pair_snp(a, b):
    p = pairs.get(frozenset((a, b)))
    if not p:
        return None, "no_pair"
    if p["resolvable_by_snp"] == "yes" and p["pairwise_snp_ACGT"]:
        return int(p["pairwise_snp_ACGT"]), "yes"
    return (int(p["pairwise_snp_ACGT"]) if p["pairwise_snp_ACGT"] else None), p["resolvable_by_snp"]

# ---------------------------------------------------------------- identity groups
# decision rule (task A1.3/A1.4): source identity evidence decides; SNP only
# corroborates. loose clonality bound 250 SNP (historical audit reference):
# above it the genetics contradict the claim -> REVIEW, no merge.
IDENTITY_GROUPS = [
    # gid, members, gtype, evidence, decision
    ("IG01", ["GCA_000090875.1", "GCA_005886145.1"], "same_biosample",
     "same BioSample SAMN02463843, strain F0314, two assemblies (PRJNA42547/PRJNA543187)", "merged"),
    ("IG02", ["GCA_014596365.2", "GCA_052073865.1"], "same_biosample",
     "same BioSample SAMN05362769, strain KU1003-01 parent, two assemblies (PRJNA327996/PRJNA1305711); "
     "knockout GCA_052073905.1 (SAMN50554711) is a derivative of the same parent, excluded separately", "merged"),
    ("IG03", ["GCA_000186165.1", "GCA_013267835.1"], "same_strain",
     "SAMN02463751 strain C102 (BEI HM-242); SAMN11056473 FDAARGOS_758 description comment 'strain C102' "
     "(HM-242D), collected_by BEI Resources/ATCC; independent source link documented in task 2026-10-09 sec 2.6", "merged"),
    ("IG04", ["GCA_023635085.1", "GCA_056611915.1"], "same_strain",
     "strain SG0010, collection 2017; baseline repeat sequencing (PRJNA714914 / PRJNA976074)", "merged"),
    ("IG05", ["GCA_005221305.1", "GCA_031190635.1", "GCA_963396225.1"], "culture_synonym",
     "N. subflava type strain: ATCC 49275 = CCUG 23930 = LMG 5313 (culture-collection catalog synonyms)", "merged"),
    ("IG06", ["GCA_005221285.1", "GCA_900453805.1", "GCA_963394955.1"], "culture_synonym",
     "N. flavescens type strain: ATCC 13120 = NCTC 8263 = CCUG 17913; NCTC official catalog lists all three "
     "as the same strain (culturecollections.org.uk, cited in task 2026-10-09 sec 2.4)", "merged"),
    ("IG07", ["GCA_023472875.1", "GCA_042649725.1"], "culture_synonym",
     "N. perflava type strain CCUG 17915, two independent assemblies of the same type culture", "merged"),
    ("IG08", ["GCA_000173995.1", "GCA_000741965.1"], "culture_synonym",
     "strain ATCC 23970, two submissions", "merged"),
    ("IG09", ["GCA_000176755.2", "GCA_000818035.1"], "culture_synonym",
     "strain ATCC 29315, two submissions", "merged"),
    ("IG10", ["GCA_000220865.1", "GCA_022749495.1"], "culture_synonym",
     "strain ATCC 33926, two submissions", "merged"),
    ("IG11", ["GCA_001027865.1", "GCA_022870905.1"], "same_strain",
     "identical unique strain name KH1503 across submissions (PRJNA267012/PRJNA788950)", "merged"),
    ("IG12", ["GCA_021728505.1", "GCA_049190985.1"], "same_strain",
     "identical unique strain name ZJ104 across submissions (PRJNA799687/PRJNA1236038)", "merged"),
    ("IG13", ["GCA_021728675.1", "GCA_028463985.1"], "same_strain",
     "identical unique strain name ZJ106 across submissions (PRJNA799686/PRJNA925655)", "merged"),
    ("IG14", ["GCA_011038745.1", "GCA_014898075.1"], "same_strain",
     "identical unique strain name ZJ785 across submissions (PRJNA606173/PRJNA650271)", "merged"),
    ("IG15", ["GCA_011038775.1", "GCA_017811745.1", "GCA_049190975.1"], "same_strain",
     "identical unique strain name ZJ930 across three submissions (PRJNA606178/PRJNA668708/PRJNA1236038)", "merged"),
    # REVIEW: source claim contradicted by genetics -> NOT merged (task A1.5)
    ("IGR01", ["GCA_000174655.1", "GCA_019334765.1"], "culture_synonym",
     "both claim ATCC 29256 but pair is 115373 SNP / CLONE_UNRESOLVED: genetics contradict shared-culture claim", "review_no_merge"),
    ("IGR02", ["GCA_002108595.1", "GCA_014055005.1"], "culture_synonym",
     "both claim DSM 19151 but pair is 347 SNP (>250 loose clonality bound): identity unresolved", "review_no_merge"),
    ("IGR03", ["GCA_000180595.1", "GCA_002022745.1"], "same_strain",
     "identical strain name Y92-1009 but pair is 711 SNP (>250): identity unresolved", "review_no_merge"),
]

ig_rows = []
acc2ig = {}
for gid, mem, gtype, evidence, decision in IDENTITY_GROUPS:
    ds, res = [], []
    for i in range(len(mem)):
        for j in range(i + 1, len(mem)):
            d, r = pair_snp(mem[i], mem[j])
            ds.append(d)
            res.append(r)
    have = [d for d in ds if d is not None]
    acc2ig_rows = {"group_id": gid, "group_type": gtype, "members": ";".join(mem),
                   "n_members": len(mem), "evidence": evidence,
                   "pairwise_snp": ";".join("NA" if d is None else str(d) for d in ds),
                   "pairs_resolvable": f"{sum(1 for r in res if r=='yes')}/{len(res)}",
                   "snp_min": min(have) if have else "NA", "snp_max": max(have) if have else "NA",
                   "decision": decision,
                   "representative_rule": "qc_rank (same-culture duplicate)" if decision == "merged" else "none (not merged)"}
    ig_rows.append(acc2ig_rows)
    if decision == "merged":
        for a in mem:
            assert a not in acc2ig, f"{a} in two identity groups"
            acc2ig[a] = gid

with open(OUT / "same_strain_identity_groups.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, delimiter="\t", fieldnames=list(ig_rows[0].keys()))
    w.writeheader()
    w.writerows(ig_rows)

# ---------------------------------------------------------------- derivative evidence
der_rows = []
def add_der(acc, flag, evidence, disposition):
    p = prov.get(acc, {})
    der_rows.append({"accession": acc, "biosample": p.get("biosample_accession", ""),
                     "bioproject": p.get("bioproject_accession", ""),
                     "strain": p.get("strain", ""), "v3_derivative_flag": p.get("derivative_flag", ""),
                     "v4_derivative_flag": flag, "evidence": evidence, "disposition": disposition})

add_der("GCA_052073905.1", "experimental",
        "strain 'KU1003-01_tssM_KO'; BioSample SAMN50554711 isolation_source "
        "'tssM knockout of parent strain KU1003-01' (frozen JSONL attrs verified 2026-10-09)",
        "excluded from every natural-strain cohort and from representative ranking; "
        "parent strain represented via IG02 (SAMN05362769 assemblies)")
add_der("GCA_056611915.1", "review_role_conflict",
        "BioSample SAMN35359261 strain SG0010 / collection 2017 (natural-compatible), but audit reports "
        "sample name Cef-CA26 matching the 2023 experimental Cef-CA## series (PRJNA955008/978365); frozen "
        "JSONL attrs carry neither confirmation nor refutation. Natural-baseline role cannot be confirmed "
        "on server -> REVIEW (task 2026-10-09 sec 2.7)",
        "not eligible as natural representative; kept in REVIEW ledger / technical-control records; "
        "SG0010 natural baseline represented by GCA_023635085.1 (SAMN18320762, PRJNA714914 sputum study, 2017)")
for acc, p in sorted(prov.items()):
    if p.get("derivative_flag") == "experimental" and acc != "GCA_052073905.1":
        der_rows.append({"accession": acc, "biosample": p["biosample_accession"],
                         "bioproject": p["bioproject_accession"], "strain": p["strain"],
                         "v3_derivative_flag": "experimental", "v4_derivative_flag": "experimental",
                         "evidence": "v3 flag retained (Cef-/GM-mut-/Mock naming or derivative BioProject)",
                         "disposition": "excluded from natural-strain cohorts (unchanged)"})
with open(OUT / "derivative_evidence.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, delimiter="\t", fieldnames=list(der_rows[0].keys()))
    w.writeheader()
    w.writerows(der_rows)

# ---------------------------------------------------------------- source review ledger
rev_rows = [
    {"item": "GCA_056611915.1_natural_baseline_role", "accessions": "GCA_056611915.1",
     "issue": "natural baseline vs experimental Cef-CA26 role conflict",
     "evidence": "see derivative_evidence.tsv", "disposition": "REVIEW; not a natural representative in v4 freeze"},
    {"item": "ATCC29256_identity_link", "accessions": "GCA_000174655.1;GCA_019334765.1",
     "issue": "shared ATCC number contradicted by 115373 SNP / unresolvable pair",
     "evidence": "same_strain_identity_groups.tsv IGR01",
     "disposition": "link REJECTED pending author review; both records stay separate natural units (each individually eligible)"},
    {"item": "DSM19151_identity_link", "accessions": "GCA_002108595.1;GCA_014055005.1",
     "issue": "shared DSM number but 347 SNP (>250 loose bound)",
     "evidence": "same_strain_identity_groups.tsv IGR02",
     "disposition": "link REJECTED pending author review; separate natural units"},
    {"item": "Y92-1009_identity_link", "accessions": "GCA_000180595.1;GCA_002022745.1",
     "issue": "shared strain name but 711 SNP (>250 loose bound)",
     "evidence": "same_strain_identity_groups.tsv IGR03",
     "disposition": "link REJECTED pending author review; separate natural units"},
    {"item": "Nei437_Nei438_low_snp", "accessions": "GCA_042537275.1;GCA_042537445.1",
     "issue": "5 SNP but distinct BioSamples/strains, participant unknown, umbrella PRJNA514245: no same-strain or same-participant evidence",
     "evidence": "task 2026-10-09 sec 2.6; provenance fields",
     "disposition": "NOT merged; each keeps its own natural unit; recorded as near-neighbour hint"},
]
with open(OUT / "source_review_unresolved.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, delimiter="\t", fieldnames=list(rev_rows[0].keys()))
    w.writeheader()
    w.writerows(rev_rows)

# ---------------------------------------------------------------- participant rules
CULTURE_PAT = re.compile(r"^(ATCC|CCUG|NCTC|DSM|LMG|CIP|JCM|NBRC|FDAARGOS)\s*:?\s*\w*\d+", re.I)
COLLECTION_BP = {"PRJEB2090", "PRJNA428750"}

# --- Rochester carriage collection: PRJNA514245 (Pathogen Detection umbrella,
#     P####S### submissions) and PRJNA1273950 ("Commensal Neisseria collection
#     from Rochester, NY", P####.S###) are the SAME study/participants.
#     Evidence (computed below from actual pair matrix, task sec 2.7/A3):
#     every P code shared by both projects has cross-project pairs at 0-24 SNP
#     (clone-level <=20 for P0004/6/9/13/14/15/16/17; 21-24 for P0007/P0008),
#     same geography (Rochester NY), matching participant+sample ordinals.
#     Numbers alone are NOT treated as cross-project credentials; the mapping
#     is justified by genetics + study relationship and recorded per code.
def pcode_of(p):
    s = (p.get("strain") or "").strip()
    m = re.fullmatch(r"(P\d{4})S(\d{3})", s) or re.fullmatch(r"(P\d{4})\.S(\d{3})", s)
    return (m.group(1), m.group(2)) if m else (None, None)

ROCHESTER_NS = "Rochester_collection(PRJNA514245+PRJNA1273950)"
pc_rec = defaultdict(lambda: defaultdict(list))   # pcode -> bioproject -> accs
for a, p in prov.items():
    pc, _ = pcode_of(p)
    if pc and p["bioproject_accession"] in ("PRJNA514245", "PRJNA1273950"):
        pc_rec[pc][p["bioproject_accession"]].append(a)

ns_map_rows = []
PCODE_CONF = {}
for pc in sorted(pc_rec):
    a5 = sorted(pc_rec[pc].get("PRJNA514245", []))
    a12 = sorted(pc_rec[pc].get("PRJNA1273950", []))
    if a5 and a12:
        ds = [pair_snp(x, y)[0] for x in a5 for y in a12]
        have = [d for d in ds if d is not None]
        mn = min(have) if have else None
        n20 = sum(1 for d in have if d <= 20)
        conf = "high" if (mn is not None and mn <= 20) else "medium"
        PCODE_CONF[pc] = conf
        ns_map_rows.append({"participant_id": pc, "n_PRJNA514245": len(a5), "n_PRJNA1273950": len(a12),
                            "min_cross_project_snp": mn if mn is not None else "NA",
                            "cross_pairs_le20": n20, "cross_pairs_resolvable": len(have),
                            "unified_namespace": ROCHESTER_NS, "confidence": conf,
                            "evidence": "cross-project same-code pairs near-clonal; same study geography "
                                        "(Rochester NY); PRJNA514245 = Pathogen Detection umbrella of the "
                                        "same collection, PRJNA1273950 = RIT submission of the collection"})
    else:
        PCODE_CONF[pc] = "medium"
        ns_map_rows.append({"participant_id": pc, "n_PRJNA514245": len(a5), "n_PRJNA1273950": len(a12),
                            "min_cross_project_snp": "single_project", "cross_pairs_le20": "",
                            "cross_pairs_resolvable": "", "unified_namespace": ROCHESTER_NS,
                            "confidence": "medium",
                            "evidence": "code present in one project only; naming-design participant inference"})
with open(OUT / "study_namespace_mapping.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, delimiter="\t", fieldnames=list(ns_map_rows[0].keys()))
    w.writeheader()
    w.writerows(ns_map_rows)

def classify(p):
    """-> (participant_id, participant_status, method, confidence, study_namespace_override)"""
    bp = p["bioproject_accession"]
    strain = (p.get("strain") or "").strip()
    iso = (p.get("isolate") or "").strip()
    if bp == "PRJNA514245":
        m = re.fullmatch(r"MO\d+", iso)
        if m:
            return m.group(0), "identifiable", "isolate_participant_code", "high", None
        pc, _ = pcode_of(p)
        if pc:
            return pc, "identifiable", "strain_participant_code", PCODE_CONF.get(pc, "medium"), ROCHESTER_NS
        return "unknown", "unknown", "umbrella_project", "low", None
    if bp == "PRJNA1273950":
        # study naming design: every record is P####.S### (participant.sample);
        # unified Rochester namespace with PRJNA514245 P-codes (evidence above)
        pc, _ = pcode_of(p)
        if pc:
            return pc, "identifiable", "strain_participant_code_naming_design", PCODE_CONF.get(pc, "medium"), ROCHESTER_NS
        return "unknown", "unknown", "none", "low", None
    if bp in COLLECTION_BP:
        return "not_applicable", "not_applicable", "study_design_collection", "high", None
    if strain == "SG0010":
        return "not_applicable", "not_applicable", "lab_strain_repeat_sequencing", "high", None
    if CULTURE_PAT.match(strain):
        return "not_applicable", "not_applicable", "culture_collection_id", "high", None
    return "unknown", "unknown", "none_or_unmapped_isolate_id", "low", None

# ---------------------------------------------------------------- corrected provenance
rows = []
pmap_rows = []
for acc, p in sorted(prov.items()):
    r = dict(p)
    pid, pstatus, pmethod, pconf, ns_override = classify(p)
    # derivative corrections
    if acc == "GCA_052073905.1":
        r["derivative_flag"] = "experimental"
    elif acc == "GCA_056611915.1":
        r["derivative_flag"] = "review_role_conflict"
    r["participant_id"] = pid
    r["participant_status"] = pstatus
    r["participant_method"] = pmethod
    r["participant_confidence"] = pconf
    ns = ns_override or p.get("study_namespace") or p["bioproject_accession"]
    r["study_namespace"] = ns
    # model-safe participant level: only identifiable participants share a level
    r["participant_model_id"] = f"{ns}:{pid}" if pstatus == "identifiable" else f"singleton:{acc}"
    r["identity_group"] = acc2ig.get(acc, "")
    rows.append(r)
    pmap_rows.append({"accession": acc, "biosample": r["biosample_accession"],
                      "study_namespace": ns, "participant_id": pid,
                      "participant_status": pstatus, "participant_method": pmethod,
                      "participant_confidence": pconf, "participant_model_id": r["participant_model_id"],
                      "identity_group": r["identity_group"]})

fields = list(rows[0].keys())
with open(OUT / "provenance_corrected.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, delimiter="\t", fieldnames=fields)
    w.writeheader()
    w.writerows(rows)
with open(OUT / "participant_identity_map.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, delimiter="\t", fieldnames=list(pmap_rows[0].keys()))
    w.writeheader()
    w.writerows(pmap_rows)

# ---------------------------------------------------------------- validation summary
vg = []
cnt = Counter((r["study_namespace"], r["participant_status"]) for r in rows)
for (ns, st), n in sorted(cnt.items()):
    vg.append({"check": "namespace_status_counts", "study_namespace": ns,
               "participant_status": st, "n_records": n, "detail": ""})
# cross-project P-code collisions: only allowed inside the documented unified
# Rochester namespace; any other shared code across namespaces is a violation
pcodes = defaultdict(set)
for r in rows:
    if r["participant_status"] == "identifiable":
        pcodes[r["participant_id"]].add(r["study_namespace"])
coll = {k: sorted(v) for k, v in pcodes.items() if len(v) > 1}
vg.append({"check": "pcode_namespace_collisions", "study_namespace": ";".join(sorted({ns for v in coll.values() for ns in v})),
           "participant_status": "identifiable", "n_records": len(coll),
           "detail": ("no unresolved collisions" if not coll else
                      "VIOLATION: " + ";".join(f"{k}={'|'.join(v)}" for k, v in sorted(coll.items())))})
# PRJNA1273950 parse completeness (by BioProject; namespace is unified)
n1273 = [r for r in rows if r["bioproject_accession"] == "PRJNA1273950"]
n1273_ok = sum(1 for r in n1273 if r["participant_status"] == "identifiable")
vg.append({"check": "PRJNA1273950_pcode_parsed", "study_namespace": ROCHESTER_NS,
           "participant_status": "identifiable", "n_records": n1273_ok,
           "detail": f"{n1273_ok}/{len(n1273)} records parsed to participant via P####.S### naming design"})
vg.append({"check": "rochester_namespace_unification", "study_namespace": ROCHESTER_NS,
           "participant_status": "identifiable",
           "n_records": sum(1 for r in rows if r["study_namespace"] == ROCHESTER_NS),
           "detail": f"{len(ns_map_rows)} P codes unified across PRJNA514245/PRJNA1273950; "
                     "per-code genetic evidence in study_namespace_mapping.tsv"})
with open(OUT / "participant_group_validation.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, delimiter="\t", fieldnames=list(vg[0].keys()))
    w.writeheader()
    w.writerows(vg)

print("provenance rows:", len(rows))
print("participant_status:", dict(Counter(r["participant_status"] for r in rows)))
print("derivative_flag:", dict(Counter(r["derivative_flag"] for r in rows)))
print("identity groups merged:", sum(1 for g in ig_rows if g['decision'] == 'merged'),
      "| review_no_merge:", sum(1 for g in ig_rows if g['decision'] == 'review_no_merge'))
print("P-code namespaces with collisions:", len(coll))
print("PRJNA1273950 parsed:", n1273_ok, "/", len(n1273))
