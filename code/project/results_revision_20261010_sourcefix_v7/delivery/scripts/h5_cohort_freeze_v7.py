#!/usr/bin/env python3
"""Freeze version 7 cohorts from corrected source, QC and canonical-unit tables.

Representatives follow collection-date and quality rules. Separate cohort
outputs retain source eligibility and representative-selection provenance.
"""
import os
import csv
import hashlib
import json
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(os.environ["NEISSERIA_PROJECT_ROOT"])
OLD = ROOT / "results_revision_20261001"
FIX = ROOT / "results_revision_20261006_fix"
F2 = ROOT / "results_revision_20261008_fix2"
F3 = ROOT / "results_revision_20261009_fix3"
SUP = ROOT / "results_revision_20261010_sourcefix_v7"
V6 = ROOT / "results_revision_20261010_nsfix_v6"
OUT = SUP / "04_cohorts_v7"
OUT.mkdir(parents=True, exist_ok=True)

QC_COMP_MIN, QC_CONT_MAX, QC_N50_MIN = 90.0, 5.0, 10000
MARGIN_CONFIDENT = 0.50

# ---------- inputs ----------
rescreen = {r["accession"]: r for r in csv.DictReader(open(F3 / "01_reference_v4/rescreen_v4.tsv"), delimiter="\t")}
mag = {r["accession"]: r for r in csv.DictReader(open(SUP / "02_source_evidence/mag_source_audit_v7.tsv"), delimiter="\t")}
prov = {r["assembly_accession"]: r for r in csv.DictReader(open(SUP / "inputs_v7/provenance_corrected_v7.tsv"), delimiter="\t")}
qc = {r["accession"]: (float(r["Completeness"] or 0), float(r["Contamination"] or 0))
      for r in csv.DictReader(open(OLD / "01_retrieval_taxonomy/checkm2_merged.tsv"), delimiter="\t")}
hist_meta = {r["assembly_accession"]: r for r in csv.DictReader(
    open(ROOT / "neisseria_subflava_associated_group_metadata.tsv", encoding="utf-8-sig"), delimiter="\t")}
n50, nctg, gsize = {}, {}, {}
for r in csv.DictReader(open(F2 / "02_assembly_covariates/assembly_covariates_fasta.tsv"), delimiter="\t"):
    if r["n_contigs"] != "NA":
        n50[r["accession"]] = int(r["n50_bp"])
        nctg[r["accession"]] = int(r["n_contigs"])
        gsize[r["accession"]] = int(r["genome_size_bp"])

units = {u["unit_id"]: u for u in csv.DictReader(open(SUP / "03_units_v7/canonical_units_s20_v7.tsv"), delimiter="\t")}
unit_of = {}
for uid, u in units.items():
    for m in u["members"].split(";"):
        assert m not in unit_of, f"duplicate membership {m}"
        unit_of[m] = uid
rep_rank = {}
ledger_rows = {r["accession"]: r for r in csv.DictReader(open(SUP / "03_units_v7/representative_ledger_v7.tsv"), delimiter="\t")}
for a, r in ledger_rows.items():
    rp = r["rank_position"]
    rep_rank[a] = (r["unit_id"], int(rp) if rp not in ("NA", "") else 10**9)

def qc_of(acc):
    if acc in qc:
        return qc[acc][0], qc[acc][1], "checkm2_1.1.0_20261001"
    h = hist_meta.get(acc, {})
    if h.get("checkm_completeness"):
        return float(h["checkm_completeness"]), float(h.get("checkm_contamination", 0) or 0), "metadata_checkm"
    return None, None, "missing"

def qc_pass(acc):
    c, x, _ = qc_of(acc)
    if c is None:
        return False
    return c >= QC_COMP_MIN and x <= QC_CONT_MAX and n50.get(acc, 0) >= QC_N50_MIN

ledger = []
def exclude(acc, cohort, reason, detail=""):
    ledger.append({"accession": acc, "cohort": cohort, "reason": reason, "detail": detail})

def record_level_ok(acc, cid):
    d = mag.get(acc, {}).get("decision", "")
    if d == "EXCLUDE_confirmed_MAG":
        exclude(acc, cid, "MAG_confirmed", mag[acc]["decision_reason"]); return False
    if d == "REVIEW_undecided_excluded_from_main":
        exclude(acc, cid, "source_undecided", mag[acc]["decision_reason"]); return False
    df = prov.get(acc, {}).get("derivative_flag", "")
    if df == "experimental":
        exclude(acc, cid, "experimental_derivative", prov[acc]["strain"]); return False
    if df == "review_role_conflict":
        exclude(acc, cid, "source_role_review",
                f"source/role conflict unresolved: {prov[acc].get('strain', '')} "
                f"({prov[acc].get('supplement_note', 'natural-baseline role unconfirmed')})")
        return False
    return True

fallback_notes = []
def unit_representative(accs, cid):
    uid = unit_of.get(accs[0])
    ranked = sorted(accs, key=lambda a: (rep_rank.get(a, ("", 10**9))[1], a))
    keep = ranked[0]
    frozen_rep = units[uid]["representative"] if uid else ""
    if frozen_rep and frozen_rep != keep:
        reason = next((l["reason"] for l in ledger if l["accession"] == frozen_rep and l["cohort"] == cid),
                      "not eligible in this cohort")
        fallback_notes.append({"unit": uid, "cohort": cid, "frozen_representative": frozen_rep,
                               "reason_not_used": reason, "serving_representative": keep})
    return keep

def freeze_from_pool(pool, cid):
    members = []
    for acc in sorted(pool):
        if not record_level_ok(acc, cid):
            continue
        if not qc_pass(acc):
            exclude(acc, cid, "QC_fail"); continue
        members.append(acc)
    by_unit = {}
    for acc in members:
        by_unit.setdefault(unit_of.get(acc, f"singleton_{acc}"), []).append(acc)
    out = []
    for uid, accs in sorted(by_unit.items()):
        accs = sorted(accs)
        keep = accs[0] if uid.startswith("singleton_") else unit_representative(accs, cid)
        out.append(keep)
        for a in accs:
            if a != keep:
                exclude(a, cid, "non_representative_duplicate", f"unit={uid}")
    return sorted(set(out))

# ---------- pools ----------
pool_v4 = {a for a, r in rescreen.items() if r["associated_group_member_v4"] == "True"}
pool_exp = {a for a, r in rescreen.items() if r["associated_group_member_expanded"] == "True"}
for a in sorted(pool_exp - pool_v4):
    exclude(a, "revision_primary_v7", "ANI_membership_loss_type_anchor",
            "passes only via historical non-type perflava reference; moved to expanded_reference_sensitivity_v7")

revision_primary = freeze_from_pool(pool_v4, "revision_primary_v7")
expanded = freeze_from_pool(pool_exp, "expanded_reference_sensitivity_v7")

# ---------- legacy cohorts ----------
legacy139 = [r["assembly_accession"] for r in csv.DictReader(
    open(ROOT / "results_main139/00_cohort/neisseria_subflava_manuscript_main_139_metadata.csv"))]
boundary6 = [l.strip() for l in open(ROOT / "metadata/cohorts/manuscript_main_excluded_boundary_samples.txt") if l.strip()]
bs2acc = {v["biosample_accession"]: k for k, v in hist_meta.items()}
legacy145 = legacy139 + [bs2acc[b] for b in boundary6 if b in bs2acc]

def dereplicate(accs, cid):
    by_unit = {}
    for acc in sorted(set(accs)):
        if not record_level_ok(acc, cid):
            continue
        by_unit.setdefault(unit_of.get(acc, f"singleton_{acc}"), []).append(acc)
    out = []
    for uid, mem in sorted(by_unit.items()):
        mem = sorted(mem)
        keep = mem[0] if uid.startswith("singleton_") else unit_representative(mem, cid)
        out.append(keep)
        for a in mem:
            if a != keep:
                exclude(a, cid, "non_representative_duplicate", f"unit={uid}")
    return sorted(set(out))

legacy139_v7 = dereplicate(legacy139, "legacy139_dereplicated_v7")
legacy145_v7 = dereplicate(legacy145, "legacy_plus6_dereplicated_v7")

# ---------- assignment_confident_v7 ----------
assignment_confident = []
for acc in revision_primary:
    r = rescreen[acc]
    try:
        mg = float(r["margin_v4"])
    except (ValueError, TypeError):
        mg = None
    if mg is not None and mg >= MARGIN_CONFIDENT and r["group_for_testing_v4"] not in ("ambiguous", "", "ANI_UNRESOLVED"):
        assignment_confident.append(acc)

# ---------- one_per_participant_v7 ----------
by_pp = {}
for acc in revision_primary:
    p = prov.get(acc, {})
    by_pp.setdefault((p.get("study_namespace", "unknown"), p.get("participant_id", "unknown"),
                      p.get("participant_status", "unknown")), []).append(acc)
one_pp = []
import re as _re
for (ns, pid, status), accs in sorted(by_pp.items()):
    if status != "identifiable":
        one_pp.extend(accs)
        continue
    def date_key(a):
        d = hist_meta.get(a, {}).get("collection_date", "") or prov.get(a, {}).get("collection_date", "")
        m = _re.search(r"(\d{4})(?:[-/](\d{1,2}))?(?:[-/](\d{1,2}))?", d or "")
        if not m or (d or "").strip().lower() in ("unknown", "na", "none", "missing",
                                                  "not applicable", "not provided", "not collected"):
            return ((9999, 0, 0),)
        return ((int(m.group(1)), int(m.group(2) or 0), int(m.group(3) or 0)),)
    ranked = sorted(accs, key=lambda a: (date_key(a), rep_rank.get(a, ("", 10**9))[1], a))
    one_pp.append(ranked[0])
    for a in ranked[1:]:
        exclude(a, "one_per_participant_v7", "participant_collapsed", f"{ns}:{pid} kept {ranked[0]}")
one_pp = sorted(set(one_pp))

# ---------- sg0010 control (restated unchanged) ----------
sg0010 = [r["assembly_accession"] for r in csv.DictReader(open(OLD / "02_independence/cohort_sg0010_control.tsv"), delimiter="\t")]

# ---------- annotation union ----------
annotation_union = sorted(set(revision_primary) | set(legacy145_v7) | set(sg0010))

# ---------- write cohorts + metadata ----------
ledger_ncbi = {r["accession"]: r for r in csv.DictReader(open(OLD / "01_retrieval_taxonomy/full_ledger.tsv"), delimiter="\t")}
COLS = ["strain", "biosample", "group_for_testing", "unit_id", "unit_basis", "identity_group",
        "study_namespace", "submitted_bioproject", "true_source_bioproject", "source_role",
        "source_evidence_level", "participant_id", "participant_status", "participant_model_id",
        "host", "isolation_source", "host_disease", "geo_loc_name", "collection_date",
        "n_contigs", "genome_size", "n50_bp", "assembly_level"]

def meta_row(acc):
    h = hist_meta.get(acc, {})
    p = prov.get(acc, {})
    rs = rescreen.get(acc, {})
    lg = ledger_ncbi.get(acc, {})
    uid = unit_of.get(acc, "")
    ub = units[uid]["unit_basis"] if uid else ""
    return {"strain": acc, "biosample": p.get("biosample_accession") or h.get("biosample_accession") or lg.get("biosample", ""),
            "group_for_testing": rs.get("group_for_testing_v4", ""),
            "unit_id": uid, "unit_basis": ub, "identity_group": p.get("identity_group", ""),
            "study_namespace": p.get("study_namespace", ""),
            "submitted_bioproject": p.get("bioproject_accession", ""),
            "true_source_bioproject": p.get("true_source_bioproject", ""),
            "source_role": p.get("source_role", ""),
            "source_evidence_level": p.get("source_evidence_level", ""),
            "participant_id": p.get("participant_id", "unknown"),
            "participant_status": p.get("participant_status", "unknown"),
            "participant_model_id": p.get("participant_model_id", f"singleton:{acc}"),
            "host": h.get("host", "") or p.get("host", ""),
            "isolation_source": h.get("isolation_source", "") or p.get("isolation_source", ""),
            "host_disease": h.get("host_disease", "") or p.get("host_disease", ""),
            "geo_loc_name": h.get("geo_loc_name", "") or p.get("geo_loc_name", ""),
            "collection_date": h.get("collection_date", "") or p.get("collection_date", ""),
            "n_contigs": nctg.get(acc, "NA"), "genome_size": gsize.get(acc, "NA"),
            "n50_bp": n50.get(acc, "NA"),
            "assembly_level": h.get("assembly_level") or lg.get("assembly_level", "")}

cohorts = {"revision_primary_v7": revision_primary,
           "legacy139_dereplicated_v7": legacy139_v7,
           "legacy_plus6_dereplicated_v7": legacy145_v7,
           "assignment_confident_v7": assignment_confident,
           "one_per_participant_v7": one_pp,
           "expanded_reference_sensitivity_v7": expanded,
           "sg0010_control": sg0010,
           "annotation_union_v7": annotation_union}
for cid, accs in cohorts.items():
    with open(OUT / f"cohort_{cid}.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["accession"])
        for a in sorted(accs):
            w.writerow([a])
    with open(OUT / f"metadata_{cid}.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=COLS)
        w.writeheader()
        for a in sorted(accs):
            w.writerow(meta_row(a))
    print(cid, len(accs))

with open(OUT / "exclusion_ledger_v7.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, delimiter="\t", fieldnames=["accession", "cohort", "reason", "detail"])
    w.writeheader()
    w.writerows(ledger)
with open(OUT / "unit_fallback_notes_v7.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, delimiter="\t",
                       fieldnames=["unit", "cohort", "frozen_representative", "reason_not_used", "serving_representative"])
    w.writeheader()
    w.writerows(fallback_notes)

# ---------- cohort representatives per unit (S4 deliverable) ----------
crep = []
for cid in ("revision_primary_v7", "one_per_participant_v7", "expanded_reference_sensitivity_v7",
            "annotation_union_v7"):
    in_c = set(cohorts[cid])
    by_u = {}
    for a in in_c:
        by_u.setdefault(unit_of.get(a, f"singleton_{a}"), []).append(a)
    for uid, accs in sorted(by_u.items()):
        urow = units.get(uid)
        frozen = urow["representative"] if urow else ""
        fr = ledger_rows.get(frozen, {})
        if cid == "annotation_union_v7":
            # union is an annotation-prep superset: >1 assembly per unit allowed by design
            crep.append({
                "unit_id": uid, "cohort": cid,
                "frozen_representative": frozen,
                "frozen_source_eligible": fr.get("source_eligible", ""),
                "frozen_qc_pass": fr.get("qc_pass", ""),
                "frozen_membership_pass": fr.get("membership_pass", ""),
                "frozen_eligible_for_cohort": fr.get("eligible_for_cohort", ""),
                "serving_representative": ";".join(sorted(accs)),
                "fallback_reason": "annotation union: multiple assemblies per unit allowed; not a natural-isolate denominator",
            })
            continue
        serving = sorted(accs)[0]
        crep.append({
            "unit_id": uid, "cohort": cid,
            "frozen_representative": frozen,
            "frozen_source_eligible": fr.get("source_eligible", ""),
            "frozen_qc_pass": fr.get("qc_pass", ""),
            "frozen_membership_pass": fr.get("membership_pass", ""),
            "frozen_eligible_for_cohort": fr.get("eligible_for_cohort", ""),
            "serving_representative": serving,
            "fallback_reason": next((n["reason_not_used"] for n in fallback_notes
                                     if n["unit"] == uid and n["cohort"] == cid), ""),
        })
with open(OUT / "representatives_by_cohort_v7.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, delimiter="\t", fieldnames=list(crep[0].keys()))
    w.writeheader()
    w.writerows(crep)

# ---------- roster diff vs fix3 v4 ----------
v4map = {"revision_primary_v7": "cohort_revision_primary_v4.tsv",
         "legacy139_dereplicated_v7": "cohort_legacy139_dereplicated_v4.tsv",
         "legacy_plus6_dereplicated_v7": "cohort_legacy_plus6_dereplicated_v4.tsv",
         "assignment_confident_v7": "cohort_assignment_confident_v4.tsv",
         "one_per_participant_v7": "cohort_one_per_participant_v4.tsv",
         "expanded_reference_sensitivity_v7": "cohort_expanded_reference_sensitivity_v4.tsv",
         "annotation_union_v7": "cohort_annotation_union_v4.tsv"}

def driver(acc, cid, change):
    df = prov.get(acc, {}).get("derivative_flag", "")
    if df == "review_role_conflict":
        return "source_role_review_ATCC29256" if acc == "GCA_019334765.1" else "source_role_review"
    if df == "experimental":
        return "derivative_exclusion"
    ig = prov.get(acc, {}).get("identity_group", "")
    if ig.startswith("IGR1"):
        return "rochester_source_identity_reclassification"
    if ig:
        return "source_identity_group"
    p = prov.get(acc, {})
    if cid == "one_per_participant_v7":
        return "participant_rule_correction"
    return "unit_redefinition" if change else "rules_cascade"

with open(OUT / "roster_diff_v7_vs_fix3.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["cohort_v7", "accession", "change_vs_fix3", "primary_driver"])
    for cid, f in v4map.items():
        old = {r["accession"] for r in csv.DictReader(open(F3 / "04_cohorts_v4" / f), delimiter="\t")}
        new = set(cohorts[cid])
        for a in sorted(old - new):
            w.writerow([cid, a, "removed_vs_fix3", driver(a, cid, True)])
        for a in sorted(new - old):
            w.writerow([cid, a, "added_vs_fix3", driver(a, cid, True)])
        print(cid, "removed:", len(old - new), "added:", len(new - old))

with open(OUT / "roster_diff_expanded_vs_primary_v7.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["accession", "in_expanded", "in_primary", "note"])
    for a in sorted(set(expanded) | set(revision_primary)):
        w.writerow([a, a in expanded, a in revision_primary,
                    "historical-reference-only member" if a in expanded and a not in revision_primary else ""])

# ---------- FASTA hashes ----------
FASTA = {}
for d in (OLD / "01_retrieval_taxonomy/fasta", ROOT / "Genomes"):
    for p in d.glob("*.fna"):
        FASTA[p.stem] = p
STEM2ACC = {}
for r in csv.DictReader(open(OLD / "02_independence/retrieval_manifest.tsv"), delimiter="\t"):
    STEM2ACC[r["assembly_accession"]] = r["assembly_accession"]
    STEM2ACC[r["biosample"]] = r["assembly_accession"]
for r in hist_meta.values():
    acc = (r.get("assembly_accession") or "").strip()
    bs = (r.get("biosample_accession") or "").strip()
    fn = (r.get("fasta_filename") or "").strip()
    if acc:
        STEM2ACC[acc] = acc
        if bs:
            STEM2ACC[bs] = acc
        if fn and (ROOT / "Genomes" / fn).exists():
            FASTA.setdefault(acc, ROOT / "Genomes" / fn)

def fasta_of(acc):
    if acc in FASTA:
        return FASTA[acc]
    for stem, a in STEM2ACC.items():
        if a == acc and stem in FASTA:
            return FASTA[stem]
    return None

hash_targets = sorted(set(annotation_union) | set(expanded))
hash_rows = []
missing_fasta = []
for acc in hash_targets:
    p = fasta_of(acc)
    if p is None:
        missing_fasta.append(acc)
        hash_rows.append({"accession": acc, "fasta_path": "MISSING", "sha256": "MISSING"})
        continue
    hash_rows.append({"accession": acc, "fasta_path": str(p),
                      "sha256": hashlib.sha256(p.read_bytes()).hexdigest()})
with open(OUT / "fasta_hashes_v7.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, delimiter="\t", fieldnames=["accession", "fasta_path", "sha256"])
    w.writeheader()
    w.writerows(hash_rows)
print("fasta missing:", missing_fasta)

# ---------- freeze manifest ----------
def sha_file(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

roster_hashes = {cid: hashlib.sha256(("\n".join(sorted(accs))).encode()).hexdigest()
                 for cid, accs in cohorts.items()}
meta_hashes = {cid: sha_file(OUT / f"metadata_{cid}.tsv") for cid in cohorts}
input_files = [
    F3 / "01_reference_v4/rescreen_v4.tsv",
    SUP / "inputs_v7/provenance_corrected_v7.tsv",
    SUP / "inputs_v7/same_strain_identity_groups_v7.tsv",
    SUP / "03_units_v7/canonical_units_s20_v7.tsv",
    SUP / "03_units_v7/representative_ledger_v7.tsv",
    SUP / "02_source_evidence/source_role_corrections_v7.tsv",
    SUP / "02_source_evidence/umbrella_leaf_project_map_v7.tsv",
    SUP / "02_source_evidence/mag_source_audit_v7.tsv",
    SUP / "02_source_evidence/mag_exclusion_corrections_v7.tsv",
    V6 / "01_rochester_evidence/participant_source_evidence.tsv",
    V6 / "02_atcc29256/culture_identity_resolution.tsv",
    OLD / "01_retrieval_taxonomy/checkm2_merged.tsv",
    F2 / "02_assembly_covariates/assembly_covariates_fasta.tsv",
    F2 / "01_canonical_units/pairs_snp_completed.tsv",
    ROOT / "neisseria_subflava_associated_group_metadata.tsv",
]
manifest = {"freeze_id": "", "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "phase": "v7_phaseA_sourcefix_cohort_freeze",
            "supersedes": "v6freeze_20261010_09c8b78719ce (namespace fix); v7 supplement: "
                          "2 experimental transformants + 2 confirmed MAGs excluded, umbrella "
                          "leaf-project mapping 34 records, MAG audit rebound 2026-10-10",
            "cohort_sizes": {cid: len(accs) for cid, accs in cohorts.items()},
            "roster_hashes": roster_hashes, "metadata_hashes": meta_hashes,
            "input_hashes": {str(p): sha_file(p) for p in input_files},
            "reference_rule": "type-material-anchored bidirectional mean ANI>=95 & min AF>=0.70 "
                              "(subflava ATCC49275 / flavescens ATCC13120 / perflava CCUG17915 two "
                              "same-culture assemblies); historical non-type perflava 27098_8_142 = "
                              "control/sensitivity only",
            "qc_rule": "completeness>=90, contamination<=5, N50>=10000",
            "snp_clone_rule": "20 SNP complete linkage within identifiable-participant groups only; "
                              "source identity groups (34 Rochester duplicate-submission groups) "
                              "contracted separately; IGR135 MO05 same-isolate exception withdrawn "
                              "(participant code only) -> ordinary SNP rule; no 250-SNP pre-merge",
            "source_corrections": ["Rochester namespaces -> Regan2026_Rochester_study (Table S1 + SRA library "
                                   "evidence); MO -> ITM_PRJNA703317/PRJNA794044; GCA_019334765.1 ATCC29256 "
                                   "identity unresolved -> REVIEW, excluded from formal natural cohorts",
                                   "GCA_042536805.1/GCA_042536845.1 = experimental transformants "
                                   "(doi:10.3389/fmicb.2022.855482, PRJNA794044) -> experimental, "
                                   "excluded from all natural cohorts incl. legacy sensitivity",
                                   "25 umbrella-PRJNA514245 records mapped to true leaf projects "
                                   "(18 frozen BioSample + 7 live NCBI elink 2026-10-10)"],
            "fasta_hashes": "fasta_hashes_v7.tsv"}
freeze_id = "v7freeze_20261010_" + roster_hashes["revision_primary_v7"][:12]
manifest["freeze_id"] = freeze_id
with open(OUT / "freeze_manifest_v7.json", "w") as fh:
    json.dump(manifest, fh, indent=2)

# ---------- freeze report ----------
rep = []
ok = True
pri_excl = {l["accession"] for l in ledger if l["cohort"] == "revision_primary_v7"
            and l["reason"] != "ANI_membership_loss_type_anchor"}
conserved = (len(pool_v4) - len(pri_excl) == len(revision_primary)) and not (pool_v4 - pri_excl - set(revision_primary))
rep.append(f"primary arithmetic: pool_v4={len(pool_v4)} - exclusions={len(pri_excl)} = {len(revision_primary)} conserved={conserved} "
           f"(32 historical-only non-members logged separately as ANI_membership_loss_type_anchor)")
ok &= conserved
exp_excl = {l["accession"] for l in ledger if l["cohort"] == "expanded_reference_sensitivity_v7"}
cons2 = (len(pool_exp) - len(exp_excl) == len(expanded)) and not (pool_exp - exp_excl - set(expanded))
rep.append(f"expanded arithmetic: pool_exp={len(pool_exp)} - exclusions={len(exp_excl)} = {len(expanded)} conserved={cons2}")
ok &= cons2
rep.append(f"units among primary unique: {len({unit_of.get(a) for a in revision_primary}) == len([a for a in revision_primary if unit_of.get(a)])}")
rep.append(f"knockout GCA_052073905.1 in any natural cohort: "
           f"{any('GCA_052073905.1' in v for k, v in cohorts.items() if k != 'sg0010_control')}")
rep.append(f"ATCC29256 GCA_019334765.1 in any natural cohort: "
           f"{any('GCA_019334765.1' in v for k, v in cohorts.items() if k != 'sg0010_control')}")
rep.append(f"fasta missing for union+expanded: {missing_fasta}")
rep.append(f"freeze_id: {freeze_id}")
rep.append(f"cohort sizes: {json.dumps(manifest['cohort_sizes'])}")
rep.append(f"exclusions by reason: {dict(Counter(l['reason'] for l in ledger))}")
ok &= not missing_fasta
ok &= not any('GCA_052073905.1' in v for k, v in cohorts.items() if k != 'sg0010_control')
ok &= not any('GCA_019334765.1' in v for k, v in cohorts.items() if k != 'sg0010_control')
with open(OUT / "cohort_freeze_v7.md", "w") as fh:
    fh.write("# cohort freeze v5 (Phase A supplement)\n\n" + "\n".join(f"- {x}" for x in rep)
             + f"\n\nALL_CHECKS: {'PASS' if ok else 'FAIL'}\n")
print("\n".join(rep))
print("ALL_CHECKS:", "PASS" if ok else "FAIL")
