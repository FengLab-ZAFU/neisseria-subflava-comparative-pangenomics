#!/usr/bin/env python3
"""fix3-H5 (Phase A4): cohort freeze v4 on corrected units / participants /
type-anchored reference rule.

Record-level rules (unchanged project conventions): MAG audit decisions,
source REVIEW records stay out of the main freeze, experimental derivatives
and role-conflict records excluded from natural cohorts; QC: completeness
>=90%, contamination <=5%, N50 >=10,000 bp (CheckM2, historical-metadata
fallback). Membership: v4 type-material-anchored rule (h2); the historical
non-type perflava reference only defines expanded_reference_sensitivity_v4.

Cohorts (sizes NOT preset; they fall out of rules + data):
  revision_primary_v4, legacy139_dereplicated_v4, legacy_plus6_dereplicated_v4,
  assignment_confident_v4, one_per_participant_v4,
  expanded_reference_sensitivity_v4, annotation_union_v4,
  sg0010_control (restated unchanged - same inputs)

Outputs (results_revision_20261009_fix3/04_cohorts_v4/):
  cohort_<id>.tsv, metadata_<id>.tsv, exclusion_ledger_v4.tsv,
  roster_diff_v4_vs_v3.tsv, unit_fallback_notes_v4.tsv,
  fasta_hashes_v4.tsv, freeze_manifest_v4.json, cohort_freeze_v4.md
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
OUT = F3 / "04_cohorts_v4"
OUT.mkdir(parents=True, exist_ok=True)

QC_COMP_MIN, QC_CONT_MAX, QC_N50_MIN = 90.0, 5.0, 10000
MARGIN_CONFIDENT = 0.50

# ---------- inputs ----------
rescreen = {r["accession"]: r for r in csv.DictReader(open(F3 / "01_reference_v4/rescreen_v4.tsv"), delimiter="\t")}
mag = {r["accession"]: r for r in csv.DictReader(open(FIX / "01_mag_audit/mag_source_audit.tsv"), delimiter="\t")}
prov = {r["assembly_accession"]: r for r in csv.DictReader(open(F3 / "02_provenance_v4/provenance_corrected.tsv"), delimiter="\t")}
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

units = {u["unit_id"]: u for u in csv.DictReader(open(F3 / "03_units_v4/canonical_units_s20_v4.tsv"), delimiter="\t")}
unit_of = {}
for uid, u in units.items():
    for m in u["members"].split(";"):
        assert m not in unit_of, f"duplicate membership {m}"
        unit_of[m] = uid
rep_rank = {}   # accession -> (unit_id, rank_position among natural-eligible)
for r in csv.DictReader(open(F3 / "03_units_v4/representative_ledger_v4.tsv"), delimiter="\t"):
    rp = r["rank_position"]
    rep_rank[r["accession"]] = (r["unit_id"], int(rp) if rp not in ("NA", "") else 10**9)

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
        exclude(acc, cid, "source_role_review", "natural-baseline role unconfirmed (SG0010/Cef-CA26)"); return False
    return True

fallback_notes = []
def unit_representative(accs, cid):
    """choose among record-level+QC-passing members of one unit by the frozen
    representative ledger rank (date/QC rule from h4)."""
    uid = unit_of.get(accs[0])
    ranked = sorted(accs, key=lambda a: (rep_rank.get(a, ("", 10**9))[1], a))
    keep = ranked[0]
    frozen_rep = units[uid]["representative"] if uid else ""
    if frozen_rep and frozen_rep != keep:
        reason = next((l["reason"] for l in ledger if l["accession"] == frozen_rep and l["cohort"] == cid),
                      "not eligible in this cohort")
        fallback_notes.append({"unit": uid, "frozen_representative": frozen_rep,
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
    exclude(a, "revision_primary_v4", "ANI_membership_loss_type_anchor",
            "passes only via historical non-type perflava reference; moved to expanded_reference_sensitivity_v4")

revision_primary = freeze_from_pool(pool_v4, "revision_primary_v4")
expanded = freeze_from_pool(pool_exp, "expanded_reference_sensitivity_v4")

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

legacy139_v4 = dereplicate(legacy139, "legacy139_dereplicated_v4")
legacy145_v4 = dereplicate(legacy145, "legacy_plus6_dereplicated_v4")

# ---------- assignment_confident_v4 ----------
assignment_confident = []
for acc in revision_primary:
    r = rescreen[acc]
    try:
        mg = float(r["margin_v4"])
    except (ValueError, TypeError):
        mg = None
    if mg is not None and mg >= MARGIN_CONFIDENT and r["group_for_testing_v4"] not in ("ambiguous", "", "ANI_UNRESOLVED"):
        assignment_confident.append(acc)

# ---------- one_per_participant_v4 ----------
by_pp = {}
for acc in revision_primary:
    p = prov.get(acc, {})
    by_pp.setdefault((p.get("study_namespace", "unknown"), p.get("participant_id", "unknown"),
                      p.get("participant_status", "unknown")), []).append(acc)
one_pp = []
for (ns, pid, status), accs in sorted(by_pp.items()):
    if status != "identifiable":
        one_pp.extend(accs)
        continue
    def date_key(a):
        d = hist_meta.get(a, {}).get("collection_date", "") or prov.get(a, {}).get("collection_date", "")
        import re as _re
        m = _re.search(r"(\d{4})(?:[-/](\d{1,2}))?(?:[-/](\d{1,2}))?", d or "")
        if not m or (d or "").strip().lower() in ("unknown", "na", "none", "missing",
                                                  "not applicable", "not provided", "not collected"):
            return ((9999, 0, 0),)
        return ((int(m.group(1)), int(m.group(2) or 0), int(m.group(3) or 0)),)
    ranked = sorted(accs, key=lambda a: (date_key(a), rep_rank.get(a, ("", 10**9))[1], a))
    one_pp.append(ranked[0])
    for a in ranked[1:]:
        exclude(a, "one_per_participant_v4", "participant_collapsed", f"{ns}:{pid} kept {ranked[0]}")
one_pp = sorted(set(one_pp))

# ---------- sg0010 control (restated unchanged) ----------
sg0010 = [r["assembly_accession"] for r in csv.DictReader(open(OLD / "02_independence/cohort_sg0010_control.tsv"), delimiter="\t")]

# ---------- annotation union ----------
annotation_union = sorted(set(revision_primary) | set(legacy145_v4) | set(sg0010))

# ---------- write cohorts + metadata ----------
ledger_ncbi = {r["accession"]: r for r in csv.DictReader(open(OLD / "01_retrieval_taxonomy/full_ledger.tsv"), delimiter="\t")}
COLS = ["strain", "biosample", "group_for_testing", "unit_id", "unit_basis", "identity_group",
        "study_namespace", "participant_id", "participant_status", "participant_model_id",
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

cohorts = {"revision_primary_v4": revision_primary,
           "legacy139_dereplicated_v4": legacy139_v4,
           "legacy_plus6_dereplicated_v4": legacy145_v4,
           "assignment_confident_v4": assignment_confident,
           "one_per_participant_v4": one_pp,
           "expanded_reference_sensitivity_v4": expanded,
           "sg0010_control": sg0010,
           "annotation_union_v4": annotation_union}
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

with open(OUT / "exclusion_ledger_v4.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, delimiter="\t", fieldnames=["accession", "cohort", "reason", "detail"])
    w.writeheader()
    w.writerows(ledger)
with open(OUT / "unit_fallback_notes_v4.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, delimiter="\t",
                       fieldnames=["unit", "frozen_representative", "reason_not_used", "serving_representative"])
    w.writeheader()
    w.writerows(fallback_notes)

# ---------- roster diff vs v3 ----------
v3map = {"revision_primary_v4": "cohort_revision_primary_v3.tsv",
         "legacy139_dereplicated_v4": "cohort_legacy139_dereplicated_v3.tsv",
         "legacy_plus6_dereplicated_v4": "cohort_legacy_plus6_dereplicated_v3.tsv",
         "assignment_confident_v4": "cohort_assignment_confident_v3.tsv",
         "one_per_participant_v4": "cohort_one_per_participant_v3.tsv",
         "annotation_union_v4": "cohort_annotation_union_v3.tsv"}
v3units = {}
for u in csv.DictReader(open(F2 / "01_canonical_units/canonical_units_s20.tsv"), delimiter="\t"):
    for m in u["members"].split(";"):
        v3units[m] = u["unit_id"]

def driver(acc, cid, change):
    df = prov.get(acc, {}).get("derivative_flag", "")
    if df == "experimental":
        return "derivative_exclusion"
    if df == "review_role_conflict":
        return "source_role_review"
    if acc in pool_exp and acc not in pool_v4:
        return "reference_anchor_rule"
    ig = prov.get(acc, {}).get("identity_group", "")
    if ig:
        return "source_identity_group"
    p = prov.get(acc, {})
    if cid == "one_per_participant_v4":
        return "participant_rule_correction"
    if v3units.get(acc) and unit_of.get(acc) and v3units[acc] != unit_of[acc]:
        return "unit_redefinition"
    return "rules_cascade"

with open(OUT / "roster_diff_v4_vs_v3.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["cohort_v4", "accession", "change", "primary_driver"])
    for cid, f in v3map.items():
        old = {r["accession"] for r in csv.DictReader(open(F2 / "03_cohorts" / f), delimiter="\t")}
        new = set(cohorts[cid])
        for a in sorted(old - new):
            w.writerow([cid, a, "removed_vs_v3", driver(a, cid, "removed")])
        for a in sorted(new - old):
            w.writerow([cid, a, "added_vs_v3", driver(a, cid, "added")])
        print(cid, "removed:", len(old - new), "added:", len(new - old))

# expanded sensitivity never existed in v3: diff vs v3 primary for context
with open(OUT / "roster_diff_expanded_vs_primary_v4.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["accession", "in_expanded", "in_primary", "note"])
    for a in sorted(set(expanded) | set(revision_primary)):
        w.writerow([a, a in expanded, a in revision_primary,
                    "historical-reference-only member" if a in expanded and a not in revision_primary else ""])

# ---------- FASTA hashes for the annotation union (superset of all cohorts) ----------
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
with open(OUT / "fasta_hashes_v4.tsv", "w", newline="") as fh:
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
    F3 / "02_provenance_v4/provenance_corrected.tsv",
    F3 / "02_provenance_v4/same_strain_identity_groups.tsv",
    F3 / "03_units_v4/canonical_units_s20_v4.tsv",
    F3 / "03_units_v4/representative_ledger_v4.tsv",
    FIX / "01_mag_audit/mag_source_audit.tsv",
    OLD / "01_retrieval_taxonomy/checkm2_merged.tsv",
    F2 / "02_assembly_covariates/assembly_covariates_fasta.tsv",
    F2 / "01_canonical_units/pairs_snp_completed.tsv",
    ROOT / "neisseria_subflava_associated_group_metadata.tsv",
]
manifest = {"freeze_id": "", "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "phase": "v4_phaseA_cohort_freeze",
            "cohort_sizes": {cid: len(accs) for cid, accs in cohorts.items()},
            "roster_hashes": roster_hashes, "metadata_hashes": meta_hashes,
            "input_hashes": {str(p): sha_file(p) for p in input_files},
            "reference_rule": "type-material-anchored bidirectional mean ANI>=95 & min AF>=0.70 "
                              "(subflava ATCC49275 / flavescens ATCC13120 / perflava CCUG17915 two "
                              "same-culture assemblies); historical non-type perflava 27098_8_142 = "
                              "control/sensitivity only",
            "qc_rule": "completeness>=90, contamination<=5, N50>=10000",
            "snp_clone_rule": "20 SNP complete linkage within identifiable-participant groups only; "
                              "source identity groups separate; no 250-SNP pre-merge",
            "fasta_hashes": "fasta_hashes_v4.tsv"}
freeze_id = "v4freeze_20261009_" + roster_hashes["revision_primary_v4"][:12]
manifest["freeze_id"] = freeze_id
with open(OUT / "freeze_manifest_v4.json", "w") as fh:
    json.dump(manifest, fh, indent=2)

# ---------- freeze report ----------
rep = []
ok = True
pri_excl = {l["accession"] for l in ledger if l["cohort"] == "revision_primary_v4"
            and l["reason"] != "ANI_membership_loss_type_anchor"}
conserved = (len(pool_v4) - len(pri_excl) == len(revision_primary)) and not (pool_v4 - pri_excl - set(revision_primary))
rep.append(f"primary arithmetic: pool_v4={len(pool_v4)} - exclusions={len(pri_excl)} = {len(revision_primary)} conserved={conserved} "
           f"(32 historical-only non-members logged separately as ANI_membership_loss_type_anchor)")
ok &= conserved
exp_excl = {l["accession"] for l in ledger if l["cohort"] == "expanded_reference_sensitivity_v4"}
cons2 = (len(pool_exp) - len(exp_excl) == len(expanded)) and not (pool_exp - exp_excl - set(expanded))
rep.append(f"expanded arithmetic: pool_exp={len(pool_exp)} - exclusions={len(exp_excl)} = {len(expanded)} conserved={cons2}")
ok &= cons2
rep.append(f"units among primary unique: {len({unit_of.get(a) for a in revision_primary}) == len([a for a in revision_primary if unit_of.get(a)])}")
rep.append(f"knockout GCA_052073905.1 in any natural cohort: "
           f"{any('GCA_052073905.1' in v for k, v in cohorts.items() if k != 'sg0010_control')}")
rep.append(f"fasta missing for union+expanded: {missing_fasta}")
rep.append(f"freeze_id: {freeze_id}")
rep.append(f"cohort sizes: {json.dumps(manifest['cohort_sizes'])}")
rep.append(f"exclusions by reason: {dict(Counter(l['reason'] for l in ledger))}")
ok &= not missing_fasta
ok &= not any('GCA_052073905.1' in v for k, v in cohorts.items() if k != 'sg0010_control')
with open(OUT / "cohort_freeze_v4.md", "w") as fh:
    fh.write("# cohort freeze v4 (Phase A)\n\n" + "\n".join(f"- {x}" for x in rep)
             + f"\n\nALL_CHECKS: {'PASS' if ok else 'FAIL'}\n")
print("\n".join(rep))
print("ALL_CHECKS:", "PASS" if ok else "FAIL")
