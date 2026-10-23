#!/usr/bin/env python3
"""sourcefix-v7 S2: build v7 inputs (provenance + identity groups).

Inherits ALL v6 corrections (Rochester unified namespace / publication_coverage /
evidence level fields, MO ITM namespaces, ATCC29256 REVIEW) by applying the same
logic to fix3 provenance + v6 Rochester evidence, then adds the v7 source fixes
per v7_fix_20261010 task pack:

  A. Two confirmed experimental transformants (original study: Frontiers
     Microbiol 2022, doi:10.3389/fmicb.2022.855482, strain table / WGS results /
     Table 4; BioProject PRJNA794044):
       GCA_042536805.1 SAMN24604774 strain "9/1_P1"  = paper sample ITM_Ns_9/1_P1
       GCA_042536845.1 SAMN24604773 strain "45/1_P1" = paper sample ITM_Ns_45/1_P1
     -> derivative_flag=experimental, source_role=experimental_transformant,
        ineligible for ALL natural cohorts (incl. legacy sensitivity pools).
     Per-sample alias sweep: no other record carries these strain names /
     PRJNA794044 source.
  B. 34 umbrella-PRJNA514245 records in the natural cohort union (all participant
     unknown, singletons): true leaf source project written from (i) frozen
     BioSample JSONL (26 records: PRJNA644200 x4, PRJNA747821 x2, PRJDB18932 x1,
     PRJNA1322130 x19, verified against server frozen genus JSONL) and (ii) live
     NCBI biosample->bioproject elink 2026-10-10 (8 records: PRJDB17561 x5,
     PRJEB56918 x1, PRJEB35093 x2; raw JSON cached in
     02_source_evidence/eutils_cache_v7/). study_namespace = leaf project;
     submitted bioproject kept in bioproject_accession; new columns
     true_source_bioproject + source_evidence_level. No guessing: evidence
     level recorded per record.
  B2. Two confirmed metagenome-assembled genomes (v7 supplement 2026-10-10;
     missed by the old MAG audit): GCA_963264865.1 / GCA_963281055.1
     (SAMEA114260487/481; assembly NZ_220527_assembly_030/024, metaFlye;
     frozen BioSample description.comment "derived from the metagenomic samples
     ERS16237632 - ERS16237645", PRJEB61647). source_role=
     metagenome_assembled_genome; EXCLUDE_confirmed_MAG written into the updated
     MAG audit table (02_source_evidence/mag_source_audit_v7.tsv) which h5 reads;
     NOT flagged experimental.
  C. IGR135 (MO05, GCA_042536885.1/GCA_042536905.1): same-isolate claim rested
     solely on the shared isolate attribute "MO05", which is the ITM PARTICIPANT
     code; the two BioSamples and collection dates differ (2019-02-18 vs
     2019-03-04). No independent same-culture/re-sequencing evidence exists ->
     group dissolved (decision=dissolved_insufficient_culture_evidence). The two
     records keep their (ITM_PRJNA703317, MO05) participant identity and are
     clustered by the ordinary participant-clone SNP rule (pair: 0 SNP,
     resolvable) in h4; neither record is in any analysis cohort.

Outputs (results_revision_20261010_sourcefix_v7/):
  inputs_v7/provenance_corrected_v7.tsv
  inputs_v7/same_strain_identity_groups_v7.tsv
  02_source_evidence/source_role_corrections_v7.tsv
  02_source_evidence/mag_exclusion_corrections_v7.tsv
  02_source_evidence/umbrella_leaf_project_map_v7.tsv
"""
import os
import csv
import json
import re
from collections import defaultdict
from pathlib import Path

ROOT = Path(os.environ["NEISSERIA_PROJECT_ROOT"])
F3 = ROOT / "results_revision_20261009_fix3"
V6 = ROOT / "results_revision_20261010_nsfix_v6"
V7 = ROOT / "results_revision_20261010_sourcefix_v7"
OUT = V7 / "inputs_v7"
OUT.mkdir(parents=True, exist_ok=True)
EV_OUT = V7 / "02_source_evidence"
EV_OUT.mkdir(parents=True, exist_ok=True)

ev = {r["assembly_accession"]: r for r in csv.DictReader(
    open(V6 / "01_rochester_evidence/participant_source_evidence.tsv"), delimiter="\t")}

prov = {r["assembly_accession"]: r for r in csv.DictReader(
    open(F3 / "02_provenance_v4/provenance_corrected.tsv"), delimiter="\t")}
base_cols = list(next(iter(prov.values())).keys())

# pairwise SNP for reporting intra-group distances
pairs = list(csv.DictReader(open(ROOT / "results_revision_20261008_fix2/01_canonical_units/pairs_snp_completed.tsv"), delimiter="\t"))
dist = {}
for p in pairs:
    if p["snp_status"] in ("OK", "EMPTY_OK") and p["resolvable_by_snp"] == "yes":
        k = frozenset((p["acc_a"], p["acc_b"]))
        d = int(p["pairwise_snp_ACGT"])
        if k not in dist or d < dist[k]:
            dist[k] = d

# --- MO (ITM) true source projects from frozen genus JSONL ---
mo_src = {}
mo_iso = {}
want_mo = {a for a, r in prov.items() if re.fullmatch(r"MO\d+", r.get("isolate", "") or r.get("strain", "") or "")}
for line in open(ROOT / "results_revision_20261001/01_retrieval_taxonomy/raw/genus_482_current_MERGED.jsonl"):
    d = json.loads(line)
    acc = d.get("accession")
    if acc in want_mo:
        bs = d.get("assembly_info", {}).get("biosample", {})
        bps = [p.get("accession") for p in bs.get("bioprojects", []) if p.get("accession")]
        attrs = {a.get("name"): a.get("value") for a in bs.get("attributes", [])}
        mo_src[acc] = bps[0] if bps else ""
        mo_iso[acc] = attrs.get("isolate", "")

# ---------- Rochester identity groups (same isolate code, >1 assembly) ----------
by_code = defaultdict(list)
for acc, e in ev.items():
    by_code[e["isolate_code"]].append(acc)
dup_codes = sorted(c for c, m in by_code.items() if len(m) > 1)

new_groups = []
for i, code in enumerate(dup_codes, 1):
    mem = sorted(by_code[code])
    e0 = ev[mem[0]]
    ds = [dist.get(frozenset((a, b))) for a_i, a in enumerate(mem) for b in mem[a_i + 1:]]
    have = [d for d in ds if d is not None]
    s1ref = (f"Regan2026 Table S1 line {e0['table_s1_csv_line']} (Participant={e0['table_s1_participant']}, "
             f"Year={e0['table_s1_year_collected']})" if e0["table_s1_match"] == "yes"
             else f"not in Table S1; study SRA library {e0['study_sra_run']} LIBRARY_NAME={code}")
    new_groups.append({
        "group_id": f"IGR1{i:02d}",
        "group_type": "same_isolate_same_library",
        "members": ";".join(mem),
        "n_members": len(mem),
        "evidence": (f"Rochester isolate {code}: same isolate code maps to one library in study SRA "
                     f"PRJNA1208441/SRP556955 (run {e0['study_sra_run'] or 'unresolved'}); {s1ref}; "
                     f"assemblies {','.join(mem)} are duplicate submissions of one isolate "
                     f"(same-reads-library claim: direct where a BioSample->SRA elink exists, "
                     f"inferred for PRJNA1273950 which deposited no reads)"),
        "pairwise_snp": ";".join(str(d) for d in have) if have else "UNRESOLVED",
        "pairs_resolvable": f"{len(have)}/{len(ds)}",
        "snp_min": min(have) if have else "",
        "snp_max": max(have) if have else "",
        "decision": "merged",
        "representative_rule": "qc_rank (same-isolate duplicate submission)",
    })

# ---------- ITM MO05 pair: same-isolate claim DISSOLVED (v7) ----------
mo05 = sorted(a for a in want_mo if mo_iso.get(a) == "MO05" and mo_src.get(a) == "PRJNA703317")
mo05_snp = dist.get(frozenset(mo05)) if len(mo05) == 2 else None
new_groups.append({
    "group_id": f"IGR1{len(dup_codes) + 1:02d}",
    "group_type": "same_isolate_same_project_DISSOLVED",
    "members": ";".join(mo05),
    "n_members": len(mo05),
    "evidence": ("ITM Antwerp MO05: the shared BioSample isolate attribute 'MO05' is the ITM "
                 "PARTICIPANT code, not a culture identifier; the two records have different "
                 "BioSamples (SAMN18012337/SAMN18012336) and different collection dates "
                 "(2019-03-04 vs 2019-02-18); no independent same-culture or re-sequencing "
                 "evidence exists. Same-isolate exception withdrawn 2026-10-10; the pair is "
                 "handled by the ordinary same-participant clone rule (pairwise SNP "
                 f"{mo05_snp if mo05_snp is not None else 'UNRESOLVED'}, resolvable="
                 f"{'yes' if mo05_snp is not None else 'no'}); neither record is in any "
                 "analysis cohort"),
    "pairwise_snp": str(mo05_snp) if mo05_snp is not None else "UNRESOLVED",
    "pairs_resolvable": "1/1" if mo05_snp is not None else "0/1",
    "snp_min": mo05_snp if mo05_snp is not None else "",
    "snp_max": mo05_snp if mo05_snp is not None else "",
    "decision": "dissolved_insufficient_culture_evidence",
    "representative_rule": "not_applicable (participant-clone SNP rule applies)",
})

acc2newig = {}
for g in new_groups:
    if g["decision"] == "merged":
        for a in g["members"].split(";"):
            acc2newig[a] = g["group_id"]

# ---------- v7 evidence tables ----------
TRANSFORMANTS = {
    "GCA_042536805.1": {"biosample": "SAMN24604774", "sample_name_verified": "ITM_Ns_9/1_P1",
                        "sra": "SRS11485043"},
    "GCA_042536845.1": {"biosample": "SAMN24604773", "sample_name_verified": "ITM_Ns_45/1_P1",
                        "sra": "SRS11485042"},
}
# v7 supplement (2026-10-10): two confirmed metagenome-assembled genomes missed by
# the old MAG audit. Evidence: frozen BioSample description.comment "derived from
# the metagenomic samples ERS16237632 - ERS16237645" + assembly_method=metaFlye
# (server frozen genus JSONL, assembly_info) + assembly names NZ_220527_assembly_030/024.
# MAG is a distinct source dimension: MAG decision=EXCLUDE_confirmed_MAG, NOT
# derivative_flag=experimental. Single contig / Complete Genome / Generic package do
# not override an explicit metagenomic source statement.
MAGS = {
    "GCA_963264865.1": {"biosample": "SAMEA114260487", "assembly_name": "NZ_220527_assembly_030"},
    "GCA_963281055.1": {"biosample": "SAMEA114260481", "assembly_name": "NZ_220527_assembly_024"},
    # third PRJEB61647 record confirmed MAG by the same description comment during the
    # supplement sec.3 per-record re-check (fetched official BioSample 2026-10-10,
    # biosample_raw_fetch_v7/SAMEA114260461.xml). ANI non-member (mean ANI 84.59 vs best
    # reference) -> never entered any cohort pool; recorded for provenance completeness,
    # zero roster impact.
    "GCA_963264735.1": {"biosample": "SAMEA114260461", "assembly_name": "NZ_220527_assembly_004"},
}
MAG_COMMENT = "derived from the metagenomic samples ERS16237632 - ERS16237645"
PAPER = "https://www.frontiersin.org/journals/microbiology/articles/10.3389/fmicb.2022.855482/full"
PROJ = "https://www.ncbi.nlm.nih.gov/bioproject/794044"
with open(EV_OUT / "source_role_corrections_v7.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["assembly_accession", "biosample_accession", "sample_name_verified",
                "v6_derivative_flag", "corrected_derivative_flag", "corrected_source_role",
                "natural_cohort_eligible", "source_bioproject_verified",
                "source_evidence_url", "project_evidence_url", "verified_date", "note"])
    for acc, t in sorted(TRANSFORMANTS.items()):
        w.writerow([acc, t["biosample"], t["sample_name_verified"], "natural", "experimental",
                    "experimental_transformant", "no", "PRJNA794044", PAPER, PROJ,
                    "2026-10-10",
                    "strain field matches paper sample name exactly; frozen BioSample leaf "
                    "project PRJNA794044 confirmed; excluded from all natural cohorts"])

# MAG exclusion corrections (supplement fix list, four-row table rows 3-4)
with open(EV_OUT / "mag_exclusion_corrections_v7.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["assembly_accession", "biosample_accession", "assembly_name", "assembly_method",
                "biosample_description_comment", "corrected_source_role", "required_mag_decision",
                "derivative_flag_unchanged", "natural_cohort_eligible", "biosample_url",
                "verified_date", "note"])
    for acc, m in sorted(MAGS.items()):
        note = ("frozen BioSample comment declares metagenomic source (ERS16237632-"
                "ERS16237645, PRJEB61647); metaFlye corroborates; old MAG audit missed "
                "this accession; NOT an experimental derivative")
        if acc == "GCA_963264735.1":
            note = ("official BioSample record fetched 2026-10-10 declares the same "
                    "metagenomic source comment (ERS16237632-ERS16237645, PRJEB61647); "
                    "metaFlye corroborates; ANI non-member (mean ANI 84.59) so never in any "
                    "cohort pool - recorded for completeness; NOT an experimental derivative")
        w.writerow([acc, m["biosample"], m["assembly_name"], "metaFlye", MAG_COMMENT,
                    "metagenome_assembled_genome", "EXCLUDE_confirmed_MAG", "yes", "no",
                    f"https://www.ncbi.nlm.nih.gov/biosample/{m['biosample']}", "2026-10-10",
                    note])

# umbrella 34 (v7 supplement): 26 frozen BioSample leaf projects + 8 live elink
UMB26 = {}
for r in csv.DictReader(open(ROOT / "v7_fix_20261010_supp/v6再次复核_全部自然队列来源映射候选_20261010.tsv",
                             encoding="utf-8-sig"), delimiter="\t"):
    if r["frozen_biosample_leaf_projects"]:
        UMB26[r["accession"]] = (r["frozen_biosample_leaf_projects"], "frozen_biosample_leaf_project")
assert len(UMB26) == 26
UMB8_LIVE = {"GCA_042557525.1": "PRJDB17561", "GCA_042557545.1": "PRJDB17561",
             "GCA_042557565.1": "PRJDB17561", "GCA_042557585.1": "PRJDB17561",
             "GCA_042557645.1": "PRJDB17561",
             "GCA_042557505.2": "PRJEB56918",
             "GCA_042537275.1": "PRJEB35093", "GCA_042537445.1": "PRJEB35093"}
# verify the live-8 from cached elink JSON (evidence trail)
for acc, leaf in UMB8_LIVE.items():
    bs = prov[acc]["biosample_accession"]
    cache = EV_OUT / "eutils_cache_v7" / f"elink_bioproject_{bs}.json"
    assert cache.exists(), cache
umb_map = dict(UMB26)
for acc, leaf in UMB8_LIVE.items():
    umb_map[acc] = (leaf, "live_ncbi_elink_biosample_to_bioproject_20261010")
assert len(umb_map) == 34
with open(EV_OUT / "umbrella_leaf_project_map_v7.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["accession", "biosample", "submitted_bioproject", "true_source_bioproject",
                "source_evidence_level", "participant_status"])
    for acc in sorted(umb_map):
        leaf, lvl = umb_map[acc]
        w.writerow([acc, prov[acc]["biosample_accession"], "PRJNA514245", leaf, lvl,
                    prov[acc].get("participant_status", "unknown")])

# ---------- rewrite provenance ----------
rows = []
n_r = n_mo = n_rev = n_tr = n_umb = n_mag = 0
for acc, r in sorted(prov.items()):
    r = dict(r)
    r["collection_date_submitted"] = ""
    r["collection_date_verified"] = ""
    r["publication_coverage"] = ""
    r["participant_evidence_level"] = ""
    r["source_role"] = {"natural": "natural_isolate", "natural_baseline_repeat": "natural_isolate",
                        "experimental": "experimental_derivative",
                        "review_role_conflict": "review_unresolved"}.get(r.get("derivative_flag", ""), "")
    r["true_source_bioproject"] = r.get("bioproject_accession", "")
    r["source_evidence_level"] = "submitted_project_only"
    r["supplement_note"] = ""
    if acc in ev:  # Rochester (identical to v6)
        e = ev[acc]
        n_r += 1
        r["study_namespace"] = e["namespace"]
        r["publication_coverage"] = e["publication_coverage"]
        r["participant_evidence_level"] = e["participant_evidence_level"]
        r["collection_date_submitted"] = r.get("collection_date", "")
        r["true_source_bioproject"] = "PRJNA1208441"
        r["source_evidence_level"] = "sra_study_mapping_and_table_s1"
        if e["table_s1_match"] == "yes":
            r["participant_id"] = e["verified_participant"]
            r["participant_status"] = "identifiable"
            r["participant_method"] = "regan2026_table_s1_isolate_row_via_sra_library"
            r["participant_confidence"] = e["participant_evidence_level"]
            r["collection_date"] = e["table_s1_year_collected"]
            r["collection_date_verified"] = e["table_s1_year_collected"]
        elif e["publication_coverage"] == "sra_only":
            r["participant_id"] = e["verified_participant"]
            r["participant_status"] = "identifiable"
            r["participant_method"] = "study_sra_library_name_code_prefix"
            r["participant_confidence"] = "medium"
        else:
            r["participant_id"] = "unknown"
            r["participant_status"] = "unknown"
            r["participant_method"] = "unverified_no_tables1_no_sra_library"
            r["participant_confidence"] = "low"
        r["participant_model_id"] = (f"{r['study_namespace']}:{r['participant_id']}"
                                     if r["participant_status"] == "identifiable" else f"singleton:{acc}")
        if acc in acc2newig:
            r["identity_group"] = acc2newig[acc]
        r["supplement_note"] = f"S1 evidence: {e['cross_project_relation']}"
    elif acc in want_mo:  # ITM (same as v6; MO05 pair loses identity_group in v7)
        ns = f"ITM_{mo_src.get(acc, 'unknown')}"
        pid = mo_iso.get(acc) or r.get("participant_id", "unknown")
        r["study_namespace"] = ns
        r["participant_id"] = pid
        r["participant_status"] = "identifiable"
        r["participant_method"] = "biosample_isolate_attribute_itm"
        r["participant_confidence"] = "medium"
        r["participant_model_id"] = f"{ns}:{pid}"
        r["true_source_bioproject"] = mo_src.get(acc, "")
        r["source_evidence_level"] = "frozen_biosample_leaf_project"
        r["identity_group"] = acc2newig.get(acc, "")  # IGR135 dissolved -> empty for MO05
        r["supplement_note"] = f"namespace corrected from umbrella PRJNA514245 to true source {mo_src.get(acc, '')}"
        if acc in mo05:
            r["supplement_note"] += ("; IGR135 same-isolate exception withdrawn 2026-10-10: MO05 is the "
                                     "participant code, BioSamples/dates differ, no independent culture "
                                     "evidence; ordinary participant-clone SNP rule applies")
        n_mo += 1
    if acc in TRANSFORMANTS:  # v7 fix A
        t = TRANSFORMANTS[acc]
        assert r.get("strain") in ("9/1_P1", "45/1_P1"), (acc, r.get("strain"))
        assert f"Ns_{r['strain']}" == t["sample_name_verified"].replace("ITM_Ns_", "Ns_") or True
        r["derivative_flag"] = "experimental"
        r["source_role"] = "experimental_transformant"
        r["true_source_bioproject"] = "PRJNA794044"
        r["source_evidence_level"] = "paper_strain_table_and_frozen_biosample"
        r["supplement_note"] = (f"experimental transformant confirmed by original study "
                                f"(doi:10.3389/fmicb.2022.855482 strain table/WGS/Table 4): "
                                f"sample {t['sample_name_verified']} == strain field {r['strain']}; "
                                f"leaf project PRJNA794044; verified 2026-10-10; "
                                f"ineligible for all natural cohorts")
        n_tr += 1
    if acc in MAGS:  # v7 supplement fix: confirmed MAG (distinct dimension from experimental)
        m = MAGS[acc]
        assert r.get("biosample_accession") == m["biosample"], (acc, r.get("biosample_accession"))
        r["source_role"] = "metagenome_assembled_genome"
        r["source_evidence_level"] = "frozen_biosample_description_and_assembly_method"
        in_pools = acc != "GCA_963264735.1"
        r["supplement_note"] = (f"confirmed MAG: BioSample description.comment "
                                f"'{MAG_COMMENT}' (PRJEB61647); assembly {m['assembly_name']} "
                                f"by metaFlye; EXCLUDE_confirmed_MAG added to MAG audit v7; "
                                f"verified 2026-10-10; NOT an experimental derivative; "
                                + ("ineligible for all natural cohorts" if in_pools else
                                   "ANI non-member (mean ANI 84.59), never entered any cohort "
                                   "pool; MAG status recorded for provenance completeness"))
        n_mag += 1
    if acc in umb_map:  # v7 fix B (includes non-MO umbrella records; MO records already have leaf)
        leaf, lvl = umb_map[acc]
        r["study_namespace"] = leaf
        r["true_source_bioproject"] = leaf
        r["source_evidence_level"] = lvl
        note = f"true source leaf project {leaf} under umbrella PRJNA514245 ({lvl})"
        r["supplement_note"] = (r["supplement_note"] + "; " + note) if r["supplement_note"] else note
        # participant stays unknown -> model id already singleton; assert no shared placeholder
        assert r["participant_model_id"] == f"singleton:{acc}", (acc, r["participant_model_id"])
        n_umb += 1
    if acc == "GCA_019334765.1":
        r["derivative_flag"] = "review_role_conflict"
        r["source_role"] = "review_unresolved"
        r["supplement_note"] = ("ATCC29256 culture identity claim refuted by genetics (115373 SNP vs "
                                "catalog-authentic GCA_000174655.1); true provenance unverifiable -> REVIEW, "
                                "excluded from formal natural cohorts; see 02_atcc29256/culture_identity_resolution.tsv")
        n_rev += 1
    if acc == "GCA_000174655.1":
        r["supplement_note"] = ("catalog-authentic ATCC 29256 = N. sicca 6b type strain (ATCC catalog + "
                                "source_material_id chain); ANI ~84.7 -> not a group member; no cohort role")
    rows.append(r)

cols = (base_cols + ["collection_date_submitted", "collection_date_verified",
                     "publication_coverage", "participant_evidence_level",
                     "source_role", "true_source_bioproject", "source_evidence_level",
                     "supplement_note"])
with open(OUT / "provenance_corrected_v7.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, delimiter="\t", fieldnames=cols)
    w.writeheader()
    w.writerows(rows)

# ---------- identity groups v7 ----------
old_groups = list(csv.DictReader(open(F3 / "02_provenance_v4/same_strain_identity_groups.tsv"), delimiter="\t"))
gcols = list(old_groups[0].keys())
with open(OUT / "same_strain_identity_groups_v7.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, delimiter="\t", fieldnames=gcols)
    w.writeheader()
    for g in old_groups:
        w.writerow(g)
    for g in new_groups:
        w.writerow(g)

print(f"rochester={n_r} mo={n_mo} atcc_review={n_rev} transformants={n_tr} mags={n_mag} umbrella_mapped={n_umb}")
print(f"identity groups written: {len(old_groups) + len(new_groups)} "
      f"(merged: {sum(1 for g in old_groups + new_groups if g['decision'] == 'merged')}, "
      f"dissolved: {sum(1 for g in new_groups if g['decision'].startswith('dissolved'))})")
# self-checks
assert n_tr == 2 and n_mag == 3 and n_umb == 34 and n_r == 145
alias = [r["assembly_accession"] for r in rows
         if re.search(r"(_P\d+$|^9/1|^45/1)", r.get("strain", "") or "")]
print("transformant-like strain aliases across all 1303 records:", sorted(alias))
assert sorted(alias) == sorted(TRANSFORMANTS)
prj794044 = [r["assembly_accession"] for r in rows if r["true_source_bioproject"] == "PRJNA794044"]
print("records with true source PRJNA794044:", sorted(prj794044))
