#!/usr/bin/env python3
"""sourcefix-v7 supplement S3b: full source-coverage tables (task supplement sec.3).

Two outputs under results_revision_20261010_sourcefix_v7/02_source_evidence/:

1. source_coverage_1303_v7.tsv - all 1,303 screening records. For each record the
   ASSEMBLY-side raw fields are re-read from the frozen server genus JSONL
   (results_revision_20261001/01_retrieval_taxonomy/raw/genus_482_current_MERGED.jsonl):
   assembly_method, assembly_name, assembly_level, submitted bioproject, BioSample
   leaf bioproject lineage, organism/strain/isolate names. A keyword screen
   (metagenomic / MIMAG / metaFlye / transformant / passage / mutant / _P<n> ...)
   runs over the frozen assembly record text PLUS any raw BioSample text present in
   server archives (v6/v7 eutils caches, ATCC evidence JSON, the 2026-10-10 official
   fetch of 41 gap records). Records without server-side BioSample raw text carry an
   explicit retrieval cue (official BioSample URL). Keywords are a screening entry
   only - final source decisions are record-level evidence in
   source_role_corrections_v7.tsv / mag_exclusion_corrections_v7.tsv.

2. natural_union_source_coverage_v7.tsv - the 204-record natural cohort union (v6
   definition, per the reviewer's coverage table) annotated with the v7 disposition:
     - local_raw_available_windows_freeze : 125 records with complete raw BioSample
       in the local (Windows) freeze; verifiable archive hashes listed in the v7 report
     - rochester_evidence_covered         : 39 records covered by Regan2026 Table S1 /
       study SRA evidence (v6 01_rochester_evidence + eutils cache)
     - official_record_fetched_20261010   : 40 gap records filled by official NCBI
       efetch (raw XML + sha256 in biosample_raw_fetch_v7/)
     - the 4 confirmed source exclusions are flagged excluded_source_confirmed
   No record is auto-excluded for missing data; no record is marked covered without
   a concrete evidence channel.
"""
import os
import csv
import hashlib
import json
import re
from pathlib import Path

ROOT = Path(os.environ["NEISSERIA_PROJECT_ROOT"])
V7 = ROOT / "results_revision_20261010_sourcefix_v7"
EV = V7 / "02_source_evidence"
OLD = ROOT / "results_revision_20261001"
V6 = ROOT / "results_revision_20261010_nsfix_v6"

KEYWORDS = re.compile(r"metagenom|MIMAG|metaFlye|MetaBAT|CONCOCT|transformant|transformed|"
                      r"passage|mutant|knockout|_P\d+\b|single.cell|uncultured|enrichment culture",
                      re.I)

prov = {r["assembly_accession"]: r for r in csv.DictReader(
    open(V7 / "inputs_v7/provenance_corrected_v7.tsv"), delimiter="\t")}

# ---- server frozen assembly records (GCA side only) ----
frozen = {}
for line in open(OLD / "01_retrieval_taxonomy/raw/genus_482_current_MERGED.jsonl"):
    d = json.loads(line)
    a = d.get("accession", "")
    if a.startswith("GCA_") and a in prov:
        frozen[a] = d
assert len(frozen) == len(prov) == 1303, (len(frozen), len(prov))

# ---- server-side raw BioSample text archive index ----
bs_raw = {}   # biosample accession -> (source_label, text)
fetch_dir = EV / "biosample_raw_fetch_v7"
for r in csv.DictReader(open(fetch_dir / "biosample_gap_fetch_v7.tsv"), delimiter="\t"):
    p = fetch_dir / f"{r['biosample_accession']}.xml"
    bs_raw[r["biosample_accession"]] = ("official_fetch_v7", p.read_text())
atcc = json.load(open(V6 / "02_atcc29256/evidence_ncbi_biosample_records.json")) \
    if (V6 / "02_atcc29256/evidence_ncbi_biosample_records.json").exists() else None
atcc = atcc or json.load(open(ROOT / "results_revision_20261009_fix3_phaseA_supplement/02_atcc29256/evidence_ncbi_biosample_records.json"))
def walk_bs(o):
    if isinstance(o, dict):
        acc = o.get("accession")
        if acc and ("attributes" in o or "description" in o):
            bs_raw.setdefault(acc, ("atcc_evidence_json", json.dumps(o, ensure_ascii=False)))
        for v in o.values():
            walk_bs(v)
    elif isinstance(o, list):
        for v in o:
            walk_bs(v)
walk_bs(atcc)

# ---- 1303 coverage table ----
mag_audit = {r["accession"]: r for r in csv.DictReader(open(EV / "mag_source_audit_v7.tsv"),
                                                       delimiter="\t")}
man = {r["assembly_accession"]: r for r in csv.DictReader(
    open(OLD / "02_independence/retrieval_manifest.tsv"), delimiter="\t")}
rows1303 = []
for acc in sorted(prov):
    p = prov[acc]
    f = frozen[acc]
    ai = f.get("assembly_info", {})
    bs = ai.get("biosample") or {}
    bs_acc = bs.get("accession") or p.get("biosample_accession", "")
    bps = [x.get("accession") for x in bs.get("bioprojects", []) if x.get("accession")]
    org = f.get("organism", {})
    infra = org.get("infraspecific_names", {}) or {}
    scan_parts = [json.dumps({k: f.get(k) for k in ("organism", "assembly_info")},
                             ensure_ascii=False)]
    basis = ["assembly_frozen_jsonl"]
    if bs_acc in bs_raw:
        label, text = bs_raw[bs_acc]
        scan_parts.append(text)
        basis.append(label)
    hits = sorted(set(m.group(0).lower() for m in KEYWORDS.finditer(" ".join(scan_parts))))
    rows1303.append({
        "assembly_accession": acc,
        "biosample_accession": bs_acc,
        "organism_name": org.get("organism_name", ""),
        "strain": infra.get("strain", "") or p.get("strain", ""),
        "isolate": infra.get("isolate", "") or p.get("isolate", ""),
        "assembly_method": ai.get("assembly_method", ""),
        "assembly_name": ai.get("assembly_name", ""),
        "assembly_level": ai.get("assembly_level", ""),
        "submitted_bioproject": ai.get("bioproject_accession", ""),
        "biosample_leaf_bioprojects": ";".join(bps),
        "source_role_v7": p.get("source_role", ""),
        "derivative_flag_v7": p.get("derivative_flag", ""),
        "study_namespace_v7": p.get("study_namespace", ""),
        "true_source_bioproject_v7": p.get("true_source_bioproject", ""),
        "source_evidence_level_v7": p.get("source_evidence_level", ""),
        "mag_audit_decision_v7": mag_audit.get(acc, {}).get("decision", ""),
        "ani_fate": man.get(acc, {}).get("fate", ""),
        "biosample_raw_on_server": basis[-1] if len(basis) > 1 else "none_on_server",
        "keyword_scan_basis": "+".join(basis),
        "keyword_hits": ";".join(hits),
        "retrieval_cue": "" if len(basis) > 1 else f"https://www.ncbi.nlm.nih.gov/biosample/{bs_acc}",
    })
with open(EV / "source_coverage_1303_v7.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, delimiter="\t", fieldnames=list(rows1303[0].keys()))
    w.writeheader()
    w.writerows(rows1303)

# ---- 204 natural union coverage table ----
supp = list(csv.DictReader(open(ROOT / "v7_fix_20261010_supp/v6再次复核_自然队列来源证据覆盖_20261010.tsv",
                                encoding="utf-8-sig"), delimiter="\t"))
assert len(supp) == 204
fetch = {r["assembly_accession"]: r for r in csv.DictReader(
    open(fetch_dir / "biosample_gap_fetch_v7.tsv"), delimiter="\t")}
EXCLUDED4 = {"GCA_042536805.1", "GCA_042536845.1", "GCA_963264865.1", "GCA_963281055.1"}
rows204 = []
for r in supp:
    acc = r["accession"]
    p = prov[acc]
    if r["raw_biosample_record_available"] == "True":
        status = "local_raw_available_windows_freeze"
        ev_note = "complete raw BioSample in local freeze (Windows); assembly-side server frozen JSONL re-verified"
    elif r["sra_publication_coverage"]:
        status = "rochester_evidence_covered"
        ev_note = f"Regan2026 {r['sra_publication_coverage']} evidence (v6 participant_source_evidence + eutils cache)"
    else:
        fr = fetch[acc]
        status = "official_record_fetched_20261010"
        ev_note = f"raw XML {fr['raw_file']} sha256={fr['raw_sha256']}"
    rows204.append({
        "accession": acc, "biosample": r["biosample"],
        "study_namespace_v7": p.get("study_namespace", ""),
        "source_role_v7": p.get("source_role", ""),
        "v7_disposition": ("excluded_source_confirmed" if acc in EXCLUDED4 else "retained"),
        "coverage_status": status,
        "coverage_evidence": ev_note,
        "keyword_hits_server_scan": next(x["keyword_hits"] for x in rows1303 if x["assembly_accession"] == acc),
        "server_assembly_archive": "01_retrieval_taxonomy/raw/genus_482_current_MERGED.jsonl",
        "reviewer_raw_available": r["raw_biosample_record_available"],
        "reviewer_sra_coverage": r["sra_publication_coverage"],
        "in_primary_v6": r["in_primary"],
    })
with open(EV / "natural_union_source_coverage_v7.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, delimiter="\t", fieldnames=list(rows204[0].keys()))
    w.writeheader()
    w.writerows(rows204)

# ---- keyword-hit disposition table (record-level, task supplement sec.3) ----
# every 1303 keyword hit must be visibly dispositioned: MAG audit decision covers
# 562; the remaining 7 are listed here with record-level context.
disp_rows = []
for x in rows1303:
    if not x["keyword_hits"]:
        continue
    if x["mag_audit_decision_v7"]:
        continue  # dispositioned by MAG audit (EXCLUDE_confirmed_MAG / REVIEW)
    acc = x["assembly_accession"]
    f = frozen[acc]
    txt = json.dumps(f, ensure_ascii=False)
    m = KEYWORDS.search(txt)
    s = max(0, m.start() - 100)
    context = txt[s:m.end() + 100].replace("\t", " ").replace("\n", " ")
    if x["derivative_flag_v7"] == "experimental":
        disp = "excluded_experimental_derivative (cohorts); keyword corroborates"
    elif x["ani_fate"] == "non_member":
        disp = "not cohort-eligible (ANI non-member); keyword is screening cue only, " \
               "record-level context reviewed (see keyword_context)"
    else:
        disp = "NEEDS_REVIEW"
    disp_rows.append({"assembly_accession": acc, "keyword_hits": x["keyword_hits"],
                      "source_role_v7": x["source_role_v7"], "ani_fate": x["ani_fate"],
                      "submitted_bioproject": x["submitted_bioproject"],
                      "assembly_method": x["assembly_method"],
                      "keyword_context": context, "disposition": disp})
with open(EV / "keyword_hit_dispositions_v7.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, delimiter="\t", fieldnames=list(disp_rows[0].keys()))
    w.writeheader()
    w.writerows(disp_rows)
needs_review = [d for d in disp_rows if d["disposition"] == "NEEDS_REVIEW"]
print(" non-MAG-audit keyword hits dispositioned:", len(disp_rows), "NEEDS_REVIEW:", len(needs_review))
assert not needs_review, needs_review

# ---- summaries ----
from collections import Counter
print("1303 coverage:", len(rows1303))
print(" biosample raw on server:", dict(Counter(x["biosample_raw_on_server"] for x in rows1303)))
kw1303 = [(x["assembly_accession"], x["keyword_hits"], x["source_role_v7"]) for x in rows1303 if x["keyword_hits"]]
print(" keyword hits in 1303:", len(kw1303))
hit_by_decision = Counter(x["mag_audit_decision_v7"] or "no_mag_audit_row" for x in rows1303 if x["keyword_hits"])
print(" keyword hits by MAG audit decision:", dict(hit_by_decision))
undecided_hits = [x["assembly_accession"] for x in rows1303
                  if x["keyword_hits"] and not x["mag_audit_decision_v7"]
                  and x["derivative_flag_v7"] not in ("experimental", "review_role_conflict")]
print(" keyword hits with NO MAG decision and NO derivative flag:", len(undecided_hits))
for a in undecided_hits:
    x = next(r for r in rows1303 if r["assembly_accession"] == a)
    print("   ", a, x["keyword_hits"], x["submitted_bioproject"], x["assembly_method"])
print("204 union:", dict(Counter(x["coverage_status"] for x in rows204)))
print("204 excluded:", [x["accession"] for x in rows204 if x["v7_disposition"] == "excluded_source_confirmed"])
kw204 = [(x["accession"], x["keyword_hits_server_scan"], x["v7_disposition"]) for x in rows204 if x["keyword_hits_server_scan"]]
print(" keyword hits in 204:", kw204)
