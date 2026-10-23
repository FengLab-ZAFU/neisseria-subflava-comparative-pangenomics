#!/usr/bin/env python3
"""A3: extend provenance.tsv with new in-scope records (not among historical 189).

Participant inference stays conservative: unknown unless strain/culture evidence
or unambiguous same-subject metadata exists. New experimental-derivative patterns
matched by strain name (GM-mut/Cef/Mock).

Appends rows to 02_independence/provenance.tsv (replaces any previous appended rows).
"""
import csv
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent.parent
RT = ROOT / "results_revision_20261001" / "01_retrieval_taxonomy"
IND = ROOT / "results_revision_20261001" / "02_independence"

DERIVATIVE_STRAIN_PAT = re.compile(r"^(GM-mut-|Cef-c-|Cef-|Mock)", re.I)

hist = {r["assembly_accession"] for r in csv.DictReader(
    open(ROOT / "neisseria_subflava_associated_group_metadata.tsv", encoding="utf-8-sig"), delimiter="\t")}
scope = [r["accession"] for r in csv.DictReader(open(RT / "screening_scope_accessions.tsv"), delimiter="\t")]
new_accs = [a for a in scope if a not in hist]

# metadata from frozen JSONL
import json
meta = {}
with open(RT / "raw/genus_482_current_MERGED.jsonl") as fh:
    for line in fh:
        try:
            r = json.loads(line)
        except Exception:
            continue
        if r["accession"] in set(new_accs):
            meta[r["accession"]] = r

base_rows = [r for r in csv.DictReader(open(IND / "provenance.tsv"), delimiter="\t")
             if r["assembly_accession"] in hist]
new_rows = []
for acc in new_accs:
    r = meta.get(acc, {})
    ai = r.get("assembly_info", {})
    bs = ai.get("biosample", {}) or {}
    org = r.get("organism", {}) or {}
    strain = (bs.get("strain") or (org.get("infraspecific_names") or {}).get("strain") or "").strip()
    iso = (bs.get("isolate") or (org.get("infraspecific_names") or {}).get("isolate") or "").strip()
    deriv = "experimental" if DERIVATIVE_STRAIN_PAT.match(strain) else "natural"
    attrs = bs.get("attributes") or []
    attr = {a.get("name"): a.get("value") for a in attrs if isinstance(a, dict)}
    # participant inference within study namespace (same rules as a3_provenance.py)
    participant, pmethod, pconf = "unknown", "none", "low"
    if strain == "SG0010":
        participant, pmethod, pconf = "SG0010", "strain_name", "high"
    if ai.get("bioproject_accession") == "PRJNA514245":
        m = re.fullmatch(r"MO\d+", iso)
        if m:
            participant, pmethod, pconf = m.group(0), "isolate_participant_code", "high"
        else:
            m2 = re.match(r"^(P\d{4})S\d{3}$", strain)
            if m2:
                participant, pmethod, pconf = m2.group(1), "strain_participant_code", "high"
    new_rows.append({
        "assembly_accession": acc,
        "biosample_accession": bs.get("accession", ""),
        "bioproject_accession": ai.get("bioproject_accession", ""),
        "bioproject_title": "",
        "lineage_titles": "",
        "complex_member": "",
        "organism_name": org.get("organism_name", ""),
        "strain": strain,
        "isolate": iso,
        "sample_name": attr.get("sample_name", ""),
        "sample_ids": "",
        "host": bs.get("host", "") or attr.get("host", ""),
        "isolation_source": bs.get("isolation_source", "") or attr.get("isolation_source", ""),
        "host_disease": attr.get("disease", ""),
        "collection_date": bs.get("collection_date", "") or attr.get("collection_date", ""),
        "geo_loc_name": bs.get("geo_loc_name", "") or attr.get("geo_loc_name", ""),
        "study_namespace": ai.get("bioproject_accession", ""),
        "participant_id": participant,
        "participant_method": pmethod,
        "participant_confidence": pconf,
        "derivative_flag": deriv,
        "biomaterial_provider": bs.get("biomaterial_provider", ""),
        "collected_by": attr.get("collected_by", ""),
        "submitter": ai.get("submitter", ""),
    })

fields = list(base_rows[0].keys())
with open(IND / "provenance.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, delimiter="\t", fieldnames=fields)
    w.writeheader()
    w.writerows(base_rows + new_rows)
from collections import Counter
print("provenance rows:", len(base_rows) + len(new_rows), "| new:", len(new_rows))
print("new derivative flags:", dict(Counter(r["derivative_flag"] for r in new_rows)))
