#!/usr/bin/env python3
"""A3: provenance table for all historical + new records.

Participant namespace rule (task doc A3): study_id:participant_id; missing -> unknown.
SG0010 derivative handling: experimental derivatives (GM-mut/Cef/Mock strains, 2023,
PRJNA955008/977800/978365) flagged derivative=experimental; natural baselines
SAMN18320762 / SAMN35359261 are same-strain repeat sequencing of SG0010.

Output: 02_independence/provenance.tsv
"""
import csv
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent.parent
OUT = ROOT / "results_revision_20261001" / "02_independence"
OUT.mkdir(parents=True, exist_ok=True)
META189 = ROOT / "neisseria_subflava_associated_group_metadata.tsv"

SG0010_BASELINES = {"SAMN18320762", "SAMN35359261"}
DERIVATIVE_BIOPROJECTS = {"PRJNA955008", "PRJNA977800", "PRJNA978365"}
DERIVATIVE_STRAIN_PAT = re.compile(r"^(GM-mut-|Cef-c-|Cef-|Mock)")


def attr_json(r):
    try:
        return json.loads(r.get("biosample_attributes_json") or "{}")
    except Exception:
        return {}


def sample_ids(r):
    try:
        return json.loads(r.get("biosample_sample_ids_json") or "[]")
    except Exception:
        return []


def lineage_titles(r):
    try:
        lin = json.loads(r.get("bioproject_lineage_json") or "[]")
    except Exception:
        return []
    titles = []
    for node in lin:
        for bp in node.get("bioprojects", []):
            if bp.get("title"):
                titles.append(bp["title"])
    return titles


def infer_participant(r):
    """Return (participant_id, method, confidence). Conservative: never guess across studies."""
    bs = r["biosample_accession"]
    strain = (r.get("biosample_strain") or "").strip()
    iso = (r.get("biosample_isolate") or "").strip()
    sname = (r.get("biosample_sample_name") or "").strip()
    bp = r["bioproject_accession"]
    attrs = attr_json(r)

    # SG0010 baselines: same strain, same collection year 2017, two sequencing efforts
    if bs in SG0010_BASELINES or strain == "SG0010":
        return ("SG0010", "strain_name", "high")

    # PRJNA714914 Singapore bronchiectasis: isolate IDs HP/TT/CG + 4 digits; participant
    # identity within isolate ID is study-specific knowledge -> unknown unless confirmed
    if bp == "PRJNA714914":
        if iso:
            return (iso, "isolate_id", "unknown_participant_mapping")
        return ("unknown", "none", "low")

    # Pathogen Detection umbrella: isolate=MO## and strain=P####S### are participant codes
    # within this study namespace (EnteroBase Neisseria carriage studies).
    if bp == "PRJNA514245":
        m = re.fullmatch(r"MO\d+", iso)
        if m:
            return (m.group(0), "isolate_participant_code", "high")
        m2 = re.match(r"^(P\d{4})S\d{3}$", strain)
        if m2:
            return (m2.group(1), "strain_participant_code", "high")
        return ("unknown", "umbrella_project", "low")

    # Maiden diversity / Diverse Set: no participant concept
    if bp in ("PRJEB2090", "PRJNA428750"):
        return ("not_applicable_collection", "study_design", "high")

    # named type strains / culture collections: not subject-derived
    if strain or sname:
        return ("not_applicable_culture", "strain_name", "high")

    return ("unknown", "none", "low")


rows = list(csv.DictReader(open(META189, encoding="utf-8-sig"), delimiter="\t"))
out = []
for r in rows:
    bs = r["biosample_accession"]
    strain = (r.get("biosample_strain") or "").strip()
    participant, method, conf = infer_participant(r)
    deriv = ""
    if r["bioproject_accession"] in DERIVATIVE_BIOPROJECTS or DERIVATIVE_STRAIN_PAT.match(strain):
        deriv = "experimental"
    if bs in SG0010_BASELINES:
        deriv = "natural_baseline_repeat"
    sid_list = ";".join(f"{x.get('label','')}:{x.get('value','')}" for x in sample_ids(r))
    out.append({
        "assembly_accession": r["assembly_accession"],
        "biosample_accession": bs,
        "bioproject_accession": r["bioproject_accession"],
        "bioproject_title": r.get("bioproject_title", ""),
        "lineage_titles": " | ".join(lineage_titles(r))[:300],
        "complex_member": r["complex_member"],
        "organism_name": r["organism_name"],
        "strain": strain,
        "isolate": (r.get("biosample_isolate") or "").strip(),
        "sample_name": (r.get("biosample_sample_name") or "").strip(),
        "sample_ids": sid_list[:200],
        "host": r.get("host", ""),
        "isolation_source": r.get("isolation_source", ""),
        "host_disease": r.get("host_disease", ""),
        "collection_date": r.get("collection_date", ""),
        "geo_loc_name": r.get("geo_loc_name", ""),
        "study_namespace": r["bioproject_accession"],
        "participant_id": participant,
        "participant_method": method,
        "participant_confidence": conf,
        "derivative_flag": deriv or "natural",
        "biomaterial_provider": r.get("biomaterial_provider", ""),
        "collected_by": r.get("collected_by", ""),
        "submitter": r.get("submitter", ""),
    })

fields = list(out[0].keys())
with open(OUT / "provenance.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, delimiter="\t", fieldnames=fields)
    w.writeheader()
    w.writerows(out)

from collections import Counter
print("derivative flags:", dict(Counter(o["derivative_flag"] for o in out)))
print("participants:", dict(Counter(o["participant_confidence"] for o in out)))
print("wrote", OUT / "provenance.tsv", len(out), "rows")
