#!/usr/bin/env python3
"""fix3-supp S1a: map every P-coded screen record's BioSample to its SRA
run(s) via NCBI eutils, so Rochester cross-project identity can be decided on
library-level evidence (Regan 2026 Table S1 SRR column) instead of codes+SNP.

Inputs: provenance_corrected.tsv (fix3), frozen genus JSONL (biosample acc).
Output: results_revision_20261009_fix3_phaseA_supplement/01_rochester_evidence/biosample_sra_map.tsv
        + raw esummary JSON cached per biosample (evidence archive)
"""
import csv
import json
import os
import time
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(os.environ["NEISSERIA_PROJECT_ROOT"])
OUT = ROOT / "results_revision_20261009_fix3_phaseA_supplement/01_rochester_evidence"
CACHE = OUT / "eutils_cache"
CACHE.mkdir(parents=True, exist_ok=True)
KEY = os.environ.get("NCBI_API_KEY", "")
BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

# target records: P-coded strains under PRJNA514245 umbrella or PRJNA1273950
prov = list(csv.DictReader(open(ROOT / "results_revision_20261009_fix3/02_provenance_v4/provenance_corrected.tsv"), delimiter="\t"))
targets = {}
import re
for r in prov:
    s = r.get("strain", "")
    if re.fullmatch(r"P\d{4}[.]?S\d{3}", s or "") and r["bioproject_accession"] in ("PRJNA514245", "PRJNA1273950"):
        targets[r["biosample_accession"]] = r["assembly_accession"]
print("target biosamples:", len(targets))

def get(url, retries=4):
    for i in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=60) as fh:
                return fh.read().decode("utf-8", "replace")
        except Exception as e:
            time.sleep(2 * (i + 1))
    raise RuntimeError(f"failed: {url}")

def cached(name, fetch):
    p = CACHE / name
    if p.exists():
        return p.read_text()
    txt = fetch()
    p.write_text(txt)
    time.sleep(0.12 if KEY else 0.4)
    return txt

rows = []
for i, (bs, acc) in enumerate(sorted(targets.items())):
    # 1) biosample uid
    es = cached(f"{bs}_esearch.json", lambda: get(
        f"{BASE}/esearch.fcgi?db=biosample&term={bs}%5Baccn%5D&retmode=json"
        + (f"&api_key={KEY}" if KEY else "")))
    ids = json.loads(es).get("esearchresult", {}).get("idlist", [])
    if not ids:
        rows.append({"assembly_accession": acc, "biosample": bs, "biosample_uid": "",
                     "sra_runs": "", "sra_studies": "", "note": "biosample not found"})
        continue
    uid = ids[0]
    # 2) elink biosample -> sra
    el = cached(f"{bs}_elink.json", lambda: get(
        f"{BASE}/elink.fcgi?dbfrom=biosample&db=sra&id={uid}&retmode=json"
        + (f"&api_key={KEY}" if KEY else "")))
    linksets = json.loads(el).get("linksets", [])
    sra_ids = []
    if linksets and linksets[0].get("linksetdbs"):
        for ldb in linksets[0]["linksetdbs"]:
            if ldb.get("dbto") == "sra":
                sra_ids.extend(ldb.get("links", []))
    runs, studies = [], []
    if sra_ids:
        # 3) esummary sra
        su = cached(f"{bs}_sra_sum.json", lambda: get(
            f"{BASE}/esummary.fcgi?db=sra&id={','.join(sra_ids)}&retmode=json"
            + (f"&api_key={KEY}" if KEY else "")))
        sj = json.loads(su).get("result", {})
        for sid in sj.get("uids", []):
            doc = sj.get(sid, {})
            expxml = doc.get("expxml", "")
            runs_xml = doc.get("runs", "")
            runs += re.findall(r'acc="(SRR\d+)"', runs_xml)
            studies += re.findall(r'study_acc="(SRP\d+|PRJNA\d+|ERP\d+|PRJEB\d+)"', runs_xml)
            if not studies:
                studies += re.findall(r'<Study[^>]*acc="([^"]+)"', expxml)
    rows.append({"assembly_accession": acc, "biosample": bs, "biosample_uid": uid,
                 "sra_runs": ";".join(sorted(set(runs))), "sra_studies": ";".join(sorted(set(studies))),
                 "note": ""})
    if (i + 1) % 25 == 0:
        print(f"  {i+1}/{len(targets)}")

with open(OUT / "biosample_sra_map.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, delimiter="\t", fieldnames=list(rows[0].keys()))
    w.writeheader()
    w.writerows(rows)
n_ok = sum(1 for r in rows if r["sra_runs"])
print(f"mapped with SRA runs: {n_ok}/{len(rows)}")
print("run counts distribution:", {})
from collections import Counter
print(Counter(len(r["sra_runs"].split(";")) if r["sra_runs"] else 0 for r in rows))
