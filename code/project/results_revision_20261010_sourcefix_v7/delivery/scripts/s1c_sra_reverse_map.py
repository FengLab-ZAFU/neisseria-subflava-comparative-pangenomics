#!/usr/bin/env python3
"""fix3-supp S1c: reverse map SRA runs -> BioSample for the Rochester source
projects (PRJNA1273950 whole-study 2025 submission; PRJNA1208441 2024 batch),
because biosample->sra elink returns no links for the PRJNA1273950 records.

Output: 01_rochester_evidence/sra_run_biosample_map.tsv (run, biosample, bioproject, study)
Raw esearch/esummary JSON cached in eutils_cache/.
"""
import json
import os
import re
import time
import urllib.request
from pathlib import Path

ROOT = Path(os.environ["NEISSERIA_PROJECT_ROOT"])
OUT = ROOT / "results_revision_20261009_fix3_phaseA_supplement/01_rochester_evidence"
CACHE = OUT / "eutils_cache"
CACHE.mkdir(parents=True, exist_ok=True)
KEY = os.environ.get("NCBI_API_KEY", "")
BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

def get(url, retries=4):
    for i in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=90) as fh:
                return fh.read().decode("utf-8", "replace")
        except Exception:
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
for proj in ("PRJNA1273950", "PRJNA1208441"):
    es = cached(f"{proj}_sra_esearch.json", lambda: get(
        f"{BASE}/esearch.fcgi?db=sra&term={proj}%5BAll%20Fields%5D&retmax=10000&retmode=json"
        + (f"&api_key={KEY}" if KEY else "")))
    ids = json.loads(es).get("esearchresult", {}).get("idlist", [])
    print(proj, "sra uids:", len(ids))
    for i in range(0, len(ids), 200):
        chunk = ids[i:i + 200]
        su = cached(f"{proj}_sra_sum_{i}.json", lambda: get(
            f"{BASE}/esummary.fcgi?db=sra&id={','.join(chunk)}&retmode=json"
            + (f"&api_key={KEY}" if KEY else "")))
        sj = json.loads(su).get("result", {})
        for sid in sj.get("uids", []):
            doc = sj.get(sid, {})
            expxml = doc.get("expxml", "")
            runs_xml = doc.get("runs", "")
            bs = re.findall(r"<Biosample>(SAMN\d+)</Biosample>", expxml, re.I) or re.findall(r'sample_acc="(SAMN\d+)"', expxml)
            runs = re.findall(r'acc="(SRR\d+)"', runs_xml)
            study = re.findall(r'study_acc="(\w+)"', runs_xml)
            for r in runs:
                rows.append({"bioproject_query": proj, "sra_uid": sid, "run": r,
                             "biosample": bs[0] if bs else "", "study": study[0] if study else ""})

import csv
with open(OUT / "sra_run_biosample_map.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, delimiter="\t", fieldnames=["bioproject_query", "sra_uid", "run", "biosample", "study"])
    w.writeheader()
    w.writerows(rows)
print("rows:", len(rows))
from collections import Counter
print(Counter(r["bioproject_query"] for r in rows))
missing_bs = sum(1 for r in rows if not r["biosample"])
print("rows missing biosample:", missing_bs)
