#!/usr/bin/env python3
"""sourcefix-v7 supplement S3a: fetch official BioSample records for the 40
coverage-gap accessions (+ GCA_963264735.1, the third PRJEB61647 record that
needs per-record verification under the mixed-project rule).

Task (supplement sec.3): the 40 records without local raw BioSample and without
Rochester paper/SRA evidence must be filled from server archives first, else
from official BioSample records; raw responses, access date and hashes must be
kept. No record may be auto-excluded for missing data, and no record may be
marked PASS without data.

Server archives were checked first: the download zips under
results_revision_20261001/01_retrieval_taxonomy/zips*/ contain only
assembly_data_report.jsonl (assembly side, no BioSample description), and the
eutils caches (fix3/v6/v7) cover only Rochester + umbrella records. None of the
41 BioSamples is present in any server archive -> official NCBI eutils fetch
(biosample esearch + efetch XML), raw responses saved here with sha256.

Outputs (results_revision_20261010_sourcefix_v7/02_source_evidence/biosample_raw_fetch_v7/):
  esearch_biosample_gap41.json / efetch_biosample_gap41.xml  (raw responses)
  <biosample>.xml                                           (per-record block)
  biosample_gap_fetch_v7.tsv                                (parsed summary)
  sha256 manifest line per file (access date 2026-10-10)
"""
import os
import csv
import hashlib
import json
import re
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path

ROOT = Path(os.environ["NEISSERIA_PROJECT_ROOT"])
V7 = ROOT / "results_revision_20261010_sourcefix_v7"
OUT = V7 / "02_source_evidence/biosample_raw_fetch_v7"
OUT.mkdir(parents=True, exist_ok=True)
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

# 40 gap records from the reviewer's coverage table + the third PRJEB61647 record
rows = list(csv.DictReader(open(ROOT / "v7_fix_20261010_supp/v6再次复核_自然队列来源证据覆盖_20261010.tsv",
                                encoding="utf-8-sig"), delimiter="\t"))
gap = [(r["accession"], r["biosample"]) for r in rows
       if r["raw_biosample_record_available"] == "False" and r["sra_publication_coverage"] == ""]
assert len(gap) == 40
gap.append(("GCA_963264735.1", "SAMEA114260461"))  # mixed-project per-record check (sec.3)
acc_of = {bs: acc for acc, bs in gap}

def get(url, params, dest):
    q = urllib.parse.urlencode(params)
    req = urllib.request.Request(f"{url}?{q}", headers={"User-Agent": "phaseA-sourcefix-v7"})
    with urllib.request.urlopen(req, timeout=60) as r:
        body = r.read()
    dest.write_bytes(body)
    return body

# esearch (batched OR query)
term = " OR ".join(f"{bs}[Accession]" for _, bs in gap)
raw_search = get(f"{EUTILS}/esearch.fcgi",
                 {"db": "biosample", "term": term, "retmode": "json", "retmax": 100},
                 OUT / "esearch_biosample_gap41.json")
time.sleep(0.4)
uids = json.loads(raw_search)["esearchresult"]["idlist"]
assert len(uids) == 41, f"expected 41 uids, got {len(uids)}"

# efetch full XML (one batch)
raw_xml = get(f"{EUTILS}/efetch.fcgi", {"db": "biosample", "id": ",".join(uids),
                                        "rettype": "xml", "retmode": "xml"},
              OUT / "efetch_biosample_gap41.xml")

KEYWORDS = re.compile(r"metagenom|MIMAG|metaFlye|MetaBAT|concoct|transform|passage|mutant|"
                      r"knockout|_P\d+\b|single.cell|uncultured|enrichment", re.I)

tree = ET.fromstring(raw_xml)
out_rows = []
seen = set()
for bs in tree.iter("BioSample"):
    accession = bs.get("accession", "")
    if accession not in acc_of:
        continue
    seen.add(accession)
    desc = bs.find("Description")
    title = (desc.findtext("Title") or "") if desc is not None else ""
    comment = ""
    if desc is not None:
        c = desc.find("Comment")
        if c is not None and c.findtext("Paragraph"):
            comment = " | ".join(p.text or "" for p in c.findall("Paragraph"))
    package = (bs.findtext("Package") or "")
    attrs = {a.get("attribute_name") or a.get("harmonized_name") or "": (a.text or "")
             for a in bs.iter("Attribute")}
    links = [l.get("target") for l in bs.iter("Link") if l.get("target", "").startswith(("PRJ", "PRJDB", "PRJEB"))]
    block = ET.tostring(bs, encoding="unicode")
    block_path = OUT / f"{accession}.xml"
    block_path.write_text(block)
    scan_text = " ".join([title, comment, package] + [f"{k}={v}" for k, v in attrs.items()])
    hits = sorted(set(m.group(0).lower() for m in KEYWORDS.finditer(scan_text)))
    out_rows.append({
        "assembly_accession": acc_of[accession],
        "biosample_accession": accession,
        "fetch_status": "official_record_fetched",
        "package": package,
        "description_title": title,
        "description_comment": comment,
        "attributes_json": json.dumps(attrs, ensure_ascii=False, sort_keys=True),
        "bioproject_links": ";".join(sorted(set(links))),
        "keyword_hits": ";".join(hits),
        "raw_file": f"biosample_raw_fetch_v7/{accession}.xml",
        "raw_sha256": hashlib.sha256(block.encode()).hexdigest(),
        "access_date": "2026-10-10",
    })

assert seen == set(acc_of), f"missing: {sorted(set(acc_of) - seen)}"
with open(OUT / "biosample_gap_fetch_v7.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, delimiter="\t", fieldnames=list(out_rows[0].keys()))
    w.writeheader()
    w.writerows(sorted(out_rows, key=lambda r: r["assembly_accession"]))

with open(OUT / "sha256_manifest.txt", "w") as fh:
    for p in sorted(OUT.iterdir()):
        if p.name == "sha256_manifest.txt":
            continue
        fh.write(f"{hashlib.sha256(p.read_bytes()).hexdigest()}  {p.name}\n")
print(f"fetched {len(out_rows)}/41 biosample records -> {OUT}")
for r in out_rows:
    if r["keyword_hits"]:
        print("KEYWORD:", r["assembly_accession"], r["biosample_accession"], r["keyword_hits"],
              "|", (r["description_comment"] or r["description_title"])[:100])
