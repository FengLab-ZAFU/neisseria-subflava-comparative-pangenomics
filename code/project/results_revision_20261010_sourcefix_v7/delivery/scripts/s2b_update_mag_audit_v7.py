#!/usr/bin/env python3
"""sourcefix-v7 supplement S2b: updated MAG source audit table.

The 2026-10-06 MAG audit (results_revision_20261006_fix/01_mag_audit/
mag_source_audit.tsv) missed two confirmed metagenome-assembled genomes
because it never scanned the BioSample description.comment channel:

  GCA_963264865.1 SAMEA114260487 NZ_220527_assembly_030 metaFlye
  GCA_963281055.1 SAMEA114260481 NZ_220527_assembly_024 metaFlye

Evidence (frozen records, verified 2026-10-10 against the server frozen genus
JSONL results_revision_20261001/01_retrieval_taxonomy/raw/genus_482_current_MERGED.jsonl
and the reviewer's attached raw records v6再次复核_两条MAG冻结原始记录_20261010.json):
  biosample description.comment = "derived from the metagenomic samples
  ERS16237632 - ERS16237645" (project PRJEB61647); assembly_method=metaFlye
  corroborates. Complete Genome / single contig / Generic package do NOT
  override an explicit metagenomic source statement.

Output: results_revision_20261010_sourcefix_v7/02_source_evidence/mag_source_audit_v7.tsv
= old 920 rows + these 2 rows with decision=EXCLUDE_confirmed_MAG.
h5_cohort_freeze_v7.py reads THIS file (not the old hardcoded path).
"""
import os
import csv
import json
from pathlib import Path

ROOT = Path(os.environ["NEISSERIA_PROJECT_ROOT"])
OLD = ROOT / "results_revision_20261006_fix/01_mag_audit/mag_source_audit.tsv"
V7 = ROOT / "results_revision_20261010_sourcefix_v7"
OUT = V7 / "02_source_evidence/mag_source_audit_v7.tsv"

NEW = {
    "GCA_963264865.1": {"biosample": "SAMEA114260487", "assembly_name": "NZ_220527_assembly_030",
                        "level": "Complete Genome", "organism": "Neisseria flavescens"},
    "GCA_963281055.1": {"biosample": "SAMEA114260481", "assembly_name": "NZ_220527_assembly_024",
                        "level": "Complete Genome", "organism": "Neisseria flavescens"},
    # found by the v7 supplement sec.3 per-record re-check of mixed project PRJEB61647;
    # official BioSample fetched 2026-10-10 carries the same metagenomic comment.
    # ANI non-member (84.59) -> in_screen1303 True but never in any cohort pool.
    "GCA_963264735.1": {"biosample": "SAMEA114260461", "assembly_name": "NZ_220527_assembly_004",
                        "level": "Contig", "organism": "Neisseria sicca"},
}
COMMENT = "derived from the metagenomic samples ERS16237632 - ERS16237645"

# verify against server frozen JSONL (assembly side) before writing
frozen = {}
for line in open(ROOT / "results_revision_20261001/01_retrieval_taxonomy/raw/genus_482_current_MERGED.jsonl"):
    d = json.loads(line)
    if d.get("accession") in NEW and d["accession"].startswith("GCA_"):
        frozen[d["accession"]] = d
for acc, m in NEW.items():
    d = frozen[acc]
    ai = d["assembly_info"]
    assert ai.get("assembly_method") == "metaFlye", (acc, ai.get("assembly_method"))
    assert ai.get("assembly_name") == m["assembly_name"], (acc, ai.get("assembly_name"))
    assert ai.get("bioproject_accession") == "PRJEB61647", (acc, ai.get("bioproject_accession"))
    assert (ai.get("biosample") or {}).get("accession") == m["biosample"]

man = {r["assembly_accession"]: r for r in csv.DictReader(
    open(ROOT / "results_revision_20261001/02_independence/retrieval_manifest.tsv"), delimiter="\t")}

rows = list(csv.DictReader(open(OLD), delimiter="\t"))
cols = list(rows[0].keys())
assert not any(r["accession"] in NEW for r in rows), "already present in old audit?"

for acc, m in sorted(NEW.items()):
    mm = man[acc]
    reason = ("BioSample description.comment declares metagenomic source samples "
              "ERS16237632-ERS16237645 (PRJEB61647); metaFlye assembly method "
              "corroborates; missed by 2026-10-06 audit (description channel not "
              "scanned); added in v7 supplement review 2026-10-10")
    if acc == "GCA_963264735.1":
        reason = ("Official BioSample record (fetched 2026-10-10, biosample_raw_fetch_v7) "
                  "declares metagenomic source samples ERS16237632-ERS16237645 (PRJEB61647); "
                  "metaFlye corroborates; ANI non-member (84.59) so never entered any cohort "
                  "pool; added in v7 supplement review for audit completeness")
    rows.append({
        "accession": acc, "current_accession": acc, "biosample": m["biosample"],
        "organism_name": m["organism"], "assembly_name": m["assembly_name"],
        "assembly_level": m["level"], "bioproject": "PRJEB61647",
        "submitter": "WASHINGTON UNIVERSITY CENTER FOR GENOME SCIENCES",
        "biosample_package": "Generic.1.0", "isolation_source": "", "host": "",
        "evidence_channels": "E4_v7_description_comment",
        "evidence_detail": (f'E4_biosample_description_comment: "{COMMENT}" || '
                            "assembly_method=metaFlye (frozen assembly record)"),
        "decision": "EXCLUDE_confirmed_MAG",
        "decision_reason": reason,
        "in_screen1303": "True",
        "in_ani358": "True" if mm.get("fate", "") not in ("", "non_member") else "False",
        "in_revision_primary287": mm.get("in_revision_primary", ""),
        "ani_group": mm.get("group_for_testing", ""),
    })

rows.sort(key=lambda x: (x["decision"], x["accession"]))
with open(OUT, "w", newline="") as fh:
    w = csv.DictWriter(fh, delimiter="\t", fieldnames=cols)
    w.writeheader()
    w.writerows(rows)
print(f"written {OUT} rows={len(rows)} (old {len(rows)-3} + 3 new MAG)")
