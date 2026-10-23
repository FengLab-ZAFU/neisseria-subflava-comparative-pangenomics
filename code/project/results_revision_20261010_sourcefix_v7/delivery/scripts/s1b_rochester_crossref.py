#!/usr/bin/env python3
"""fix3-supp S1b (v2): Rochester cross-project participant/isolate evidence.

Evidence chain per assembly record (P-coded, PRJNA514245-umbrella or PRJNA1273950):
  1. strain code PxxxxSxxx (BioSample strain field, normalized)
  2. SRA library with LIBRARY_NAME == code in the study SRA project
     PRJNA1208441 / SRP556955 (one run per code; PRJNA1273950 itself has ZERO
     SRA runs -> its assemblies have no own sequencing library)
  3. Regan 2026 Table S1 row with Sample == code -> authoritative Participant
     and Year.Collected; Table S1 accession column == the SRA run (verified
     166/166 consistent, zero mismatch)
  4. for the 55 umbrella records: biosample->sra elink run (raw JSON cached)

Outputs:
  01_rochester_evidence/participant_source_evidence.tsv
  01_rochester_evidence/s1_findings_summary.md
"""
import os
import csv
import json
import re
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(os.environ["NEISSERIA_PROJECT_ROOT"])
# output dir: default = v6 ns-fix supplement; pass argv[1] to override.
# 2026-10-10 ns-fix: one real study = one namespace; Table S1 coverage and
# participant evidence level are separate columns, never separate namespaces.
import sys
OUT = Path(sys.argv[1]) if len(sys.argv) > 1 else \
    ROOT / "results_revision_20261010_nsfix_v6/01_rochester_evidence"
OUT.mkdir(parents=True, exist_ok=True)
CACHE_V5 = ROOT / "results_revision_20261009_fix3_phaseA_supplement/01_rochester_evidence"
S1CSV = ROOT / "fix3_supplement_20261009/evidence/regan2026_suppl/spectrum.00725-26-s0002.csv"

def norm_code(s):
    s = (s or "").strip().upper().replace(".", "")
    m = re.fullmatch(r"(P\d{4})S(\d{3})", s)
    return f"{m.group(1)}S{m.group(2)}" if m else ""

# --- Table S1 ---
allr = list(csv.reader(open(S1CSV, newline="")))
hdr = allr[1]
s1_by_code = {}
for ln, r in enumerate(allr[2:], start=3):
    if r and r[0].strip().startswith("P"):
        d = dict(zip(hdr, r)); d["_line"] = ln
        s1_by_code[d["Sample"].strip()] = d
print("Table S1 isolate rows:", len(s1_by_code))

# --- SRA library map (s1c) ---
lib = {}
for r in csv.DictReader(open(CACHE_V5 / "sra_run_library_names.tsv"), delimiter="\t"):
    lib[r["library_name"]] = r
print("SRA libraries:", len(lib))

# --- elink map (s1a) ---
elink = {r["biosample"]: r for r in csv.DictReader(open(CACHE_V5 / "biosample_sra_map.tsv"), delimiter="\t")}

# --- true source project per biosample (frozen JSONL extract) ---
src_proj = {}
for line in open(ROOT / "results_revision_20261009_fix3_phaseA_supplement/01_rochester_evidence/biosample_source_projects_rochester.jsonl"):
    d = json.loads(line)
    if d.get("biosample"):
        src_proj[d["biosample"]] = d.get("biosample_bioprojects") or ""

# --- provenance P records ---
prov = list(csv.DictReader(open(ROOT / "results_revision_20261009_fix3/02_provenance_v4/provenance_corrected.tsv"), delimiter="\t"))
prec = []
for r in prov:
    c = norm_code(r.get("strain", ""))
    if c and r["bioproject_accession"] in ("PRJNA514245", "PRJNA1273950"):
        r["_code"] = c
        prec.append(r)
print("P-coded provenance records:", len(prec))

records = []
for p in sorted(prec, key=lambda r: r["assembly_accession"]):
    code = p["_code"]
    bs = p["biosample_accession"]
    el_run = (elink.get(bs, {}).get("sra_runs") or "")
    l = lib.get(code)
    s1 = s1_by_code.get(code)
    rec = {
        "assembly_accession": p["assembly_accession"],
        "biosample": bs,
        "isolate_code": code,
        "submitted_bioproject": p["bioproject_accession"],
        "true_source_bioproject": src_proj.get(bs, ""),
        "original_strain_field": p.get("strain", ""),
        "biosample_collection_date": p.get("collection_date", ""),
        "biosample_sra_elink_run": el_run,
        "study_sra_project": "PRJNA1208441 (SRP556955)",
        "study_sra_run": l["run"] if l else "",
        "study_sra_run_biosample": l["sra_biosample"] if l else "",
        "table_s1_match": "yes" if s1 else "no",
        "table_s1_csv_line": s1["_line"] if s1 else "",
        "table_s1_sample": s1["Sample"] if s1 else "",
        "table_s1_participant": s1["Participant"] if s1 else "",
        "table_s1_isolate": s1["Isolate"] if s1 else "",
        "table_s1_year_collected": s1["Year.Collected"] if s1 else "",
        "table_s1_accession_srr": s1["accession"] if s1 else "",
    }
    # evidence verdicts
    srr_consistent = s1 and l and s1["accession"].strip() == l["run"]
    rec["srr_library_tables1_consistent"] = str(bool(srr_consistent))
    # 2026-10-10 ns-fix: publication coverage and participant evidence level are
    # separate fields; the real-study namespace is unified for all records that
    # map to PRJNA1208441/SRP556955.
    if s1:
        rec["namespace"] = "Regan2026_Rochester_study"
        rec["publication_coverage"] = "table_s1"
        rec["participant_evidence_level"] = "high" if l else "medium"
        rec["verified_participant"] = s1["Participant"]
        rec["verified_year"] = s1["Year.Collected"]
        rec["evidence_summary"] = (
            f"Regan2026 Table S1 line {s1['_line']}: Sample={s1['Sample']} Participant={s1['Participant']} "
            f"Year.Collected={s1['Year.Collected']}; Table S1 SRR={s1['accession']} == study SRA library "
            f"{l['run'] if l else 'MISSING'} (LIBRARY_NAME={code})")
    elif l:
        rec["namespace"] = "Regan2026_Rochester_study"
        rec["publication_coverage"] = "sra_only"
        rec["participant_evidence_level"] = "medium"
        rec["verified_participant"] = code.split("S")[0]  # participant prefix from code
        rec["verified_year"] = ""
        rec["evidence_summary"] = (
            f"Not published in Regan2026 Table S1, but SRA library {l['run']} with LIBRARY_NAME={code} "
            f"exists in the same study SRA project PRJNA1208441/SRP556955; participant prefix "
            f"{rec['verified_participant']} inferred from isolate code (code prefix == Participant in "
            f"166/166 Table S1 rows), confidence medium; year unknown")
    else:
        rec["namespace"] = f"{p['bioproject_accession']}_unverified"
        rec["publication_coverage"] = "none"
        rec["participant_evidence_level"] = "low"
        rec["verified_participant"] = "unknown"
        rec["verified_year"] = ""
        rec["evidence_summary"] = "No Table S1 row and no SRA library with this code; provenance unverified"
    rec["confidence"] = rec["participant_evidence_level"]
    records.append(rec)

# --- isolate-level cross-project relations ---
by_code = defaultdict(list)
for r in records:
    by_code[r["isolate_code"]].append(r)
for r in records:
    grp = by_code[r["isolate_code"]]
    projs = sorted({g["true_source_bioproject"] or g["submitted_bioproject"] for g in grp})
    if len(grp) > 1:
        same_lib = len({g["study_sra_run"] for g in grp}) == 1 and grp[0]["study_sra_run"]
        rel = "same_isolate_duplicate_submission" if same_lib else "same_isolate_library_unresolved"
        if len(projs) > 1:
            rel += ";cross_project(" + ",".join(projs) + ")"
        else:
            rel += ";within_project"
        r["cross_project_relation"] = rel
        # layered evidence wording (2026-10-10 ns-fix): duplicate submission is
        # supported by study source + unique isolate code + Table S1/SRA mapping;
        # "same reads library" is direct only where a biosample->SRA link exists,
        # otherwise inferred (PRJNA1273950 deposited no reads of its own).
        direct = bool(r["biosample_sra_elink_run"])
        r["same_reads_library_evidence"] = ("direct_biosample_sra_link" if direct else
                                            "inferred_via_unique_isolate_code_no_reads_in_PRJNA1273950")
        r["disposition"] = ("source_identity_duplicate: merge assemblies of this isolate code into one "
                            "source-identity group; QC-ranked representative; NOT a same-participant "
                            "resampling clone")
    else:
        r["cross_project_relation"] = "unique_isolate"
        r["same_reads_library_evidence"] = "not_applicable"
        r["disposition"] = ("participant-level grouping only within namespace " + r["namespace"] +
                            "; cross-isolate merging requires same verified participant AND SNP rule")

cols = ["assembly_accession", "biosample", "isolate_code", "submitted_bioproject", "true_source_bioproject",
        "original_strain_field", "biosample_collection_date", "biosample_sra_elink_run",
        "study_sra_project", "study_sra_run", "study_sra_run_biosample", "srr_library_tables1_consistent",
        "table_s1_match", "table_s1_csv_line", "table_s1_sample", "table_s1_participant", "table_s1_isolate",
        "table_s1_year_collected", "table_s1_accession_srr", "namespace", "publication_coverage",
        "participant_evidence_level", "verified_participant",
        "verified_year", "cross_project_relation", "same_reads_library_evidence", "evidence_summary",
        "disposition", "confidence"]
with open(OUT / "participant_source_evidence.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, delimiter="\t", fieldnames=cols)
    w.writeheader()
    w.writerows(records)

# --- findings summary ---
c_rel = Counter(r["cross_project_relation"].split(";")[0] for r in records)
c_ns = Counter(r["namespace"] for r in records)
dup_codes = {k: v for k, v in by_code.items() if len(v) > 1}
xp = {k: v for k, v in dup_codes.items()
      if len({g["true_source_bioproject"] or g["submitted_bioproject"] for g in v}) > 1}
yr = defaultdict(set)
for r in records:
    if r["verified_participant"] != "unknown" and r["verified_year"]:
        yr[r["verified_year"]].add(r["verified_participant"])

with open(OUT / "s1_findings_summary.md", "w") as fh:
    fh.write(f"""# S1 findings: Rochester cross-project participant evidence (2026-10-10, ns-fix)

Evidence base: Regan et al. 2026, Microbiol Spectr, doi:10.1128/spectrum.00725-26
(PMC13340140, open access). Table S1 = spectrum.00725-26-s0002.csv (evidence copy with
checksums in fix3_supplement_20261009/evidence/). SRA evidence: eutils esearch/elink/
esummary raw JSON in results_revision_20261009_fix3_phaseA_supplement/01_rochester_evidence/eutils_cache/
(accessed 2026-10-09; zero-run status of PRJNA1273950 re-confirmed online 2026-10-10).

## Study structure (empirical)
- One study, one SRA project: PRJNA1208441 / SRP556955 "Commensal Neisseria collection
  from Rochester, NY" (RIT). 214 runs, every run's LIBRARY_NAME is a unique PxxxxSxxx
  isolate code (1 run per code).
- Assembly BioProjects are organizational views, not sequencing projects:
  PRJNA514245 (umbrella) -> BioSample-level source project is PRJNA1208441 (2024 batch,
  SKESA assemblies); PRJNA1273950 (2025, SPAdes assemblies) has ZERO SRA runs of its own.
- Table S1 covers 166 isolates, 31 participants; 2024 participants P0001-P0017 segment,
  2025 participants P0018-P0037 segment, ZERO overlap in participants or sample codes
  across years. No participant was sampled in both years.
- Table S1 accession column == SRA run with LIBRARY_NAME == Sample code: 166/166
  consistent, zero mismatches.
- Table S1 Sample code prefix == Participant column: 166/166 consistent -> the code
  prefix is an empirically validated participant key within this study (basis for the
  medium-confidence participant inference of sra_only records).

## Screen-scope records ({len(records)} P-coded assemblies in the ANI screen)
- namespaces: {dict(c_ns)}
- publication coverage: {dict(Counter(r["publication_coverage"] for r in records))}
- participant evidence levels: {dict(Counter(r["participant_evidence_level"] for r in records))}
- relations: {dict(c_rel)}
- distinct isolate codes: {len(by_code)}; codes with >1 assembly: {len(dup_codes)},
  ALL of them cross-project (PRJNA1208441-sourced + PRJNA1273950): {len(xp)}
- of the duplicate codes: {sum(1 for k in dup_codes if k in s1_by_code)} in Table S1,
  {sum(1 for k in dup_codes if k not in s1_by_code and k in lib)} only in study SRA,
  {sum(1 for k in dup_codes if k not in s1_by_code and k not in lib)} unverified
- verified participants by year: 2024 -> {len(yr.get('2024', set()))}, 2025 -> {len(yr.get('2025', set()))},
  overlap -> {sorted(yr.get('2024', set()) & yr.get('2025', set()))}

## Answer to the review question (2024 vs 2025 same P/S codes) - layered evidence wording
Same PxxxxSxxx code in the two assembly BioProjects = THE SAME ISOLATE submitted twice
(PRJNA1208441 2024 SKESA batch, visible under umbrella PRJNA514245; and the 2025
whole-study SPAdes re-assembly PRJNA1273950). The duplicate-submission conclusion rests
on: (i) the same study source (all reads under SRP556955), (ii) the unique isolate code
(one Table S1 row / one SRA LIBRARY_NAME per code), and (iii) the Table S1/SRA mapping.
NOTE: the stronger claim that both assemblies use the SAME reads library is DIRECT only
where a BioSample->SRA elink exists (55 umbrella records); for PRJNA1273950 assemblies
it is INFERRED (that project deposited no reads, and the isolate code maps to exactly
one library in the study SRA project). This distinction does not change the
source-identity duplicate handling.
=> same-isolate DUPLICATE SUBMISSION -> handled by the source-identity duplicate rule
   (QC-ranked representative), NOT by the same-participant SNP-clone rule.
It is NOT same-participant resampling (sample codes are isolate-level; the two years'
participant segments do not overlap) and NOT code reuse (one Table S1 row / one library
per code).

## Namespace policy (2026-10-10 fix)
One real study = one namespace: all records mapping to PRJNA1208441/SRP556955 carry
namespace Regan2026_Rochester_study. Whether an isolate appears in Table S1 is recorded
in publication_coverage (table_s1 / sra_only / none), and participant evidence strength
in participant_evidence_level (high = Table S1 row; medium = SRA library code prefix,
validated by the 166/166 prefix==Participant concordance; low = unverified). Neither
field splits the study namespace. The 7 records with neither Table S1 row nor study SRA
library stay unknown with independent model ids and are never merged.
""")
print("relations:", dict(c_rel))
print("namespaces:", dict(c_ns))
print("wrote participant_source_evidence.tsv and s1_findings_summary.md")
