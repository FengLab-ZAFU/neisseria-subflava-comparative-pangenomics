#!/usr/bin/env python3
"""A2: tier the frozen Neisseria genus metadata and produce the screening scope.

Inputs:
  raw/genus_482_current_MERGED.jsonl   - frozen genus metadata (2026-10-01)
  raw/gono_meno_descendants_485.jsonl      - N. gonorrhoeae subtree taxids
  raw/gono_meno_descendants_487.jsonl      - N. meningitidis subtree taxids

Outputs (results_revision_20261001/01_retrieval_taxonomy/):
  tier_summary.tsv            - reproducible exclusion counts per source rule
  screening_scope_accessions.tsv - in-scope GenBank accessions to download+screen
  out_of_frame_gono_meno.tsv  - excluded gonorrhoeae/meningitidis subtree records
  full_ledger.tsv             - every frozen record with exclusion category
"""
import csv
import json
from collections import Counter
from pathlib import Path

RAW = Path(__file__).resolve().parent.parent.parent / "results_revision_20261001" / "01_retrieval_taxonomy" / "raw"
OUT = RAW.parent
META189 = Path(__file__).resolve().parent.parent.parent / "neisseria_subflava_associated_group_metadata.tsv"

# --- load gono/meno descendant taxids ---
excl_taxids = {485: None, 487: None}
for tid, fn in [(485, "gono_meno_descendants_485.jsonl"), (487, "gono_meno_descendants_487.jsonl")]:
    n = 0
    with open(RAW / fn) as fh:
        for line in fh:
            r = json.loads(line)
            tx = r["taxonomy"]
            tid_rec = tx["tax_id"]
            excl_taxids[tid_rec] = tid_rec
            parents = tx.get("parents", [])
            assert 482 in parents or tid_rec in (485, 487), f"unexpected lineage for {tid_rec}: {parents}"
            n += 1
    print(f"taxid {tid} subtree: {n} named taxa (+itself)")

records = []
bad = 0
with open(RAW / "genus_482_current_MERGED.jsonl") as fh:
    for line in fh:
        try:
            records.append(json.loads(line))
        except json.JSONDecodeError:
            bad += 1
print(f"frozen records: {len(records)} (unparseable lines: {bad})")

# distinct assembly_type values (for MAG detection audit)
type_counts = Counter(r.get("assembly_info", {}).get("assembly_type") for r in records)
print("assembly_type values:", dict(type_counts))

ledger = []
counts = Counter()
scope = []
for r in records:
    acc = r["accession"]
    cur = r.get("current_accession", acc)
    ai = r.get("assembly_info", {})
    org = r.get("organism", {})
    taxid = org.get("tax_id")
    status = ai.get("assembly_status")
    source = r.get("source_database", "")
    atypical = bool(ai.get("atypical"))
    asm_type = ai.get("assembly_type", "") or ""
    is_mag = "metagenome" in asm_type.lower()

    cat = "IN_SCOPE"
    if status != "current":
        cat = "EXCLUDED_suppressed"
    elif source != "SOURCE_DATABASE_GENBANK":
        cat = "EXCLUDED_refseq_only"
    elif atypical:
        cat = "EXCLUDED_atypical"
    elif is_mag:
        cat = "EXCLUDED_mag"
    elif taxid in excl_taxids:
        cat = f"OUT_OF_FRAME_gono_meno"
    else:
        scope.append(r)

    counts[cat] += 1
    ledger.append({
        "accession": acc, "current_accession": cur,
        "biosample": (ai.get("biosample") or {}).get("accession", ""),
        "organism_name": org.get("organism_name", ""),
        "tax_id": taxid or "",
        "assembly_level": ai.get("assembly_level", ""),
        "assembly_status": status, "source_database": source,
        "atypical": atypical, "assembly_type": asm_type,
        "assembly_name": ai.get("assembly_name", ""),
        "release_date": ai.get("release_date", ""),
        "total_sequence_length": (r.get("assembly_stats") or {}).get("total_sequence_length", ""),
        "exclusion_category": cat,
    })

# tier summary (reproducible counts)
with open(OUT / "tier_summary.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["step", "rule", "n_remaining", "n_excluded_at_step"])
    total = len(records)
    w.writerow(["0_frozen", "Neisseria taxid 482, --assembly-version current, no source/MAG/atypical filters", total, 0])
    rem = counts["IN_SCOPE"] + sum(v for k, v in counts.items() if k.startswith("OUT_OF_FRAME"))
    w.writerow(["1_status", "assembly_status == current (excl. suppressed)", rem, counts["EXCLUDED_suppressed"]])
    rem -= counts["EXCLUDED_refseq_only"]
    w.writerow(["2_source", "GenBank only (drop RefSeq-only)", rem, counts["EXCLUDED_refseq_only"]])
    rem -= counts["EXCLUDED_atypical"]
    w.writerow(["3_atypical", "exclude atypical", rem, counts["EXCLUDED_atypical"]])
    rem -= counts["EXCLUDED_mag"]
    w.writerow(["4_mag", "exclude MAG", rem, counts["EXCLUDED_mag"]])
    w.writerow(["5_species", "exclude gonorrhoeae/meningitidis subtree taxids (final screening scope)", counts["IN_SCOPE"], counts["OUT_OF_FRAME_gono_meno"]])

with open(OUT / "screening_scope_accessions.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["accession", "biosample", "organism_name", "tax_id", "assembly_level",
                "assembly_name", "release_date", "total_sequence_length", "n_contigs",
                "contig_n50", "gc_percent", "isolation_source", "host", "geo_loc_name", "collection_date"])
    for r in sorted(scope, key=lambda x: x["accession"]):
        ai = r.get("assembly_info", {})
        st = r.get("assembly_stats", {})
        bs = ai.get("biosample") or {}
        org = r.get("organism", {})
        w.writerow([r["accession"], bs.get("accession", ""), org.get("organism_name", ""), org.get("tax_id", ""),
                    ai.get("assembly_level", ""), ai.get("assembly_name", ""), ai.get("release_date", ""),
                    st.get("total_sequence_length", ""), st.get("number_of_contigs", ""),
                    st.get("contig_n50", ""), st.get("gc_percent", ""),
                    bs.get("isolation_source", ""), bs.get("host", ""), bs.get("geo_loc_name", ""), bs.get("collection_date", "")])

with open(OUT / "out_of_frame_gono_meno.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, delimiter="\t", fieldnames=list(ledger[0].keys()))
    w.writeheader()
    for row in ledger:
        if row["exclusion_category"] == "OUT_OF_FRAME_gono_meno":
            w.writerow(row)

with open(OUT / "full_ledger.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, delimiter="\t", fieldnames=list(ledger[0].keys()))
    w.writeheader()
    w.writerows(sorted(ledger, key=lambda x: x["accession"]))

# --- historical 189 coverage check ---
hist = {}
with open(META189, encoding="utf-8-sig") as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        hist[row["assembly_accession"]] = row["biosample_accession"]
frozen_acc = {r["accession"] for r in records}
frozen_cur = {r.get("current_accession", r["accession"]) for r in records}
missing = {a: b for a, b in hist.items() if a not in frozen_acc and a not in frozen_cur}
updated = {a: b for a, b in hist.items() if a not in frozen_acc and a in frozen_cur}
print(f"\nHistorical 189: in freeze={len(hist)-len(missing)-len(updated)}, version-updated={len(updated)}, missing-from-freeze={len(missing)}")
for a, b in list(missing.items())[:20]:
    print("  MISSING:", a, b)
for a, b in list(updated.items())[:20]:
    print("  UPDATED:", a, "->", b)

print(f"\nfinal screening scope: {counts['IN_SCOPE']} assemblies (out-of-frame gono/meno: {counts['OUT_OF_FRAME_gono_meno']})")
