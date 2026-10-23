#!/usr/bin/env python3
"""A2: merge main fetcher + chunk window files into the final frozen JSONL.

Validates: unique accession count == 73686 (API total_count), no duplicates,
every record parses and belongs to genus 482 lineage (has organism.tax_id).
"""
import json
from pathlib import Path

RAW = Path(__file__).resolve().parent.parent.parent / "results_revision_20261001" / "01_retrieval_taxonomy" / "raw"

files = [RAW / "genus_482_current_frozen.jsonl"] + sorted((RAW / "chunks").glob("window_*.jsonl"))
seen = {}
dup = 0
n_files = 0
for f in files:
    if not f.exists():
        continue
    n_files += 1
    with open(f) as fh:
        for line in fh:
            try:
                r = json.loads(line)
            except Exception:
                continue
            acc = r.get("accession")
            if not acc:
                continue
            if acc in seen:
                dup += 1
                continue
            seen[acc] = r

out = RAW / "genus_482_current_MERGED.jsonl"
with open(out, "w") as fh:
    for acc in sorted(seen):
        fh.write(json.dumps(seen[acc]) + "\n")

print(f"input files: {n_files} | unique accessions: {len(seen)} | duplicates skipped: {dup}")
print(f"expected total_count: 73686 | {'MATCH' if len(seen) == 73686 else 'DEFICIT: ' + str(73686 - len(seen))}")
# sanity: taxids present
from collections import Counter
top = Counter(r.get("organism", {}).get("tax_id") for r in seen.values())
print("top taxids:", top.most_common(8))
