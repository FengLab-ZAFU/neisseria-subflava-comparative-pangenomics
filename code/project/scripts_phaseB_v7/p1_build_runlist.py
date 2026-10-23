#!/usr/bin/env python3
"""P1 step 3a: build the Bakta run list (19 pending + cache downgrades) and the
annotation_sources.tsv used by the contig binding check and final manifest.

Outputs (01_annotation/):
  bakta_runlist.tsv       accession, fasta_path, reason
  annotation_sources.tsv  accession, annotation_source, annotation_dir
"""
import os
import csv
from pathlib import Path

ROOT = Path(os.environ["NEISSERIA_PROJECT_ROOT"])
V7 = ROOT / "results_revision_20261010_sourcefix_v7"
OUT = ROOT / "results_revision_20261010_phaseB_v7/01_annotation"
CACHE = ROOT / "results_revision_20261001/03_annotation/bakta"
NEW = OUT / "bakta_new"
PENDING = Path(os.environ["ANNOTATION_PENDING_TSV"])


def read_tsv(p):
    with open(p, newline="", encoding="utf-8-sig") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def main():
    binding = {r["accession"]: r for r in read_tsv(OUT / "fasta_binding_212.tsv")}
    cache = {r["accession"]: r for r in read_tsv(OUT / "annotation_cache_verification.tsv")}
    pending = read_tsv(PENDING)
    assert len(pending) == 19, f"pending list has {len(pending)} rows"

    runlist, sources = [], []
    for r in pending:
        acc = r["accession"]
        b = binding.get(acc)
        if b is None or b["match"] != "PASS":
            runlist.append({"accession": acc, "fasta_path": "",
                            "reason": "BLOCKED_fasta_binding_failed"})
        else:
            runlist.append({"accession": acc, "fasta_path": b["path"],
                            "reason": "pending19_" + r["existing_binding_status"]})
        sources.append({"accession": acc, "annotation_source": "new_bakta_1.12.0_db6.0",
                        "annotation_dir": str(NEW / acc)})
    for acc, c in sorted(cache.items()):
        if c["verdict"] == "REUSE_OK":
            sources.append({"accession": acc, "annotation_source": "reused_bakta_1.12.0_db6.0",
                            "annotation_dir": str(CACHE / acc)})
        else:
            b = binding.get(acc)
            if b is not None and b["match"] == "PASS":
                runlist.append({"accession": acc, "fasta_path": b["path"],
                                "reason": "cache_downgrade_" + c["reason"]})
                sources.append({"accession": acc, "annotation_source": "new_bakta_1.12.0_db6.0",
                                "annotation_dir": str(NEW / acc)})
            else:
                runlist.append({"accession": acc, "fasta_path": "",
                                "reason": "BLOCKED_cache_downgrade_and_fasta_binding_failed"})
                sources.append({"accession": acc, "annotation_source": "BLOCKED",
                                "annotation_dir": ""})

    assert len(sources) == 212, f"sources has {len(sources)} rows"
    with open(OUT / "bakta_runlist.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["accession", "fasta_path", "reason"], delimiter="\t")
        w.writeheader(); w.writerows(runlist)
    with open(OUT / "annotation_sources.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["accession", "annotation_source", "annotation_dir"], delimiter="\t")
        w.writeheader(); w.writerows(sorted(sources, key=lambda r: r["accession"]))
    runnable = [r for r in runlist if r["fasta_path"]]
    blocked = [r for r in runlist if not r["fasta_path"]]
    print(f"runlist: {len(runlist)} total, {len(runnable)} runnable, {len(blocked)} blocked")
    print(f"sources: {sum(1 for s in sources if s['annotation_source'].startswith('reused'))} reused, "
          f"{sum(1 for s in sources if s['annotation_source'].startswith('new'))} new, "
          f"{sum(1 for s in sources if s['annotation_source']=='BLOCKED')} blocked")
    for r in blocked:
        print("BLOCKED", r["accession"], r["reason"])


if __name__ == "__main__":
    main()
