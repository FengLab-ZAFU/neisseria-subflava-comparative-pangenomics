#!/usr/bin/env python3
"""M8 (main152): finalize run_info + sha256 ledger for 09_interval_overlap.

Usage: m8_finalize.py <outdir> <tables_dir>
"""
from pathlib import Path
import csv
import hashlib
import json
import os
import platform
import subprocess
import sys
from datetime import datetime, timezone

ROOT = str(Path(os.environ["NEISSERIA_PROJECT_ROOT"]))


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    outdir, tables = sys.argv[1], sys.argv[2]
    inputs = {
        "occurrence_map": f"{ROOT}/results_revision_20261015_main152_species/02_pangenome/pgap2_main152/pav_tables/occurrence_map.tsv",
        "partition_table": f"{ROOT}/results_revision_20261015_main152_species/02_pangenome/pgap2_main152/pav_tables/partition_table.tsv",
        "annotation_manifest": f"{ROOT}/results_revision_20261015_main152_species/01_annotation/annotation_manifest_main152.tsv",
        "main152_metadata": f"{ROOT}/results_revision_20261015_main152_species/00_freeze/input_packet_snapshot/00_input/cohort/main152_metadata.tsv",
        "per_tool_coverage": f"{ROOT}/results_revision_20261015_main152_species/00_freeze/input_packet_snapshot/00_input/reference_sources/per_tool_coverage_historical_main152.tsv",
    }
    info = {
        "module": "M8_interval_overlap",
        "run_id": "main152_m8_20261015",
        "hostname": platform.node(),
        "utc_finalized": datetime.now(timezone.utc).isoformat(),
        "python": sys.version.split()[0],
        "numpy": __import__("numpy").__version__,
        "master_seed": 20261015,
        "n_iter": 1000,
        "family_id": "interval_main_and_species",
        "coordinate_convention": "0-based half-open throughout; unified key (accession, frozen_contig_id)",
        "scripts": ["scripts_main152/m8_build_intervals.py",
                    "scripts_main152/m8_overlap_null.py",
                    "scripts_main152/m8_gi_cases.py",
                    "scripts_main152/m8_coordinate_audit.py",
                    "scripts_main152/m8_finalize.py"],
        "inputs_sha256": {k: sha256(v) for k, v in inputs.items()},
        "key_findings_and_decisions": [
            "occurrence_map.tsv coordinates are empirically 0-based half-open (PGAP2 total.involved_annot), NOT 1-based closed: verified against Bakta GFF for 328376/328376 rows with zero exceptions; the handoff note claiming 1-based closed was overridden by this evidence.",
            "geNomad provirus intervals were mapped to the parent contig via source_seq (1-based closed, verified by exact subsequence match); the historical 169-cohort run keyed proviruses by the nested 'seq|provirus_x_y' name so provirus-gene overlaps were invisible there -- new hit counts are therefore materially higher and correct.",
            "IntegronFinder: zero integrons in all 152 genomes (and zero in all 301 historical genomes); parser path kept and documented.",
            "MOB-suite: zero plasmid contig calls in all 152 genomes; whole-contig summary contains geNomad predictions only.",
            "GCA_901873995.1: one find_proviruses provirus (CABFLW010000022.1:83-10458) not classified as virus in the geNomad summary step; kept as local MGE from find_proviruses output (documented deviation from the old virus_summary-based parsing, +1 interval).",
            "Circularity: only the 9 contigs with explicit Bakta fna header [topology=circular] were treated as circular; all others linear/unknown.",
            "3 linear contigs have a degenerate single-legal-offset layout (IS spans the whole contig); their gene instances were excluded from the testable null set and observed recomputed on the same subset.",
            "MGE-free contigs with locatable genes were retained as fixed zero background.",
        ],
    }
    with open(os.path.join(outdir, "run_info.json")) as fh:
        prev = json.load(fh)
    prev.update(info)
    with open(os.path.join(outdir, "run_info.json"), "w") as fh:
        json.dump(prev, fh, indent=2, ensure_ascii=False)

    # sha256 ledger over all delivered files
    entries = []
    for dirpath, _, files in os.walk(outdir):
        for f in sorted(files):
            if f in ("checksums.sha256", "MANIFEST.tsv"):
                continue  # ledger files exclude themselves (self-reference)
            p = os.path.join(dirpath, f)
            entries.append((os.path.relpath(p, outdir), sha256(p), os.path.getsize(p)))
    with open(os.path.join(outdir, "checksums.sha256"), "w") as fh:
        for rel, h, sz in sorted(entries):
            fh.write(f"{h}  {rel}\n")
    with open(os.path.join(outdir, "MANIFEST.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["relpath", "bytes", "sha256"])
        for rel, h, sz in sorted(entries):
            w.writerow([rel, sz, h])
    print(f"ledger: {len(entries)} files")
    for rel, h, sz in sorted(entries):
        print(f"  {sz:>10}  {rel}")


if __name__ == "__main__":
    main()
