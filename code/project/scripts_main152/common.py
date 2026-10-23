#!/usr/bin/env python3
"""Shared paths, cohort metadata and input checks for the main152 analyses."""
from pathlib import Path
import csv
import hashlib
import os

ROOT = str(Path(os.environ["NEISSERIA_PROJECT_ROOT"]))
SOURCE_ROOT = os.environ.get("NEISSERIA_SOURCE_ROOT", ROOT)

RESULTS = os.path.join(ROOT, "results_revision_20261015_main152_species")
PACKET = os.path.join(RESULTS, "00_freeze/input_packet_snapshot")
SCRIPTS = os.path.join(ROOT, "scripts_main152")
MASTER_SEED = 20261015
RUN_ID = "main152_20261015"
PLAN_ID = "main152_species_midpoint_v3_audited"

ROSTER_TXT = os.path.join(PACKET, "00_input/cohort/main152_accessions.txt")
GROUPS = ["subflava", "perflava", "flavescens"]
GROUP_SIZES = {"subflava": 31, "perflava": 114, "flavescens": 7}


def read_roster():
    with open(ROSTER_TXT, "rb") as fh:
        accs = [ln.strip() for ln in fh.read().decode("utf-8").splitlines() if ln.strip()]
    assert len(accs) == 152 and len(set(accs)) == 152, f"roster problem: {len(accs)}"
    return accs


def read_tsv(path):
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(path, fieldnames, rows):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t",
                           extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow(r)


def sha256_file(path, _buf=1024 * 1024):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while True:
            b = fh.read(_buf)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def sha256_bytes(b):
    return hashlib.sha256(b).hexdigest()


def resolve_path(recorded):
    """Resolve an input path in the prepared or external source-data root."""
    if not recorded:
        return None, "empty"
    if os.path.exists(recorded):
        return recorded, "as_recorded"
    relative = recorded.replace("\\", "/")
    marker = "/Neisseria_subflava/"
    if marker in relative:
        relative = relative.split(marker, 1)[1]
    elif recorded.startswith(ROOT):
        relative = os.path.relpath(recorded, ROOT)
    elif os.path.isabs(recorded):
        return None, "MISSING"
    for base in [ROOT, os.environ.get("NEISSERIA_SOURCE_ROOT", ROOT)]:
        candidate = Path(base) / relative
        if candidate.exists():
            return str(candidate), "resolved_relative"
    return None, "MISSING"


def load_metadata():
    return read_tsv(os.path.join(PACKET, "00_input/cohort/main152_metadata.tsv"))


def load_group_map():
    """accession -> species_group (group_for_testing in frozen metadata)."""
    return {r["strain"]: r["group_for_testing"] for r in load_metadata()}
