#!/usr/bin/env python3
"""P1 step 2: verify the historical Bakta annotation cache for the 193 accessions
marked bound_ok_sequence_identical in v7 05_diff_ledgers/annotation_input_binding_v7.tsv.

Checks per accession:
  (a) bakta output dir exists with non-empty <acc>.gff3/.gbff/.faa/.ffn
  (b) per-genome log records bakta version=1.12.0 and the recorded database path
  (c) sequence identity: the FASTA actually annotated (bakta outdir <acc>.fna)
      has the same multiset of contig sequences as the frozen v7 FASTA
      (re-hashed/verified in step 1, fasta_binding_212.tsv)

Output: results_revision_20261010_phaseB_v7/01_annotation/annotation_cache_verification.tsv
Verdict REUSE_OK or DOWNGRADE:<reason>; downgraded accessions must be re-annotated.
"""
import os
import csv, re, sys
from pathlib import Path

ROOT = Path(os.environ["NEISSERIA_PROJECT_ROOT"])
V7 = ROOT / "results_revision_20261010_sourcefix_v7"
OUT = ROOT / "results_revision_20261010_phaseB_v7/01_annotation"
BAKTA_CACHE = ROOT / "results_revision_20261001/03_annotation/bakta"
DB_EXPECTED = os.environ["BAKTA_RECORDED_DB_PATH"]


def read_tsv(p):
    with open(p, newline="", encoding="utf-8-sig") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def seq_multiset(p):
    seqs = []
    cur = []
    with open(p) as fh:
        for line in fh:
            if line.startswith(">"):
                if cur:
                    seqs.append("".join(cur).upper())
                cur = []
            else:
                cur.append(line.strip())
    if cur:
        seqs.append("".join(cur).upper())
    return sorted(seqs)


def main():
    binding = {r["accession"]: r for r in read_tsv(OUT / "fasta_binding_212.tsv")}
    rows = []
    for r in read_tsv(V7 / "05_diff_ledgers/annotation_input_binding_v7.tsv"):
        if r["bakta_annotated"] != "True":
            continue
        acc = r["accession"]
        d = BAKTA_CACHE / acc
        rec = {"accession": acc, "bakta_dir": str(d), "dir_exists": "",
               "gff3": "", "gbff": "", "faa": "", "ffn": "",
               "log_exists": "", "log_version_1_12_0": "", "log_db_path_ok": "",
               "seq_identity_with_frozen_fasta": "",
               "verdict": "", "reason": ""}
        reasons = []
        rec["dir_exists"] = str(d.is_dir())
        if not d.is_dir():
            reasons.append("bakta_dir_missing")
        files = {}
        for ext in ("gff3", "gbff", "faa", "ffn"):
            f = d / f"{acc}.{ext}"
            ok = f.is_file() and f.stat().st_size > 0
            rec[ext] = "ok" if ok else "missing_or_empty"
            if not ok:
                reasons.append(f"{ext}_missing_or_empty")
            files[ext] = str(f) if ok else ""
        log = d / f"{acc}.log"
        if log.is_file():
            rec["log_exists"] = "yes"
            txt = log.read_text(errors="replace")
            rec["log_version_1_12_0"] = "yes" if re.search(r"version=1\.12\.0\b", txt) else "no"
            rec["log_db_path_ok"] = "yes" if f"database: type=parameter, path={DB_EXPECTED}" in txt else "no"
            if rec["log_version_1_12_0"] != "yes":
                reasons.append("log_bakta_version_not_1.12.0")
            if rec["log_db_path_ok"] != "yes":
                reasons.append("log_db_path_mismatch")
        else:
            rec["log_exists"] = "no"
            rec["log_version_1_12_0"] = "no"
            rec["log_db_path_ok"] = "no"
            reasons.append("log_missing")
        # sequence identity between annotated input and frozen fasta
        b = binding.get(acc)
        ann_fna = d / f"{acc}.fna"
        if b is None or b["match"] != "PASS":
            rec["seq_identity_with_frozen_fasta"] = "frozen_fasta_not_bound"
            reasons.append("frozen_fasta_not_bound")
        elif not ann_fna.is_file():
            rec["seq_identity_with_frozen_fasta"] = "annotated_input_fna_missing"
            reasons.append("annotated_input_fna_missing")
        else:
            same = seq_multiset(ann_fna) == seq_multiset(b["path"])
            rec["seq_identity_with_frozen_fasta"] = "identical" if same else "DIFFERENT"
            if not same:
                reasons.append("annotated_input_sequence_differs_from_frozen")
        rec["verdict"] = "REUSE_OK" if not reasons else "DOWNGRADE"
        rec["reason"] = ";".join(reasons)
        rows.append(rec)

    fields = list(rows[0].keys())
    with open(OUT / "annotation_cache_verification.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader(); w.writerows(rows)
    ok = sum(1 for r in rows if r["verdict"] == "REUSE_OK")
    print(f"cache verification: {ok}/{len(rows)} REUSE_OK, {len(rows)-ok} DOWNGRADE")
    for r in rows:
        if r["verdict"] != "REUSE_OK":
            print("DOWNGRADE", r["accession"], r["reason"])


if __name__ == "__main__":
    main()
