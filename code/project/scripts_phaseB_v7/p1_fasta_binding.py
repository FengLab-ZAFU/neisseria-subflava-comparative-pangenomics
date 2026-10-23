#!/usr/bin/env python3
"""P1 step 1: locate and hash-bind the 212 input FASTAs of the v7 annotation
preparation set (annotation_union_v7 193 U expanded_reference_sensitivity_v7 188).

Inputs (read-only):
  v7 freeze:  results_revision_20261010_sourcefix_v7/04_cohorts_v7/fasta_hashes_v7.tsv
              metadata_annotation_union_v7.tsv / metadata_expanded_reference_sensitivity_v7.tsv
  fallback search dirs: results_revision_20261001/01_retrieval_taxonomy/ (accession-named),
              Genomes/ (BioSample-named, only accepted when sha256 matches the freeze)
  prior covariates: results_revision_20261008_fix2/02_assembly_covariates/assembly_covariates_fasta.tsv

Outputs (results_revision_20261010_phaseB_v7/01_annotation/):
  fasta_binding_212.tsv, assembly_covariates_212.tsv, fasta_binding_failures.tsv
"""
import csv, hashlib, os, sys
from pathlib import Path

ROOT = Path(os.environ["NEISSERIA_PROJECT_ROOT"])
V7 = ROOT / "results_revision_20261010_sourcefix_v7"
OUT = ROOT / "results_revision_20261010_phaseB_v7/01_annotation"
OUT.mkdir(parents=True, exist_ok=True)

SEARCH_DIRS = [
    ROOT / "results_revision_20261001/01_retrieval_taxonomy",
    ROOT / "results_revision_20261001/01_retrieval_taxonomy/fasta",
]
COVAR_FILE = ROOT / "results_revision_20261008_fix2/02_assembly_covariates/assembly_covariates_fasta.tsv"

_INDEX = None

def accession_index():
    """Build once: accession.version -> [paths] for all fasta under the retrieval dir."""
    global _INDEX
    if _INDEX is None:
        _INDEX = {}
        base = SEARCH_DIRS[0]
        for p in base.rglob("*"):
            if p.is_file() and p.suffix in (".fna", ".fasta", ".fa"):
                stem = p.name.rsplit(".", 1)[0]
                _INDEX.setdefault(stem, []).append(p)
    return _INDEX


def read_tsv(p):
    with open(p, newline="", encoding="utf-8-sig") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def sha256_file(p, chunk=1 << 20):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        while True:
            b = fh.read(chunk)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def fasta_stats(p):
    """Return (n_contigs, total_len, n50, gc_percent) with a plain parser."""
    lens, gc, at = [], 0, 0
    cur = 0
    with open(p) as fh:
        for line in fh:
            if line.startswith(">"):
                if cur:
                    lens.append(cur)
                cur = 0
            else:
                s = line.strip().upper()
                cur += len(s)
                gc += s.count("G") + s.count("C")
                at += s.count("A") + s.count("T")
    if cur:
        lens.append(cur)
    lens.sort(reverse=True)
    total = sum(lens)
    acc = 0
    n50 = 0
    for L in lens:
        acc += L
        if acc >= total / 2:
            n50 = L
            break
    gc_pct = (100.0 * gc / (gc + at)) if (gc + at) else ""
    return len(lens), total, n50, gc_pct


def locate(accession, biosample, listed_path):
    """Return (path, how) or (None, reason). Never substitutes a file whose
    identity is not verified by hash (hash check happens afterwards)."""
    cands = []
    if listed_path:
        cands.append((Path(listed_path) if os.path.isabs(listed_path) else ROOT / listed_path, "freeze_listed_path"))
    for d in SEARCH_DIRS:
        cands.append((d / f"{accession}.fna", "accession_named_cache"))
        cands.append((d / f"{accession}.fasta", "accession_named_cache"))
    if biosample:
        cands.append((ROOT / "Genomes" / f"{biosample}.fna", "genomes_biosample_named"))
    seen = set()
    for p, how in cands:
        if str(p) in seen:
            continue
        seen.add(str(p))
        if p.is_file():
            return str(p), how
    # slow fallback: index of every fasta under the retrieval dir, built once
    for p in accession_index().get(accession, []):
        return str(p), "index_accession"
    return None, "not_found"


def main():
    freeze = read_tsv(V7 / "04_cohorts_v7/fasta_hashes_v7.tsv")
    assert len(freeze) == 212, f"expected 212 rows, got {len(freeze)}"

    meta = {}
    for mf in ("metadata_annotation_union_v7.tsv", "metadata_expanded_reference_sensitivity_v7.tsv"):
        for r in read_tsv(V7 / "04_cohorts_v7" / mf):
            meta.setdefault(r["strain"], r)

    prior = {}
    if COVAR_FILE.exists():
        for r in read_tsv(COVAR_FILE):
            prior[r["accession"]] = r

    rows, failures = [], []
    for fr in freeze:
        acc, expected = fr["accession"], fr["sha256"]
        m = meta.get(acc, {})
        biosample = m.get("biosample", "")
        path, how = locate(acc, biosample, fr.get("fasta_path", ""))
        rec = {
            "accession": acc, "biosample": biosample, "path": path or "",
            "locate_method": how, "expected_sha256": expected,
            "actual_sha256": "", "match": "", "n_contigs": "",
            "total_length_bp": "", "n50_bp": "", "gc_percent": "",
            "assembly_level": m.get("assembly_level", ""),
        }
        if path is None:
            rec["match"] = "FAIL_NOT_FOUND"
            failures.append({**rec, "reason": "file not found in any search location"})
        else:
            actual = sha256_file(path)
            rec["actual_sha256"] = actual
            if actual != expected:
                rec["match"] = "FAIL_HASH_MISMATCH"
                failures.append({**rec, "reason": f"sha256 mismatch (file located via {how})"})
            else:
                rec["match"] = "PASS"
                n, tot, n50, gcp = fasta_stats(path)
                rec.update(n_contigs=n, total_length_bp=tot, n50_bp=n50,
                           gc_percent=f"{gcp:.4f}" if gcp != "" else "")
        rows.append(rec)

    # covariate comparison (only for PASS rows): vs fix2 table and v7 metadata
    cov_rows = []
    for rec in rows:
        acc = rec["accession"]
        m = meta.get(acc, {})
        p = prior.get(acc, {})
        def diff(a, b):
            return "" if a == "" or b == "" else ("same" if str(a) == str(b) else "DIFF")
        cov_rows.append({
            "accession": acc,
            "binding_match": rec["match"],
            "n_contigs_p1": rec["n_contigs"],
            "n_contigs_v7meta": m.get("n_contigs", ""),
            "n_contigs_fix2": p.get("n_contigs", ""),
            "n_contigs_diff_flag": "yes" if "DIFF" in (diff(rec["n_contigs"], m.get("n_contigs", "")), diff(rec["n_contigs"], p.get("n_contigs", ""))) else "no",
            "total_length_bp_p1": rec["total_length_bp"],
            "genome_size_v7meta": m.get("genome_size", ""),
            "genome_size_fix2": p.get("genome_size_bp", ""),
            "size_diff_flag": "yes" if "DIFF" in (diff(rec["total_length_bp"], m.get("genome_size", "")), diff(rec["total_length_bp"], p.get("genome_size_bp", ""))) else "no",
            "n50_bp_p1": rec["n50_bp"],
            "n50_v7meta": m.get("n50_bp", ""),
            "n50_fix2": p.get("n50_bp", ""),
            "n50_diff_flag": "yes" if "DIFF" in (diff(rec["n50_bp"], m.get("n50_bp", "")), diff(rec["n50_bp"], p.get("n50_bp", ""))) else "no",
            "gc_percent_p1": rec["gc_percent"],
            "fasta_sha256": rec["actual_sha256"],
            "sha256_fix2": p.get("fasta_sha256", ""),
            "sha256_fix2_flag": diff(rec["actual_sha256"], p.get("fasta_sha256", "")),
            "fasta_path": rec["path"],
        })

    bfields = ["accession", "biosample", "path", "locate_method", "expected_sha256",
               "actual_sha256", "match", "n_contigs", "total_length_bp", "n50_bp",
               "gc_percent", "assembly_level"]
    with open(OUT / "fasta_binding_212.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=bfields, delimiter="\t")
        w.writeheader(); w.writerows(rows)
    with open(OUT / "assembly_covariates_212.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(cov_rows[0].keys()), delimiter="\t")
        w.writeheader(); w.writerows(cov_rows)
    ffields = bfields + ["reason"]
    with open(OUT / "fasta_binding_failures.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=ffields, delimiter="\t")
        w.writeheader(); w.writerows(failures)

    npass = sum(1 for r in rows if r["match"] == "PASS")
    ndiff = sum(1 for r in cov_rows if "yes" in (r["n_contigs_diff_flag"], r["size_diff_flag"], r["n50_diff_flag"]) or r["sha256_fix2_flag"] == "DIFF")
    print(f"binding: {npass}/212 PASS, {len(failures)} FAIL")
    print(f"covariate rows with any diff flag: {ndiff}")
    for f in failures:
        print("FAIL", f["accession"], f["biosample"], f["reason"])


if __name__ == "__main__":
    main()
