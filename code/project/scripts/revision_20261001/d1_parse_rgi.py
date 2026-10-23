#!/usr/bin/env python3
"""D1: parse RGI main outputs for tet(M) status across all genomes.

Status rules (per revision task spec):
  positive_complete  - RGI hit on tet(M) with Cut_Off in {Perfect, Strict} and
                       Percentage Length of Reference Sequence >= 95
  indeterminate      - any tet(M) hit that is partial (coverage < 95), Loose-only,
                       or conflicting (multiple conflicting hits)
  negative           - RGI ran successfully, no tet(M) hit at all
  missing            - RGI output absent / empty / run failure

Writes tetm_rgi_status.tsv: accession, status, best_cutoff, best_identity,
best_coverage_pct, n_hits, detail.
"""
import os
import csv, sys
from pathlib import Path

ROOT = Path(os.environ["NEISSERIA_PROJECT_ROOT"])
RGI_DIR = ROOT / "results_revision_20261001/08_amr_distribution/main"
LIST = sys.argv[1] if len(sys.argv) > 1 else "/tmp/rgi_list.txt"
OUT = ROOT / "results_revision_20261001/07_amr_tetm/tetm_rgi_status.tsv"

accs = [l.split()[0] for l in open(LIST) if l.strip() and not l.startswith("accession")]

def parse(acc):
    f = RGI_DIR / acc / f"{acc}.txt"
    if not f.exists() or f.stat().st_size == 0:
        return dict(status="missing", best_cutoff="", best_identity="", best_cov="", n_hits=0, detail="rgi_output_absent")
    rows = []
    with open(f) as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        if rd.fieldnames is None:
            return dict(status="missing", best_cutoff="", best_identity="", best_cov="", n_hits=0, detail="rgi_output_empty")
        for r in rd:
            fam = (r.get("AMR Gene Family") or "")
            hit = (r.get("Best_Hit_ARO") or "")
            if "tet(M)" in fam or "tet(M)" in hit:
                rows.append(r)
    if not rows:
        return dict(status="negative", best_cutoff="", best_identity="", best_cov="", n_hits=0, detail="no_tetM_hit")
    # rank hits: Perfect > Strict > Loose; then coverage
    rank = {"Perfect": 3, "Strict": 2, "Loose": 1}
    def key(r):
        cov = float(r.get("Percentage Length of Reference Sequence") or 0)
        return (rank.get(r.get("Cut_Off", ""), 0), cov, float(r.get("Best_Identities") or 0))
    rows.sort(key=key, reverse=True)
    best = rows[0]
    cov = float(best.get("Percentage Length of Reference Sequence") or 0)
    cut = best.get("Cut_Off", "")
    ident = best.get("Best_Identities", "")
    complete_hits = [r for r in rows if r.get("Cut_Off") in ("Perfect", "Strict")
                     and float(r.get("Percentage Length of Reference Sequence") or 0) >= 95]
    if complete_hits:
        status = "positive_complete"
    else:
        status = "indeterminate"
    detail = ";".join(f"{r.get('Cut_Off')}:{r.get('Best_Identities')}:{r.get('Percentage Length of Reference Sequence')}"
                      for r in rows[:4])
    return dict(status=status, best_cutoff=cut, best_identity=ident,
                best_cov=cov, n_hits=len(rows), detail=detail)

with open(OUT, "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["accession", "tetM_status", "best_cutoff", "best_identity_pct",
                "best_coverage_pct", "n_tetM_hits", "detail"])
    counts = {}
    for acc in accs:
        r = parse(acc)
        counts[r["status"]] = counts.get(r["status"], 0) + 1
        w.writerow([acc, r["status"], r["best_cutoff"], r["best_identity"],
                    r["best_cov"], r["n_hits"], r["detail"]])
print("status_counts:", counts)
print("wrote", OUT)
