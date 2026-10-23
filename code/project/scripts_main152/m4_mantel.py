#!/usr/bin/env python3
"""M4c Mantel: descriptive Pearson correlation between new ML unrooted
patristic distances and Shell+Cloud Jaccard (999 whole-genome label
permutations). all_variable feature set as separately-labelled sensitivity.

Permutation rule: permute genome labels (not independent distance entries);
upper-tail two-sided descriptive correlation; MC p=(1+b)/(B+1) with b =
count(|r_perm| >= |r_obs|); report r, p, B. Not evidence controlling for
phylogenetic confounding; no root-direction interpretation.

Run: conda run -n gubbins python3 m4_mantel.py  (dendropy for patristic)
"""
import json
import os
import sys

import dendropy
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import RESULTS, MASTER_SEED, read_roster

CORE = os.path.join(RESULTS, "03_phylogeny/strict_singlecopy_core")
STR = os.path.join(RESULTS, "05_structure")
B = 999


def read_matrix(path):
    with open(path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")[1:]
        rows = []
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            rows.append([float(x) for x in f[1:]])
        order = [ln.split("\t")[0] for ln in open(path).read().splitlines()[1:]]
    return order, np.array(rows)


def upper_tri(m, n):
    i, j = np.triu_indices(n, k=1)
    return m[i, j]


def mantel(d1, d2, rng):
    n = d1.shape[0]
    x = upper_tri(d1, n)
    y = upper_tri(d2, n)
    if np.std(x) == 0 or np.std(y) == 0:
        return None, None, "constant distance vector"
    r_obs = np.corrcoef(x, y)[0, 1]
    b = 0
    for it in range(B):
        perm = rng.permutation(n)
        yp = upper_tri(d2[np.ix_(perm, perm)], n)
        rp = np.corrcoef(x, yp)[0, 1]
        if abs(rp) >= abs(r_obs):
            b += 1
    p = (1 + b) / (B + 1)
    return r_obs, p, None


def main():
    roster = read_roster()
    tree = dendropy.Tree.get(path=os.path.join(CORE, "main152_ML_unrooted.nwk"),
                             schema="newick", preserve_underscores=True)
    dm = tree.phylogenetic_distance_matrix()
    taxa = {t.label: t for t in tree.taxon_namespace}
    assert set(taxa) == set(roster)
    n = len(roster)
    pat = np.array([[dm.distance(taxa[a], taxa[b]) for b in roster] for a in roster])

    rng = np.random.default_rng(MASTER_SEED)
    rows = []
    for fs in ["shell_cloud", "all_variable"]:
        order, jac = read_matrix(os.path.join(STR, f"jaccard_distance_{fs}.tsv"))
        assert order == roster, f"{fs} matrix order != roster"
        assert np.allclose(jac, jac.T) and np.allclose(np.diag(jac), 0)
        r, p, note = mantel(pat, jac, rng)
        rows.append(dict(analysis_id=f"mantel_patristic_vs_{fs}", feature_set=fs,
                         test="Mantel_Pearson", permutation_design="free_whole_genome_labels",
                         n=n, statistic=r if r is None else f"{r:.6f}",
                         R2="NA", p=p if p is None else f"{p:.6f}", B_used=B,
                         status="TESTED" if r is not None else "NOT_ESTIMABLE",
                         reason=note or "",
                         interpretation="descriptive correlation only; not phylogenetic-confound control"))
        print(fs, "r =", r, "p =", p, note or "")

    out = os.path.join(STR, "mantel_results.tsv")
    with open(out, "w", newline="") as fh:
        import csv
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    print("written", out)

    # merge into distance_statistics.tsv (contract F3_DISTANCE_STATS); pad
    # Mantel-specific rows to its 16 columns, fold interpretation into reason
    ds_path = os.path.join(STR, "distance_statistics.tsv")
    with open(ds_path, newline="") as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        ds_cols = rdr.fieldnames
        ds_rows = [r for r in rdr if not r["analysis_id"].startswith("mantel_")]
    for r in rows:
        merged = {c: "NA" for c in ds_cols}
        merged["analysis_id"] = r["analysis_id"]
        merged["feature_set"] = r["feature_set"]
        merged["test"] = r["test"]
        merged["permutation_design"] = r["permutation_design"]
        merged["n"] = str(r["n"])
        merged["statistic"] = r["statistic"]
        merged["R2"] = r["R2"]
        merged["p"] = r["p"]
        merged["B_used"] = str(r["B_used"])
        merged["status"] = r["status"]
        merged["reason"] = (r["reason"] + "; " if r["reason"] else "") + r["interpretation"]
        ds_rows.append(merged)
    with open(ds_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=ds_cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(ds_rows)
    print("merged", len(rows), "Mantel rows into", ds_path, "total", len(ds_rows))


if __name__ == "__main__":
    main()
