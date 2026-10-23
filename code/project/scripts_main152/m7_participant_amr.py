#!/usr/bin/env python3
"""Participant-sensitivity: AMR/tet(M) carriage on the 121-genome representative
subset (complements M7 defense sensitivity; no PGAP2/tree rerun).
"""
import json
import re
import math
import os
import sys
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import RESULTS, read_tsv, write_tsv

OUT = os.path.join(RESULTS, "08_defense_mge_vfdb")


def wilson(k, n, z=1.959964):
    if n == 0:
        return ("NA", "NA", "NA")
    p = k / n
    d = 1 + z * z / n
    c = p + z * z / (2 * n)
    s = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    return (f"{p:.6f}", f"{(c - s) / d:.6f}", f"{(c + s) / d:.6f}")


def main():
    sub = read_tsv(os.path.join(OUT, "participant_subset_manifest.tsv"))
    selected = {r["accession"] for r in sub if str(r["selected"]) == "1"}
    assert len(selected) == 121
    pg = {r["accession"]: r for r in read_tsv(os.path.join(RESULTS, "07_amr_tetm/amr_tetm_per_genome.tsv"))}
    car = read_tsv(os.path.join(RESULTS, "07_amr_tetm/carriage_by_species.tsv"))
    features = sorted({r["feature"] for r in car})

    groups = {}
    for r in sub:
        if r["accession"] in selected:
            groups.setdefault(r["species_group"], []).append(r["accession"])

    rows = []
    for feat in features:
        is_tetm = feat == "gene:tet(M)"
        for g in ("subflava", "perflava", "flavescens"):
            accs = groups.get(g, [])
            n_total = len(accs)
            if is_tetm:
                n_pos = sum(1 for a in accs if pg[a]["tetM_status"] == "positive_complete")
                n_neg = sum(1 for a in accs if pg[a]["tetM_status"] == "negative")
                n_ind = sum(1 for a in accs if pg[a]["tetM_status"] == "indeterminate")
                n_mis = sum(1 for a in accs if pg[a]["tetM_status"] == "missing")
                n_ev = n_pos + n_neg
            else:
                fam = feat.split("family:", 1)[-1]
                n_pos = n_neg = 0
                n_ind = n_mis = 0
                for a in accs:
                    st = pg[a]["RGI_run_status"]
                    if st != "ok":
                        n_mis += 1
                        continue
                    hits = pg[a].get("amr_family_hitcounts_perfect_strict", "none")
                    # family names may contain '; ' themselves — match exact
                    # name + ':<count>' token instead of naive splitting
                    present = bool(re.search(re.escape(fam) + r":\d+", hits)) if hits not in ("", "none") else False
                    if present:
                        n_pos += 1
                    else:
                        n_neg += 1
                n_ev = n_pos + n_neg
            p, lo, hi = wilson(n_pos, n_ev)
            rows.append(dict(analysis_id="participant_sensitivity_121", feature=feat,
                             species_group=g, n_total=n_total, n_evaluable=n_ev,
                             n_positive=n_pos, n_negative=n_neg, n_indeterminate=n_ind,
                             n_missing=n_mis, prevalence=p, CI_low=lo, CI_high=hi,
                             status="TESTED" if n_ev else "NOT_ESTIMABLE",
                             reason="" if n_ev else "no evaluable genomes"))
    write_tsv(os.path.join(OUT, "sensitivity_amr_tetm.tsv"),
              ["analysis_id", "feature", "species_group", "n_total", "n_evaluable",
               "n_positive", "n_negative", "n_indeterminate", "n_missing",
               "prevalence", "CI_low", "CI_high", "status", "reason"], rows)

    # quick compare vs main
    cmp_rows = []
    main = {(r["feature"], r["species_group"]): r for r in car}
    for r in rows:
        m = main.get((r["feature"], r["species_group"]))
        if not m:
            continue
        same = "yes" if (int(m["n_total"]) == int(r["n_total"])) else "no"
        cmp_rows.append(dict(feature=r["feature"], species_group=r["species_group"],
                             main_n=m["n_total"], main_prev=m["prevalence"],
                             subset_n=r["n_total"], subset_prev=r["prevalence"],
                             identical_denominator=same))
    write_tsv(os.path.join(OUT, "sensitivity_amr_tetm_vs_main.tsv"),
              ["feature", "species_group", "main_n", "main_prev", "subset_n", "subset_prev",
               "identical_denominator"], cmp_rows)
    print(json.dumps(dict(status="PASS", rows=len(rows), subset_n=len(selected),
                          groups={g: len(v) for g, v in groups.items()}), indent=1))


if __name__ == "__main__":
    main()
