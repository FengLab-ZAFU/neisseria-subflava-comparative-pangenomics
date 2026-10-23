#!/usr/bin/env python3
"""M5b: COG composition (integer detection + 1/k fractional weight) and
cluster-level Fisher enrichment for families COG_partition and
COG_species_associated. Reads fisher_output_cog.tsv produced by m5_04_fisher.R
from fisher_input_cog.tsv written here on first pass (--prepare), then
--finalize writes all tables.
"""
from pathlib import Path
import csv, json, os, sys
from collections import defaultdict

ROOT = str(Path(os.environ["NEISSERIA_PROJECT_ROOT"]))
RES = f"{ROOT}/results_revision_20261015_main152_species"
OUT = f"{RES}/06_function"
PGAP2 = f"{RES}/02_pangenome/pgap2_main152"
GROUPS = ["subflava", "perflava", "flavescens"]
PARTITIONS = ["Strict_core", "Core", "Soft_core", "Shell", "Cloud"]

def bh(pvals):
    m = len(pvals)
    order = sorted(range(m), key=lambda i: pvals[i])
    q = [0.0] * m
    prev = 1.0
    for rank in range(m, 0, -1):
        i = order[rank - 1]
        prev = min(prev, pvals[i] * m / rank)
        q[i] = prev
    return q

def load_inputs():
    partition = {}
    for r in csv.DictReader(open(f"{PGAP2}/pav_tables/partition_table.tsv"), delimiter="\t"):
        partition[r["cluster"]] = r["partition_recomputed"]
    cmap = {}
    for r in csv.DictReader(open(f"{OUT}/cluster_cog_map.tsv"), delimiter="\t"):
        cmap[r["cluster"]] = r
    # eligible Shell+Cloud (carrier>=5)
    eligible = set()
    for r in csv.DictReader(open(f"{OUT}/cluster_filter_ledger.tsv"), delimiter="\t"):
        if r["disposition"] == "tested":
            eligible.add(r["cluster"])
    cand_by_group = defaultdict(set)
    cp = f"{OUT}/associated_candidates_ranked.tsv"
    if os.path.exists(cp):
        for r in csv.DictReader(open(cp), delimiter="\t"):
            cand_by_group[r["target_group"]].add(r["cluster"])
    return partition, cmap, eligible, cand_by_group

def build_tests():
    """Return list of dicts: family_id, set_id, category, a, b, c, d, fg_ids, bg_ids, status, reason."""
    partition, cmap, eligible, cand_by_group = load_inputs()
    annotated = {cl for cl, r in cmap.items() if r["annotation_state"] == "annotated"}
    letters = {cl: set(cmap[cl]["cog_letters"]) - {"-"} for cl in annotated}

    tests = []

    # ---- family COG_partition: foreground = annotated clusters of the partition;
    # background = all other annotated main clusters.
    for p in PARTITIONS:
        fg = {cl for cl in annotated if partition[cl] == p}
        bg = annotated - fg
        cats = sorted({l for cl in fg | bg for l in letters[cl]})
        if not fg or not bg:
            tests.append(dict(family_id="COG_partition", set_id=p, category="NA",
                              a=0, b=0, c=0, d=0, fg_ids=sorted(fg), bg_ids=sorted(bg),
                              status="NOT_ESTIMABLE", reason="empty foreground or background"))
            continue
        for cat in cats:
            a = sum(1 for cl in fg if cat in letters[cl])
            c = sum(1 for cl in bg if cat in letters[cl])
            tests.append(dict(family_id="COG_partition", set_id=p, category=cat,
                              a=a, b=len(fg) - a, c=c, d=len(bg) - c,
                              fg_ids=sorted(fg), bg_ids=sorted(bg),
                              status="TESTED", reason=""))

    # ---- family COG_species_associated: background = eligible annotated
    # Shell+Cloud; foreground = associated candidates of the group (annotated).
    bg_all = eligible & annotated
    for g in GROUPS:
        fg = cand_by_group.get(g, set()) & annotated
        set_id = f"assoc_{g}"
        if not fg:
            cats = sorted({l for cl in bg_all for l in letters[cl]})
            for cat in cats:
                tests.append(dict(family_id="COG_species_associated", set_id=set_id,
                                  category=cat, a=0, b=0, c=0, d=0,
                                  fg_ids=[], bg_ids=sorted(bg_all),
                                  status="NOT_ESTIMABLE",
                                  reason="empty foreground (no annotated associated candidates)"))
            continue
        bg = bg_all - fg
        cats = sorted({l for cl in fg | bg for l in letters[cl]})
        for cat in cats:
            a = sum(1 for cl in fg if cat in letters[cl])
            c = sum(1 for cl in bg if cat in letters[cl])
            tests.append(dict(family_id="COG_species_associated", set_id=set_id,
                              category=cat, a=a, b=len(fg) - a, c=c, d=len(bg) - c,
                              fg_ids=sorted(fg), bg_ids=sorted(bg),
                              status="TESTED", reason=""))
    return tests, partition, cmap, annotated, letters, eligible, cand_by_group

def prepare():
    tests, *_ = build_tests()
    with open(f"{OUT}/fisher_input_cog.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["job_id", "a", "b", "c", "d"])
        n = 0
        for i, t in enumerate(tests):
            if t["status"] == "TESTED":
                w.writerow([f"cog_{i}", t["a"], t["b"], t["c"], t["d"]])
                n += 1
    sys.stderr.write(f"fisher jobs written: {n} / {len(tests)} rows\n")

def finalize():
    tests, partition, cmap, annotated, letters, eligible, cand_by_group = build_tests()
    fisher = {r["job_id"]: r for r in csv.DictReader(open(f"{OUT}/fisher_output_cog.tsv"), delimiter="\t")}

    # BH per family over valid tests
    fam_idx = defaultdict(list)
    for i, t in enumerate(tests):
        if t["status"] == "TESTED":
            fam_idx[t["family_id"]].append(i)
    qmap = {}
    for fam, idxs in fam_idx.items():
        qs = bh([float(fisher[f"cog_{i}"]["p"]) for i in idxs])
        for i, q in zip(idxs, qs):
            qmap[i] = q

    os.makedirs(f"{OUT}/foreground_background_sets", exist_ok=True)
    set_written = set()
    with open(f"{OUT}/COG_enrichment_all.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["family_id", "set_id", "category", "a", "b", "c", "d",
                    "OR", "CI_low", "CI_high", "p", "q", "status", "reason"])
        for i, t in enumerate(tests):
            key = (t["family_id"], t["set_id"])
            if key not in set_written:
                set_written.add(key)
                with open(f"{OUT}/foreground_background_sets/{t['family_id']}__{t['set_id']}.tsv", "w", newline="") as sf:
                    sw = csv.writer(sf, delimiter="\t", lineterminator="\n")
                    sw.writerow(["role", "cluster"])
                    for cl in t["fg_ids"]:
                        sw.writerow(["foreground", cl])
                    for cl in t["bg_ids"]:
                        sw.writerow(["background", cl])
            if t["status"] == "TESTED":
                f = fisher[f"cog_{i}"]
                w.writerow([t["family_id"], t["set_id"], t["category"], t["a"], t["b"],
                            t["c"], t["d"], f["OR"], f["CI_low"], f["CI_high"], f["p"],
                            f"{qmap[i]:.10g}", "TESTED", ""])
            else:
                w.writerow([t["family_id"], t["set_id"], t["category"], t["a"], t["b"],
                            t["c"], t["d"], "NA", "NA", "NA", "NA", "NA",
                            t["status"], t["reason"]])

    # ---- composition table
    sets = {}
    for p in PARTITIONS:
        sets[p] = sorted(cl for cl in cmap if partition[cl] == p)
    for g in GROUPS:
        sets[f"assoc_{g}"] = sorted(cand_by_group.get(g, set()))
    sets["all_main"] = sorted(cmap)
    all_cats = sorted({l for cl in annotated for l in letters[cl]})
    with open(f"{OUT}/COG_composition.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["set_id", "category", "integer_detected_clusters", "fractional_weight",
                    "annotated_background", "unannotated_count", "set_size"])
        for sid, members in sets.items():
            ann = [cl for cl in members if cl in annotated]
            unann = len(members) - len(ann)
            for cat in all_cats:
                integ = sum(1 for cl in ann if cat in letters[cl])
                frac = sum(1.0 / len(letters[cl]) for cl in ann if cat in letters[cl])
                w.writerow([sid, cat, integ, f"{frac:.6f}", len(ann), unann, len(members)])

    summary = {"n_enrichment_rows": len(tests),
               "n_tested": sum(1 for t in tests if t["status"] == "TESTED"),
               "n_not_estimable": sum(1 for t in tests if t["status"] != "TESTED"),
               "categories_detected": all_cats,
               "annotated_clusters": len(annotated),
               "composition_sets": {k: len(v) for k, v in sets.items()}}
    json.dump(summary, open(f"{OUT}/cog_enrichment_summary.json", "w"), indent=1)
    sys.stderr.write(json.dumps({k: v for k, v in summary.items() if k != "categories_detected"}) + "\n")

if __name__ == "__main__":
    if "--prepare" in sys.argv:
        prepare()
    elif "--finalize" in sys.argv:
        finalize()
    else:
        sys.exit("usage: m5_06_cog_enrichment.py --prepare|--finalize")
