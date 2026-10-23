#!/usr/bin/env python3
"""M8 (main152): cluster-level gene-MGE/GI overlap + positional null model.

Implements 执行细则 §6 exactly:
  - Local MGE union per contig (ISEScan IS + IntegronFinder integron +
    geNomad provirus), merged; whole-contig predictions excluded, reported
    separately.
  - gene-MGE two tiers: any_bp (overlap_union_bp > 0; touching endpoints do
    NOT count) and half_gene (overlap_union_bp / gene_length >= 0.5).
  - gene-GI descriptive: overlap_union_bp / gene_length >= 0.5.
  - GI-MGE descriptive: overlap_union_bp / GI_length >= 0.5.
  - Cluster-level hit = at least one qualifying instance; three denominators
    (all / localizable / testable-null) reported separately.
  - Positional null: per contig, merged MGE layout (count/lengths/relative
    arrangement fixed) shifted by a uniformly drawn integer offset; genes/GIs
    fixed. Linear/unknown contig length L, layout span [s_min, e_max): offset
    uniform on [-s_min, L - e_max] inclusive. Circular replicons (explicit
    Bakta [topology=circular] evidence only): offset uniform on 0..L-1,
    wrap-around split into two segments. Contigs independent; per-iteration
    and per-contig sub-seeds derived deterministically from master seed
    20261015 (saved per iteration).
  - Layouts with exactly one legal offset are degenerate: excluded from the
    testable instance set; observed recomputed on the same subset.
    MGE-free contigs with locatable genes are kept as fixed zero background.
  - 1000 iterations; upper-tail empirical P = (1 + #{null >= obs}) / 1001.
    BH family interval_main_and_species = pooled + 3 species strata x
    5 global partitions x 2 criteria, all testable combos.

Parallelism: multiprocessing over accessions, <=32 workers. RNG is fully
deterministic per (iteration, accession, contig) and independent of
scheduling order.

Usage: m8_overlap_null.py <tables_dir> <metadata.tsv> <outdir> [n_workers]
"""
import csv
import hashlib
import json
import os
import sys
import time
from collections import defaultdict
from multiprocessing import Pool

import numpy as np

MASTER_SEED = 20261015
N_ITER = 1000
RUN_ID = "main152_m8_20261015"
FAMILY_ID = "interval_main_and_species"
PARTITIONS = ["Strict_core", "Core", "Soft_core", "Shell", "Cloud"]
GROUPS = ["pooled", "subflava", "perflava", "flavescens"]
CRITERIA = ["any_bp", "half_gene"]


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def union_iv(intervals):
    iv = sorted(intervals)
    out = []
    for s, e in iv:
        if out and s <= out[-1][1]:
            if e > out[-1][1]:
                out[-1][1] = e
        else:
            out.append([s, e])
    return [(s, e) for s, e in out]


def iter_seed(i):
    h = hashlib.sha256(f"{MASTER_SEED}|{RUN_ID}|iter|{i}".encode()).digest()
    return int.from_bytes(h[:8], "little")


def contig_seed(isd, acc, contig):
    h = hashlib.sha256(f"{isd}|{acc}|{contig}".encode()).digest()
    return int.from_bytes(h[:8], "little")


def draw_offsets(acc, contig, kind, dmin, dmax, L):
    """1000 offsets for one contig, deterministic."""
    out = np.empty(N_ITER, dtype=np.int64)
    for i in range(1, N_ITER + 1):
        sd = contig_seed(iter_seed(i), acc, contig)
        # splitmix64-style stream from seed; single uniform draw per iteration
        rng = sd & 0xFFFFFFFFFFFFFFFF
        rng = (rng + 0x9E3779B97F4A7C15) & 0xFFFFFFFFFFFFFFFF
        z = rng
        z = ((z ^ (z >> 30)) * 0xBF58476D1CE4E5B9) & 0xFFFFFFFFFFFFFFFF
        z = ((z ^ (z >> 27)) * 0x94D049BB133111EB) & 0xFFFFFFFFFFFFFFFF
        z = z ^ (z >> 31)
        if kind == "circular":
            out[i - 1] = z % L
        else:
            out[i - 1] = dmin + z % (dmax - dmin + 1)
    return out


# ---------------- per-accession worker ----------------
def process_accession(task):
    (acc, contigs_info, genes_acc) = task
    """contigs_info: {contig: (L, topology, U_tuple, dmin, dmax, status)}
       genes_acc: list of dicts with start/end/cluster (valid instances)."""
    per_cluster = {}   # cluster -> dict with arrays
    inst_rows = []     # per-instance observed overlaps
    by_ctg = defaultdict(list)
    for g in genes_acc:
        by_ctg[g["contig_id"]].append(g)
    for ctg, glist in sorted(by_ctg.items()):
        L, topo, U, dmin, dmax, status = contigs_info[ctg]
        # observed overlaps (real positions), all instances
        for g in glist:
            gs, ge = g["start"], g["end"]
            obs_ov = 0
            for s, e in U:
                ov = min(ge, e) - max(gs, s)
                if ov > 0:
                    obs_ov += ov
            eligible = status != "DEGENERATE_LAYOUT"
            inst_rows.append((g["gene_id"], obs_ov, obs_ov > 0,
                              2 * obs_ov >= ge - gs, eligible, status))
        if status != "OK":
            continue
        d = draw_offsets(acc, ctg,
                         "circular" if topo == "circular" else "linear",
                         dmin, dmax, L)
        U_arr = np.asarray(U, dtype=np.int64)
        for g in glist:
            gs, ge = g["start"], g["end"]
            glen = ge - gs
            tot = np.zeros(N_ITER, dtype=np.int64)
            if topo == "circular":
                for s, e in U_arr:
                    ns = (s + d) % L
                    ne = ns + (e - s)
                    seg1 = np.minimum(ge, np.minimum(ne, L)) - np.maximum(gs, ns)
                    np.clip(seg1, 0, None, out=seg1)
                    seg2 = np.minimum(ge, ne - L) - gs
                    np.clip(seg2, 0, None, out=seg2)
                    # where ne <= L, ne-L <= 0 -> seg2 clipped to 0; safe
                    tot += seg1 + seg2
            else:
                for s, e in U_arr:
                    ov = np.minimum(ge, e + d) - np.maximum(gs, s + d)
                    np.clip(ov, 0, None, out=ov)
                    tot += ov
            ent = per_cluster.setdefault(g["cluster"], {
                "obs_any": False, "obs_half": False,
                "any": np.zeros(N_ITER, dtype=bool),
                "half": np.zeros(N_ITER, dtype=bool)})
            ent["any"] |= tot > 0
            ent["half"] |= 2 * tot >= glen
    # merge observed flags for eligible instances
    obs_map = {gid: (ov, a, h, el) for gid, ov, a, h, el, _st in inst_rows}
    for g in genes_acc:
        ov, a, h, el = obs_map[g["gene_id"]]
        if not el:
            continue
        ent = per_cluster.setdefault(g["cluster"], {
            "obs_any": False, "obs_half": False,
            "any": np.zeros(N_ITER, dtype=bool),
            "half": np.zeros(N_ITER, dtype=bool)})
        ent["obs_any"] |= a
        ent["obs_half"] |= h
    return acc, per_cluster, inst_rows


def main():
    tables, meta_file, outdir = sys.argv[1], sys.argv[2], sys.argv[3]
    n_workers = min(32, int(sys.argv[4]) if len(sys.argv) > 4 else 32)
    os.makedirs(outdir, exist_ok=True)
    t0 = time.time()

    meta = {r["strain"]: r for r in read_tsv(meta_file)}
    group_of = {a: meta[a]["group_for_testing"] for a in meta}

    contigs = {(r["accession"], r["contig_id"]): r for r in read_tsv(os.path.join(tables, "contigs.tsv"))}
    genes = [r for r in read_tsv(os.path.join(tables, "genes.tsv"))]
    mge = [r for r in read_tsv(os.path.join(tables, "mge_local.tsv"))]
    gis = [r for r in read_tsv(os.path.join(tables, "gi.tsv"))]
    wc = [r for r in read_tsv(os.path.join(tables, "whole_contig.tsv"))]

    # ---- contig layout info ----
    mge_by_ctg = defaultdict(list)
    for r in mge:
        mge_by_ctg[(r["accession"], r["contig_id"])].append((int(r["start"]), int(r["end"])))

    contig_info = {}   # (acc,ctg) -> (L, topo, U, dmin, dmax, status)
    degenerate = []
    for key, c in contigs.items():
        L = int(c["length_bp"])
        topo = c["topology"]
        U = union_iv(mge_by_ctg.get(key, []))
        if not U:
            contig_info[key] = (L, topo, tuple(), 0, 0, "NO_LOCAL_MGE")
            continue
        s_min = min(s for s, _ in U)
        e_max = max(e for _, e in U)
        if topo == "circular":
            n_off = L
            dmin, dmax = 0, L - 1
        else:
            dmin, dmax = -s_min, L - e_max
            n_off = dmax - dmin + 1
        if n_off <= 1:
            contig_info[key] = (L, topo, tuple(U), dmin, dmax, "DEGENERATE_LAYOUT")
            degenerate.append((key[0], key[1], topo, L, len(U), s_min, e_max, n_off))
        else:
            contig_info[key] = (L, topo, tuple(U), dmin, dmax, "OK")

    # ---- per-accession gene task payloads ----
    genes_by_acc = defaultdict(list)
    for r in genes:
        if r["valid"] != "1":
            continue
        r["start"] = int(r["start"]); r["end"] = int(r["end"])
        genes_by_acc[r["accession"]].append(r)

    tasks = []
    for acc in sorted(genes_by_acc):
        ci = {ctg: contig_info[(acc, ctg)] for ctg in
              {g["contig_id"] for g in genes_by_acc[acc]}}
        tasks.append((acc, ci, genes_by_acc[acc]))

    # ---- run workers ----
    results = {}
    with Pool(n_workers) as pool:
        for acc, per_cluster, inst_rows in pool.imap_unordered(process_accession, tasks):
            results[acc] = (per_cluster, inst_rows)

    # ---- per-instance observed table ----
    gi_by_ctg = defaultdict(list)
    for r in gis:
        gi_by_ctg[(r["accession"], r["contig_id"])].append((int(r["start"]), int(r["end"])))
    gi_union_by_ctg = {k: union_iv(v) for k, v in gi_by_ctg.items()}

    inst_lookup = {}   # (acc, gene_id) -> (obs_ov, obs_any, obs_half, eligible, status)
    for acc, (_, inst_rows) in results.items():
        for gid, ov, a, h, el, st in inst_rows:
            inst_lookup[(acc, gid)] = (ov, a, h, el, st)

    occ_fields = ["accession", "group", "gene_id", "cluster", "partition", "contig_id",
                  "coordinate_system", "start", "end", "strand", "gene_length",
                  "mge_overlap_union_bp", "mge_any_bp", "mge_half_gene",
                  "gi_overlap_union_bp", "gi_half_gene",
                  "eligible_for_null", "contig_layout_status",
                  "unmapped_invalid_reason",
                  "orig_contig_id", "orig_start", "orig_end", "orig_coordinate_system"]
    with open(os.path.join(outdir, "occurrence_overlap.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=occ_fields)
        w.writeheader()
        for r in genes:
            row = {k: r.get(k, "") for k in
                   ["accession", "gene_id", "cluster", "partition", "contig_id",
                    "start", "end", "strand", "orig_contig_id", "orig_start",
                    "orig_end", "orig_coordinate_system"]}
            row["group"] = group_of.get(r["accession"], "NA")
            row["coordinate_system"] = "0-based_half-open" if r["valid"] == "1" else "NA"
            row["unmapped_invalid_reason"] = r["invalid_reason"] if r["valid"] != "1" else ""
            if r["valid"] == "1":
                ov, a, h, el, st = inst_lookup[(r["accession"], r["gene_id"])]
                gs, ge = int(r["start"]), int(r["end"])
                row["gene_length"] = ge - gs
                row["mge_overlap_union_bp"] = ov
                row["mge_any_bp"] = "TRUE" if a else "FALSE"
                row["mge_half_gene"] = "TRUE" if h else "FALSE"
                gov = 0
                for s, e in gi_union_by_ctg.get((r["accession"], r["contig_id"]), []):
                    o = min(ge, e) - max(gs, s)
                    if o > 0:
                        gov += o
                row["gi_overlap_union_bp"] = gov
                row["gi_half_gene"] = "TRUE" if 2 * gov >= (ge - gs) else "FALSE"
                row["eligible_for_null"] = "TRUE" if el else "FALSE"
                row["contig_layout_status"] = st
            else:
                for k in ["gene_length", "mge_overlap_union_bp", "gi_overlap_union_bp"]:
                    row[k] = "NA"
                for k in ["mge_any_bp", "mge_half_gene", "gi_half_gene",
                          "eligible_for_null", "contig_layout_status"]:
                    row[k] = "NA"
            w.writerow(row)

    # ---- aggregate to cluster-scope level ----
    # scope membership: pooled = all accessions; group scope = its genomes only
    cl_scope = {}   # (cluster, scope) -> {obs_any, obs_half, any[], half[]}
    cl_partition = {}
    cl_detected = defaultdict(set)   # cluster -> set of scopes where detected (any instance)
    cl_localizable = defaultdict(set)
    cl_testable = defaultdict(set)
    for r in genes:
        cl = r["cluster"]
        cl_partition[cl] = r["partition"]
        scopes = ("pooled", group_of[r["accession"]])
        for sc in scopes:
            cl_detected[cl].add(sc)
        if r["valid"] == "1":
            for sc in scopes:
                cl_localizable[cl].add(sc)
            if inst_lookup[(r["accession"], r["gene_id"])][3]:
                for sc in scopes:
                    cl_testable[cl].add(sc)

    for acc, (per_cluster, _) in results.items():
        for cl, ent in per_cluster.items():
            for sc in ("pooled", group_of[acc]):
                agg = cl_scope.setdefault((cl, sc), {
                    "obs_any": False, "obs_half": False,
                    "any": np.zeros(N_ITER, dtype=bool),
                    "half": np.zeros(N_ITER, dtype=bool)})
                agg["obs_any"] |= ent["obs_any"]
                agg["obs_half"] |= ent["obs_half"]
                agg["any"] |= ent["any"]
                agg["half"] |= ent["half"]

    # ---- overlap_summary.tsv ----
    sum_rows = []
    for grp in GROUPS:
        for p in PARTITIONS:
            cls_all = [c for c in cl_partition
                       if cl_partition[c] == p and grp in cl_detected.get(c, ())]
            cls_loc = [c for c in cls_all if grp in cl_localizable.get(c, ())]
            cls_test = [c for c in cls_all if grp in cl_testable.get(c, ())]
            for crit in CRITERIA:
                obs_key = "obs_any" if crit == "any_bp" else "obs_half"
                for dname, dset in (("all_clusters", cls_all),
                                    ("localizable_clusters", cls_loc),
                                    ("testable_null_clusters", cls_test)):
                    if dname == "testable_null_clusters":
                        n_hit = sum(1 for c in dset
                                    if cl_scope.get((c, grp), {}).get(obs_key, False))
                    else:
                        # hit over all localizable instances (real positions)
                        n_hit = 0
                        for c in dset:
                            hit = False
                            for r in genes_by_cluster[c]:
                                if r["valid"] != "1":
                                    continue
                                if grp != "pooled" and group_of[r["accession"]] != grp:
                                    continue
                                ov, a, h, el, st = inst_lookup[(r["accession"], r["gene_id"])]
                                if (a if crit == "any_bp" else h):
                                    hit = True
                                    break
                            n_hit += int(hit)
                    n_tot = len(dset)
                    sum_rows.append({
                        "group": grp, "partition": p, "criterion": crit,
                        "denominator_type": dname, "n_hit": n_hit, "n_total": n_tot,
                        "proportion": (n_hit / n_tot) if n_tot else "NA",
                        "status": "OK" if n_tot else "NO_CLUSTERS_IN_DENOMINATOR"})

    # ---- null model per combo ----
    combo_rows = []
    perm_rows = []
    elig_rows = []
    combos_testable = []
    for grp in GROUPS:
        for p in PARTITIONS:
            cls_all = [c for c in cl_partition
                       if cl_partition[c] == p and grp in cl_detected.get(c, ())]
            cls_test = sorted(c for c in cls_all if grp in cl_testable.get(c, ()))
            n_deg = sum(1 for (a, ct, *_r) in degenerate if grp == "pooled" or group_of[a] == grp)
            for crit in CRITERIA:
                key = (grp, p, crit)
                elig_rows.append({
                    "group": grp, "partition": p, "criterion": crit,
                    "n_detected_clusters": len(cls_all),
                    "n_localizable_clusters": sum(1 for c in cls_all if grp in cl_localizable.get(c, ())),
                    "n_testable_clusters": len(cls_test),
                    "n_degenerate_layout_contigs_in_scope": n_deg,
                    "testable": "TRUE" if cls_test else "FALSE",
                    "reason": "" if cls_test else "no cluster with >=1 eligible instance in scope"})
                if not cls_test:
                    combo_rows.append({"group": grp, "partition": p, "criterion": crit,
                                       "n_eligible": 0, "observed": "NA", "null_mean": "NA",
                                       "null_q025": "NA", "null_q975": "NA",
                                       "effect_ratio": "NA", "p": "NA", "q": "NA",
                                       "family_id": FAMILY_ID, "status": "NOT_ESTIMABLE",
                                       "reason": "no testable cluster (no eligible instance in scope)"})
                    continue
                obs_key = "obs_any" if crit == "any_bp" else "obs_half"
                arr_key = "any" if crit == "any_bp" else "half"
                mat = np.stack([cl_scope[(c, grp)][arr_key] for c in cls_test])
                hits = mat.sum(axis=0)
                obs = sum(1 for c in cls_test if cl_scope[(c, grp)][obs_key])
                combos_testable.append((key, cls_test, obs, hits))

    # seeds per iteration
    seeds = {i: iter_seed(i) for i in range(1, N_ITER + 1)}
    for (grp, p, crit), cls_test, obs, hits in combos_testable:
        for i in range(1, N_ITER + 1):
            perm_rows.append({"group": grp, "partition": p, "criterion": crit,
                              "iteration": i, "seed": seeds[i],
                              "n_eligible": len(cls_test), "hits": int(hits[i - 1])})

    # stats + BH
    valid = []
    for (grp, p, crit), cls_test, obs, hits in combos_testable:
        B = N_ITER
        mean = float(hits.mean())
        srt = sorted(int(x) for x in hits)
        q025 = srt[int(0.025 * (B - 1))]
        q975 = srt[int(0.975 * (B - 1))]
        p_emp = (1 + int((hits >= obs).sum())) / (B + 1)
        if mean > 0:
            eff = obs / mean
            eff_s = repr(eff)
        elif obs > 0:
            eff_s = "Inf"
        else:
            eff_s = "NA"
        reason = ""
        if mean == 0:
            reason = "null_mean_zero;effect_ratio " + ("infinite" if obs > 0 else "undefined")
        if srt[0] == srt[-1]:
            reason = (reason + ";" if reason else "") + "degenerate_null_constant_distribution"
        row = {"group": grp, "partition": p, "criterion": crit,
               "n_eligible": len(cls_test), "observed": obs,
               "null_mean": mean, "null_q025": q025, "null_q975": q975,
               "effect_ratio": eff_s, "p": p_emp, "q": "NA",
               "family_id": FAMILY_ID, "status": "OK", "reason": reason}
        valid.append(row)
    m = len(valid)
    order = sorted(range(m), key=lambda i: valid[i]["p"])
    prev = 1.0
    qs = [1.0] * m
    for rank in range(m - 1, -1, -1):
        i = order[rank]
        prev = min(prev, valid[i]["p"] * m / (rank + 1))
        qs[i] = prev
    for i, row in enumerate(valid):
        row["q"] = qs[i]
    combo_rows.extend(valid)
    combo_rows.sort(key=lambda r: (GROUPS.index(r["group"]),
                                   PARTITIONS.index(r["partition"]),
                                   CRITERIA.index(r["criterion"])))

    # ---- write outputs ----
    def dump(name, rows, fields):
        with open(os.path.join(outdir, name), "w", newline="") as fh:
            w = csv.DictWriter(fh, delimiter="\t", fieldnames=fields)
            w.writeheader()
            for r in rows:
                w.writerow(r)

    dump("overlap_summary.tsv", sum_rows,
         ["group", "partition", "criterion", "denominator_type", "n_hit", "n_total",
          "proportion", "status"])
    dump("interval_null_results.tsv", combo_rows,
         ["group", "partition", "criterion", "n_eligible", "observed", "null_mean",
          "null_q025", "null_q975", "effect_ratio", "p", "q", "family_id", "status", "reason"])
    dump("null_permutation_counts.tsv", perm_rows,
         ["group", "partition", "criterion", "iteration", "seed", "n_eligible", "hits"])
    dump("null_eligibility.tsv", elig_rows,
         ["group", "partition", "criterion", "n_detected_clusters", "n_localizable_clusters",
          "n_testable_clusters", "n_degenerate_layout_contigs_in_scope", "testable", "reason"])

    # whole-contig summary (per prediction; contigs validated against map)
    wc_rows = []
    for r in wc:
        acc = r["accession"]
        c = contigs.get((acc, r["contig_id"]))
        wc_rows.append({"accession": acc, "group": group_of.get(acc, "NA"),
                        "contig_id": r["contig_id"],
                        "contig_length_bp": c["length_bp"] if c else "NA",
                        "contig_topology": c["topology"] if c else "NA",
                        "tool": r["tool"], "prediction": r["prediction"],
                        "detail": r["detail"],
                        "status": "OK" if c else "CONTIG_NOT_IN_MAP"})
    dump("whole_contig_summary.tsv", wc_rows,
         ["accession", "group", "contig_id", "contig_length_bp", "contig_topology",
          "tool", "prediction", "detail", "status"])

    # GI-MGE descriptive
    mge_union_by_ctg = {k: union_iv(v) for k, v in mge_by_ctg.items()}
    gi_rows = []
    for r in gis:
        key = (r["accession"], r["contig_id"])
        gs, ge = int(r["start"]), int(r["end"])
        glen = ge - gs
        ov = 0
        for s, e in mge_union_by_ctg.get(key, []):
            o = min(ge, e) - max(gs, s)
            if o > 0:
                ov += o
        gi_rows.append({"accession": r["accession"], "group": group_of.get(r["accession"], "NA"),
                        "contig_id": r["contig_id"], "gi_id": r["gi_id"],
                        "start": gs, "end": ge, "gi_length": glen,
                        "mge_overlap_union_bp": ov,
                        "overlap_over_gi_length": ov / glen if glen else "NA",
                        "gi_half_covered_by_mge": "TRUE" if 2 * ov >= glen else "FALSE"})
    dump("gi_mge_overlap.tsv", gi_rows,
         ["accession", "group", "contig_id", "gi_id", "start", "end", "gi_length",
          "mge_overlap_union_bp", "overlap_over_gi_length", "gi_half_covered_by_mge"])

    dump("null_degenerate_contigs.tsv",
         [{"accession": a, "contig_id": ct, "topology": t, "length_bp": L,
           "n_merged_mge_intervals": n, "span_start": s, "span_end": e,
           "n_legal_offsets": no}
          for a, ct, t, L, n, s, e, no in sorted(degenerate)],
         ["accession", "contig_id", "topology", "length_bp", "n_merged_mge_intervals",
          "span_start", "span_end", "n_legal_offsets"])

    run_info = {
        "run_id": RUN_ID, "module": "M8_interval_overlap",
        "master_seed": MASTER_SEED, "n_iter": N_ITER, "family_id": FAMILY_ID,
        "n_workers": n_workers,
        "seed_derivation": "iter_seed=sha256('20261015|<run_id>|iter|<i>')[:8]le; "
                           "contig_seed=sha256('<iter_seed>|<acc>|<contig>')[:8]le; "
                           "offset=splitmix64(contig_seed) mod range",
        "n_combos_total": len(GROUPS) * len(PARTITIONS) * len(CRITERIA),
        "n_combos_testable": len(valid),
        "n_combos_not_estimable": len(combo_rows) - len(valid),
        "n_degenerate_layout_contigs": len(degenerate),
        "runtime_sec": round(time.time() - t0, 1),
    }
    with open(os.path.join(outdir, "run_info.json"), "w") as fh:
        json.dump(run_info, fh, indent=2)
    print(json.dumps(run_info, indent=2))


# module-level gene index for summary hits (built in main before use)
genes_by_cluster = defaultdict(list)

if __name__ == "__main__":
    # build genes_by_cluster inside main via a preload wrapper
    _orig_main = main

    def main_wrapper():
        tables = sys.argv[1]
        for r in read_tsv(os.path.join(tables, "genes.tsv")):
            genes_by_cluster[r["cluster"]].append(r)
        _orig_main()

    main_wrapper()
