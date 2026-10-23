#!/usr/bin/env python3
"""M8 (main152): independent coordinate/overlap audit with a naive algorithm.

Independent of the numpy implementation in m8_overlap_null.py: overlaps are
recomputed with literal integer position sets (or plain sorted-merge), seeds
are re-derived from the documented hash chain, and circular shifts are
evaluated on an explicit position circle.

Coverage:
  real instances (>=20): plus/minus strand, boundary-touch non-overlap,
  1-bp overlap, >=50% threshold boundary (both sides), same annotation
  contig name across accessions, genes on circular contigs with local MGE,
  random fillers (seeded, fixed).
  synthetic unit cases (clearly labelled, not biological data): boundary
  touch, 1-bp overlap, exact 50% threshold, overlapping-MGE union,
  circular wrap-around split, degenerate-layout exclusion.
  null-model end-to-end recheck: 3 iterations of the smallest testable
  combo recomputed from the saved per-iteration seeds with naive shifting.

Usage: m8_coordinate_audit.py <tables_dir> <outdir> <metadata.tsv>
"""
import csv
import hashlib
import json
import os
import random
import sys
from collections import defaultdict

MASTER_SEED = 20261015
RUN_ID = "main152_m8_20261015"
N_ITER = 1000
OUT_FIELDS = ["case_id", "case_kind", "edge_types", "accession", "feature_id",
              "contig_id", "gene_start", "gene_end", "strand",
              "main_overlap_bp", "naive_overlap_bp",
              "main_any_bp", "naive_any_bp", "main_half_gene", "naive_half_gene",
              "match", "notes"]


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def naive_overlap_bp(gs, ge, intervals):
    """Literal position-set overlap of [gs,ge) against union of intervals."""
    pos = set()
    for s, e in intervals:
        pos.update(range(s, e))
    return sum(1 for p in range(gs, ge) if p in pos)


def iter_seed(i):
    h = hashlib.sha256(f"{MASTER_SEED}|{RUN_ID}|iter|{i}".encode()).digest()
    return int.from_bytes(h[:8], "little")


def contig_offset(isd, acc, contig, kind, dmin, dmax, L):
    h = hashlib.sha256(f"{isd}|{acc}|{contig}".encode()).digest()
    rng = int.from_bytes(h[:8], "little") & 0xFFFFFFFFFFFFFFFF
    rng = (rng + 0x9E3779B97F4A7C15) & 0xFFFFFFFFFFFFFFFF
    z = rng
    z = ((z ^ (z >> 30)) * 0xBF58476D1CE4E5B9) & 0xFFFFFFFFFFFFFFFF
    z = ((z ^ (z >> 27)) * 0x94D049BB133111EB) & 0xFFFFFFFFFFFFFFFF
    z = z ^ (z >> 31)
    if kind == "circular":
        return z % L
    return dmin + z % (dmax - dmin + 1)


def shift_naive(U, d, L, circular):
    """Shift merged union by d; circular wraps around the position circle."""
    if not circular:
        return [(s + d, e + d) for s, e in U]
    out = []
    for s, e in U:
        ns = (s + d) % L
        ne = ns + (e - s)
        if ne <= L:
            out.append((ns, ne))
        else:
            out.append((ns, L))
            out.append((0, ne - L))
    return out


def main():
    tables, outdir, meta_file = sys.argv[1], sys.argv[2], sys.argv[3]
    meta = {r["strain"]: r for r in read_tsv(meta_file)}
    group_of = {a: meta[a]["group_for_testing"] for a in meta}
    genes = [r for r in read_tsv(os.path.join(tables, "genes.tsv"))]
    valid = [r for r in genes if r["valid"] == "1"]
    mge_by_ctg = defaultdict(list)
    for r in read_tsv(os.path.join(tables, "mge_local.tsv")):
        mge_by_ctg[(r["accession"], r["contig_id"])].append((int(r["start"]), int(r["end"])))
    contigs = {(r["accession"], r["contig_id"]): r for r in read_tsv(os.path.join(tables, "contigs.tsv"))}
    occ = {(r["accession"], r["gene_id"]): r
           for r in read_tsv(os.path.join(outdir, "occurrence_overlap.tsv"))}

    def union(iv):
        iv = sorted(iv)
        out = []
        for s, e in iv:
            if out and s <= out[-1][1]:
                out[-1][1] = max(out[-1][1], e)
            else:
                out.append([s, e])
        return out

    mge_union = {k: union(v) for k, v in mge_by_ctg.items()}

    # ---------------- stratified sampling of real instances ----------------
    rng = random.Random(20261015)
    picked = {}   # (acc,gid) -> set(edge types)
    by_acc_ctg = defaultdict(list)
    for r in valid:
        by_acc_ctg[(r["accession"], r["contig_id"])].append(r)

    def add(r, etype):
        picked.setdefault((r["accession"], r["gene_id"]), set()).add(etype)

    # precompute naive overlaps only where needed; use main table to find strata
    cand_1bp, cand_touch, cand_half_lo, cand_half_hi = [], [], [], []
    cand_minus, cand_plus, cand_circ = [], [], []
    for r in valid:
        o = occ[(r["accession"], r["gene_id"])]
        ov = int(o["mge_overlap_union_bp"])
        gs, ge = int(r["start"]), int(r["end"])
        glen = ge - gs
        U = mge_union.get((r["accession"], r["contig_id"]), [])
        if not U:
            continue
        touch = any(s == ge or e == gs for s, e in U)
        if ov == 1:
            cand_1bp.append(r)
        if ov == 0 and touch:
            cand_touch.append(r)
        if 0.45 <= ov / glen < 0.5:
            cand_half_lo.append(r)
        if 0.5 <= ov / glen <= 0.55:
            cand_half_hi.append(r)
        if ov > 0 and r["strand"] == "-":
            cand_minus.append(r)
        if ov > 0 and r["strand"] == "+":
            cand_plus.append(r)
        if contigs[(r["accession"], r["contig_id"])]["topology"] == "circular":
            cand_circ.append(r)
    for pool, tag, k in ((cand_1bp, "1bp_overlap", 3),
                         (cand_touch, "boundary_touch_no_overlap", 3),
                         (cand_half_lo, "half_threshold_just_below", 3),
                         (cand_half_hi, "half_threshold_just_above", 3),
                         (cand_minus, "minus_strand", 3),
                         (cand_plus, "plus_strand", 3),
                         (cand_circ, "circular_contig_gene", 3)):
        for r in rng.sample(pool, min(k, len(pool))):
            add(r, tag)
    # same annotation contig name across accessions: use contig map to find
    # two accessions whose annotation 'contig_1' both carry local MGE
    ct1 = [key for key in mge_union if key[1] and
           contigs[key]["annotation_contig_id"] == "contig_1"]
    if len(ct1) >= 2:
        for key in ct1[:2]:
            gs_list = by_acc_ctg.get(key, [])
            if gs_list:
                add(rng.choice(gs_list), "same_contig_name_cross_accession")
    # random fillers to reach >=24
    pool = [r for r in valid if (r["accession"], r["gene_id"]) not in picked]
    for r in rng.sample(pool, max(0, 24 - len(picked))):
        add(r, "random_filler")

    rows = []
    n_ok = 0
    for i, ((acc, gid), etypes) in enumerate(sorted(picked.items()), 1):
        r = next(x for x in by_acc_ctg[(acc, occ[(acc, gid)]["contig_id"])] if x["gene_id"] == gid)
        o = occ[(acc, gid)]
        gs, ge = int(r["start"]), int(r["end"])
        glen = ge - gs
        U = mge_union.get((acc, r["contig_id"]), [])
        nov = naive_overlap_bp(gs, ge, U)
        nany = nov > 0
        nhalf = 2 * nov >= glen
        mov = int(o["mge_overlap_union_bp"])
        many = o["mge_any_bp"] == "TRUE"
        mhalf = o["mge_half_gene"] == "TRUE"
        match = (nov, nany, nhalf) == (mov, many, mhalf)
        n_ok += match
        rows.append({"case_id": f"R{i:03d}", "case_kind": "real",
                     "edge_types": ";".join(sorted(etypes)), "accession": acc,
                     "feature_id": gid, "contig_id": r["contig_id"],
                     "gene_start": gs, "gene_end": ge, "strand": r["strand"],
                     "main_overlap_bp": mov, "naive_overlap_bp": nov,
                     "main_any_bp": many, "naive_any_bp": nany,
                     "main_half_gene": mhalf, "naive_half_gene": nhalf,
                     "match": "TRUE" if match else "FALSE",
                     "notes": f"gene_len={glen};n_mge_intervals_contig={len(U)}"})

    # ---------------- synthetic unit cases (not biological data) ----------------
    syn = []
    # S1 boundary touch
    syn.append(("S01", "boundary_touch_no_overlap", 100, 200, [(200, 350)]))
    # S2 1-bp overlap
    syn.append(("S02", "1bp_overlap", 100, 200, [(199, 350)]))
    # S3 exact 50% threshold (gene len 200, union overlap exactly 100)
    syn.append(("S03", "half_threshold_exact", 100, 300, [(150, 250)]))
    # S4 just below 50%
    syn.append(("S04", "half_threshold_just_below", 100, 300, [(151, 250)]))
    # S5 overlapping MGE intervals must be unioned (no double count)
    syn.append(("S05", "mge_union_double_count_guard", 100, 300, [(120, 220), (180, 280)]))
    syn_results = []
    for cid, tag, gs, ge, ivs in syn:
        U = union(ivs)
        nov = naive_overlap_bp(gs, ge, U)
        # replicate the main formula independently (sorted-merge overlap)
        mov = sum(max(0, min(ge, e) - max(gs, s)) for s, e in U)
        ok = nov == mov
        syn_results.append({"case_id": cid, "case_kind": "synthetic",
                            "edge_types": tag, "accession": "SYNTHETIC",
                            "feature_id": cid, "contig_id": "SYNTHETIC",
                            "gene_start": gs, "gene_end": ge,
                            "strand": "NA",
                            "main_overlap_bp": mov, "naive_overlap_bp": nov,
                            "main_any_bp": mov > 0, "naive_any_bp": nov > 0,
                            "main_half_gene": 2 * mov >= ge - gs,
                            "naive_half_gene": 2 * nov >= ge - gs,
                            "match": "TRUE" if ok else "FALSE",
                            "notes": f"synthetic unit case; intervals={ivs}"})
    # S6/S7 circular wrap: L=1000 circular, U=[(100,200)], offsets that wrap
    for cid, d, gs, ge, tag in (("S06", 850, 960, 990, "circular_wrap_offset"),
                                ("S07", 850, 0, 40, "circular_wrap_origin_segment")):
        L = 1000
        U = [(100, 200)]
        shifted = shift_naive(U, d, L, True)
        nov = naive_overlap_bp(gs, ge, shifted)
        # naive position circle, fully independent
        circle = {(p + d) % L for s, e in U for p in range(s, e)}
        nov2 = sum(1 for p in range(gs, ge) if p in circle)
        ok = nov == nov2 == 30 if cid == "S06" else nov == nov2 == 40
        syn_results.append({"case_id": cid, "case_kind": "synthetic",
                            "edge_types": tag, "accession": "SYNTHETIC",
                            "feature_id": cid, "contig_id": "SYNTHETIC_CIRC_L1000",
                            "gene_start": gs, "gene_end": ge, "strand": "NA",
                            "main_overlap_bp": nov, "naive_overlap_bp": nov2,
                            "main_any_bp": nov > 0, "naive_any_bp": nov2 > 0,
                            "main_half_gene": "NA", "naive_half_gene": "NA",
                            "match": "TRUE" if ok else "FALSE",
                            "notes": f"circular L=1000 offset={d} shifted={shifted}"})
    rows.extend(syn_results)

    # ---------------- degenerate-layout check ----------------
    deg_rows = []
    for r in read_tsv(os.path.join(outdir, "null_degenerate_contigs.tsv")):
        acc, ctg = r["accession"], r["contig_id"]
        L = int(contigs[(acc, ctg)]["length_bp"])
        U = mge_union[(acc, ctg)]
        s_min = min(s for s, _ in U)
        e_max = max(e for _, e in U)
        n_off = L - (e_max - s_min) + 1
        # genes on this contig must all be ineligible
        genes_here = [g for g in by_acc_ctg.get((acc, ctg), [])]
        all_inelig = all(occ[(g["accession"], g["gene_id"])]["eligible_for_null"] == "FALSE"
                         for g in genes_here)
        ok = n_off == 1 and int(r["n_legal_offsets"]) == 1 and all_inelig
        deg_rows.append({"case_id": f"DEG_{acc}_{ctg}", "case_kind": "degenerate_check",
                         "edge_types": "single_legal_offset_exclusion",
                         "accession": acc, "feature_id": "NA", "contig_id": ctg,
                         "gene_start": "NA", "gene_end": "NA", "strand": "NA",
                         "main_overlap_bp": "NA", "naive_overlap_bp": "NA",
                         "main_any_bp": "NA", "naive_any_bp": "NA",
                         "main_half_gene": "NA", "naive_half_gene": "NA",
                         "match": "TRUE" if ok else "FALSE",
                         "notes": f"L={L};span=[{s_min},{e_max});n_legal_offsets={n_off};"
                                  f"n_genes_on_contig={len(genes_here)};all_ineligible={all_inelig}"})
    rows.extend(deg_rows)

    # ---------------- synthetic circular genome through the REAL worker -------
    # exercises the numpy circular-shift code path of m8_overlap_null directly
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import m8_overlap_null as mainimpl
    syn_contig = "SYNTH_C1"
    syn_L = 1000
    syn_U = ((100, 200), (500, 510))
    syn_genes = [{"gene_id": "SYNTH_G1", "cluster": "synth_clA", "contig_id": syn_contig,
                  "start": 150, "end": 250, "partition": "Cloud", "strand": "+"},
                 {"gene_id": "SYNTH_G2", "cluster": "synth_clB", "contig_id": syn_contig,
                  "start": 980, "end": 1000, "partition": "Cloud", "strand": "-"},
                 {"gene_id": "SYNTH_G3", "cluster": "synth_clC", "contig_id": syn_contig,
                  "start": 0, "end": 60, "partition": "Cloud", "strand": "+"},
                 {"gene_id": "SYNTH_G4", "cluster": "synth_clD", "contig_id": syn_contig,
                  "start": 700, "end": 800, "partition": "Cloud", "strand": "-"}]
    ci = {syn_contig: (syn_L, "circular", syn_U, 0, syn_L - 1, "OK")}
    _, syn_clusters, _ = mainimpl.process_accession(("SYNTH_CIRC", ci, syn_genes))
    d_vec = mainimpl.draw_offsets("SYNTH_CIRC", syn_contig, "circular", 0, syn_L - 1, syn_L)
    circ_ok = True
    mism = 0
    for it in range(N_ITER):
        d = int(d_vec[it])
        circle = set()
        for s, e in syn_U:
            circle.update((p + d) % syn_L for p in range(s, e))
        for g in syn_genes:
            gs, ge = g["start"], g["end"]
            nov = sum(1 for p in range(gs, ge) if p in circle)
            ent = syn_clusters[g["cluster"]]
            if (nov > 0) != bool(ent["any"][it]):
                circ_ok = False; mism += 1
            if (2 * nov >= ge - gs) != bool(ent["half"][it]):
                circ_ok = False; mism += 1
    rows.append({"case_id": "SYNTH_CIRC_WORKER", "case_kind": "synthetic",
                 "edge_types": "circular_wrap_all_1000_iterations_through_main_worker",
                 "accession": "SYNTH_CIRC", "feature_id": "4 genes x 1000 iters x 2 criteria",
                 "contig_id": syn_contig, "gene_start": "NA", "gene_end": "NA",
                 "strand": "NA", "main_overlap_bp": "NA", "naive_overlap_bp": "NA",
                 "main_any_bp": "NA", "naive_any_bp": "NA",
                 "main_half_gene": "NA", "naive_half_gene": "NA",
                 "match": "TRUE" if circ_ok else "FALSE",
                 "notes": f"main worker vs naive position circle; mismatches={mism} of 8000 checks"})

    # ---------------- null end-to-end recheck ----------------
    def _null_recheck(grp, part, crit, iters):
        scope_accs = [a for a in meta if grp == "pooled" or group_of[a] == grp]
        # degenerate layouts in scope (independent recomputation)
        deg_ctgs = set()
        for acc in scope_accs:
            for (a2, ctg), U in mge_union.items():
                if a2 != acc:
                    continue
                c = contigs[(acc, ctg)]
                L = int(c["length_bp"])
                if c["topology"] == "circular":
                    n_off = L
                else:
                    s_min = min(s for s, _ in U)
                    e_max = max(e for _, e in U)
                    n_off = L - (e_max - s_min) + 1
                if n_off <= 1:
                    deg_ctgs.add((acc, ctg))
        elig_cls = set()
        for g in valid:
            if g["partition"] != part:
                continue
            if grp != "pooled" and group_of[g["accession"]] != grp:
                continue
            if (g["accession"], g["contig_id"]) not in deg_ctgs:
                elig_cls.add(g["cluster"])
        for it in iters:
            isd = iter_seed(it)
            cl_hit = defaultdict(bool)
            for acc in scope_accs:
                for (a2, ctg), U in mge_union.items():
                    if a2 != acc or (acc, ctg) in deg_ctgs:
                        continue
                    c = contigs[(acc, ctg)]
                    L = int(c["length_bp"])
                    circ = c["topology"] == "circular"
                    if circ:
                        dmin, dmax = 0, L - 1
                    else:
                        s_min = min(s for s, _ in U)
                        e_max = max(e for _, e in U)
                        dmin, dmax = -s_min, L - e_max
                    d = contig_offset(isd, acc, ctg, "circular" if circ else "linear",
                                      dmin, dmax, L)
                    shifted = shift_naive(U, d, L, circ)
                    for g in by_acc_ctg.get((acc, ctg), []):
                        if g["partition"] != part:
                            continue
                        gs, ge = int(g["start"]), int(g["end"])
                        tot = 0
                        for s, e in shifted:
                            ov = min(ge, e) - max(gs, s)
                            if ov > 0:
                                tot += ov
                        hit = tot > 0 if crit == "any_bp" else 2 * tot >= ge - gs
                        if hit:
                            cl_hit[g["cluster"]] = True
            n_hits = sum(1 for c in elig_cls if cl_hit.get(c, False))
            pr = next((r for r in perm_all
                       if r["group"] == grp and r["partition"] == part
                       and r["criterion"] == crit and r["iteration"] == str(it)), None)
            ok = pr is not None and int(pr["hits"]) == n_hits \
                and int(pr["n_eligible"]) == len(elig_cls) \
                and int(pr["seed"]) == iter_seed(it)
            rows.append({"case_id": f"NULL_{grp}_{part}_{crit}_iter{it}",
                         "case_kind": "null_iteration_recheck",
                         "edge_types": f"end_to_end_seed_to_hits;combo={grp}/{part}/{crit}",
                         "accession": "MULTI", "feature_id": "NA", "contig_id": "MULTI",
                         "gene_start": "NA", "gene_end": "NA", "strand": "NA",
                         "main_overlap_bp": "NA", "naive_overlap_bp": "NA",
                         "main_any_bp": "NA", "naive_any_bp": "NA",
                         "main_half_gene": "NA", "naive_half_gene": "NA",
                         "match": "TRUE" if ok else "FALSE",
                         "notes": f"iteration={it};seed={iter_seed(it)};"
                                  f"main_hits={pr['hits'] if pr else 'NA'};naive_hits={n_hits};"
                                  f"n_eligible={len(elig_cls)}"})

    nullres = [r for r in read_tsv(os.path.join(outdir, "interval_null_results.tsv"))
               if r["status"] == "OK"]
    nullres.sort(key=lambda r: int(r["n_eligible"]))
    combo = nullres[0]
    perm_all = read_tsv(os.path.join(outdir, "null_permutation_counts.tsv"))
    _null_recheck(combo["group"], combo["partition"], combo["criterion"], (1, 500, 1000))
    _null_recheck("pooled", "Strict_core", "any_bp", (7,))  # scope covers all 9 circular contigs

    # ---------------- write ----------------
    with open(os.path.join(outdir, "coordinate_audit.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=OUT_FIELDS)
        w.writeheader()
        w.writerows(rows)
    n_real = sum(1 for r in rows if r["case_kind"] == "real")
    n_match = sum(1 for r in rows if r["match"] == "TRUE")
    summary = {
        "n_real_instances": n_real,
        "n_synthetic_cases": sum(1 for r in rows if r["case_kind"] == "synthetic"),
        "n_degenerate_checks": sum(1 for r in rows if r["case_kind"] == "degenerate_check"),
        "n_null_iteration_rechecks": sum(1 for r in rows if r["case_kind"] == "null_iteration_recheck"),
        "n_total_cases": len(rows),
        "n_match": n_match,
        "n_mismatch": len(rows) - n_match,
        "real_pool_sizes": {"1bp": len(cand_1bp), "touch": len(cand_touch),
                            "half_lo": len(cand_half_lo), "half_hi": len(cand_half_hi),
                            "minus": len(cand_minus), "plus": len(cand_plus),
                            "circular": len(cand_circ)},
        "edge_types_covered_by_real": sorted({t for r in rows if r["case_kind"] == "real"
                                              for t in r["edge_types"].split(";")}),
        "smallest_combo_rechecked": f"{combo['group']}/{combo['partition']}/{combo['criterion']}",
        "conclusion": "CONSISTENT" if n_match == len(rows) else "MISMATCH_FOUND",
    }
    with open(os.path.join(outdir, "coordinate_audit_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2, ensure_ascii=False)
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
