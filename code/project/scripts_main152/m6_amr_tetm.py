#!/usr/bin/env python3
"""M6 (main152, 2026-10-15): AMR / tet(M) per-genome status, per-species
carriage, locus evidence, and reconciliation against the display template.

Reads (all read-only):
  - 152 roster / metadata / display metadata from the frozen input packet
  - per-tool coverage index -> historical RGI raw outputs (trash-resolved)
  - phaseB_v7 tetM summary tables (hash-bound reuse, re-subset to 152)
  - phaseB_v7 tetM reference blast outputs (locus coordinate verification)
  - phaseB_v7 08_gi_mge/tables/mge.tsv (0-based half-open, annotation contig
    namespace) + 01_annotation/contig_id_maps for MGE same-contig / overlap
  - 00_preflight/resolved_paths.tsv for frozen FASTA paths (contig lengths)

tet(M) four-state rules (identical to scripts/revision_20261001/d1_parse_rgi.py):
  positive_complete : Perfect/Strict hit AND ref-length coverage >= 95
  indeterminate     : tet(M) hits exist but none meet the above
  negative          : RGI output present/parses, zero tet(M) hits
  missing           : RGI output absent/empty/unparsable

Location grading (2026-10-15 rules; the historical ">1 Mb => chromosomal_verified"
rule is NOT reused):
  predicted_plasmid_associated            : MOB-suite molecule_type==plasmid OR
                                            geNomad plasmid call on the contig
  chromosomal_verified                    : Complete Genome assembly with a single
                                            contig whose length == genome size
                                            (confirmed chromosome replicon) AND
                                            MOB-suite molecule_type==chromosome
                                            AND no plasmid evidence AND locus
                                            mapped inside that contig
  chromosome_associated_assembly_supported: draft contig > 1 Mb, no plasmid
                                            evidence (large draft contigs are NOT
                                            called chromosomal_verified)
  undetermined_contig_boundary            : locus within 2 kb of a contig end on
                                            a short contig
  undetermined_small_contig               : otherwise

Statistical unit = genome (not hit row / not copy).
Outputs land in results_revision_20261015_main152_species/07_amr_tetm/.
"""
import csv
import json
import os
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common

OUT = os.path.join(common.RESULTS, "07_amr_tetm")
PHASEB = os.path.join(common.ROOT, "results_revision_20261010_phaseB_v7")
AMR_V7 = os.path.join(PHASEB, "06_annotation_summaries/amr_tetm")
HIST_BLAST = os.path.join(common.SOURCE_ROOT,
                          "results_revision_20261001/07_amr_tetm/blast")
NEW_BLAST = os.path.join(AMR_V7, "blast_new19")
MGE_TSV = os.path.join(PHASEB, "08_gi_mge/tables/mge.tsv")
CONTIG_MAPS = os.path.join(common.RESULTS, "01_annotation/contig_id_maps")

LOCUS_MERGE_GAP = 5000
BOUNDARY_SLACK = 2000
TETM_FAMILY = "tetracycline-resistant ribosomal protection protein"
RGI_VERSION = ("RGI 6.0.3 (biocontainers docker image rgi:6.0.3--pyha8f3691_0), "
               "CARD data_version 4.0.2, `-t contig -a DIAMOND` "
               "(per results_revision_20261010_phaseB_v7/06_annotation_summaries/"
               "per_genome_tool_versions.tsv)")

HIST_FILES = [
    os.path.join(AMR_V7, "tetm_genome_classification_v7.tsv"),
    os.path.join(AMR_V7, "tetm_locus_classification_v7.tsv"),
    os.path.join(AMR_V7, "locus_evidence_grading_v7.tsv"),
    os.path.join(AMR_V7, "amr_hits_per_genome_v7.tsv"),
    os.path.join(AMR_V7, "amr_per_genome_status_v7.tsv"),
]


def parse_rgi(acc, txt_path):
    """Replicates d1_parse_rgi.py parse() plus AMR family summary."""
    if txt_path is None or not os.path.isfile(txt_path) or os.path.getsize(txt_path) == 0:
        return dict(run_status="missing", tetM_status="missing",
                    tetM_best_cutoff="NA", tetM_best_identity="NA",
                    tetM_best_cov="NA", n_tetM_hits=0, detail="rgi_output_absent",
                    fam_ps=Counter(), fam_loose=Counter(), n_ps=0, n_loose=0,
                    n_hits_all=0)
    try:
        with open(txt_path, newline="", encoding="utf-8") as fh:
            rd = csv.DictReader(fh, delimiter="\t")
            if rd.fieldnames is None:
                raise ValueError("no header")
            rows = list(rd)
    except Exception as e:
        return dict(run_status="missing", tetM_status="missing",
                    tetM_best_cutoff="NA", tetM_best_identity="NA",
                    tetM_best_cov="NA", n_tetM_hits=0,
                    detail=f"rgi_output_unparsable:{e}",
                    fam_ps=Counter(), fam_loose=Counter(), n_ps=0, n_loose=0,
                    n_hits_all=0)
    tet = [r for r in rows
           if "tet(M)" in (r.get("AMR Gene Family") or "")
           or "tet(M)" in (r.get("Best_Hit_ARO") or "")]
    fam_ps, fam_loose = Counter(), Counter()
    n_ps = n_loose = 0
    for r in rows:
        fam = (r.get("AMR Gene Family") or "NA").strip() or "NA"
        if (r.get("Cut_Off") or "") in ("Perfect", "Strict"):
            fam_ps[fam] += 1
            n_ps += 1
        else:
            fam_loose[fam] += 1
            n_loose += 1
    base = dict(run_status="ok", fam_ps=fam_ps, fam_loose=fam_loose,
                n_ps=n_ps, n_loose=n_loose, n_hits_all=len(rows))
    if not tet:
        return dict(base, tetM_status="negative", tetM_best_cutoff="NA",
                    tetM_best_identity="NA", tetM_best_cov="NA",
                    n_tetM_hits=0, detail="no_tetM_hit")
    rank = {"Perfect": 3, "Strict": 2, "Loose": 1}

    def key(r):
        return (rank.get(r.get("Cut_Off", ""), 0),
                float(r.get("Percentage Length of Reference Sequence") or 0),
                float(r.get("Best_Identities") or 0))
    tet.sort(key=key, reverse=True)
    best = tet[0]
    complete = [r for r in tet if r.get("Cut_Off") in ("Perfect", "Strict")
                and float(r.get("Percentage Length of Reference Sequence") or 0) >= 95]
    detail = ";".join(f"{r.get('Cut_Off')}:{r.get('Best_Identities')}:"
                      f"{r.get('Percentage Length of Reference Sequence')}"
                      for r in tet[:4])
    return dict(base,
                tetM_status="positive_complete" if complete else "indeterminate",
                tetM_best_cutoff=best.get("Cut_Off", ""),
                tetM_best_identity=best.get("Best_Identities", ""),
                tetM_best_cov=best.get("Percentage Length of Reference Sequence", ""),
                n_tetM_hits=len(tet), detail=detail)


def contig_lengths(fa):
    lens, name, n = {}, None, 0
    with open(fa) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    lens[name] = n
                name, n = line[1:].split()[0], 0
            else:
                n += len(line.strip())
    if name is not None:
        lens[name] = n
    return lens


def merge_hsps(hsps):
    """hsps: (qstart,qend,pident) -> (union_cov_bp, weighted_pid)."""
    hsps = sorted(hsps, key=lambda h: min(h[0], h[1]))
    merged, cov, ident_sum = [], 0, 0.0
    for qs, qe, pid in hsps:
        a, b = min(qs, qe), max(qs, qe)
        if merged and a <= merged[-1][1]:
            na = max(a, merged[-1][1] + 1)
            if na <= b:
                cov += b - na + 1
                ident_sum += pid * (b - na + 1)
            merged[-1][1] = max(merged[-1][1], b)
        else:
            merged.append([a, b])
            cov += b - a + 1
            ident_sum += pid * (b - a + 1)
    return cov, (ident_sum / cov if cov else 0.0)


def classify(covA, pidA, covD, pidD):
    if max(covA, covD) < 80:
        return "partial"
    if pidA > pidD + 0.5:
        return "American-like"
    if pidD > pidA + 0.5:
        return "Dutch-like"
    return "mixed"


def recompute_loci(acc, fa_lens):
    """Re-run p5a_03 locus grouping on the raw blast file. Returns list of
    dicts (contig, start1, end1, strand, locus_class, n_hsp, cov/pid stats)
    or None if no blast file."""
    f_hist = os.path.join(HIST_BLAST, f"{acc}.tsv")
    f_new = os.path.join(NEW_BLAST, f"{acc}.tsv")
    if os.path.isfile(f_hist):
        bf, origin = f_hist, "historical_v1_verified"
    elif os.path.isfile(f_new):
        bf, origin = f_new, "this_round_same_method"
    else:
        return None, "missing"
    hsps, qlen = [], {}
    with open(bf) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 14:
                continue
            qid, sid, pid = p[0], p[1], float(p[2])
            qs, qe, ss, se = int(p[6]), int(p[7]), int(p[8]), int(p[9])
            qlen[qid] = int(p[12])
            klass = "American" if qid.startswith("American") else "Dutch"
            hsps.append((sid, min(ss, se), max(ss, se), klass, qs, qe, pid,
                         "+" if ss <= se else "-"))
    by_contig = defaultdict(list)
    for h in hsps:
        by_contig[h[0]].append(h)
    loci = []
    for ctg, hl in by_contig.items():
        hl.sort(key=lambda h: h[1])
        cur_s, cur_e, cur = None, None, []
        for h in hl:
            if cur_s is None or h[1] > cur_e + LOCUS_MERGE_GAP:
                if cur:
                    loci.append((ctg, cur_s, cur_e, cur))
                cur_s, cur_e, cur = h[1], h[2], [h]
            else:
                cur_e = max(cur_e, h[2])
                cur.append(h)
        if cur:
            loci.append((ctg, cur_s, cur_e, cur))
    out = []
    for ctg, ls, le, lh in sorted(loci):
        per_ref = defaultdict(list)
        strands = set()
        for h in lh:
            per_ref[h[3]].append((h[4], h[5], h[6]))
            strands.add(h[7])
        st = {}
        for k, ref_h in per_ref.items():
            L = qlen[[q for q in qlen if q.startswith(k)][0]]
            c, w = merge_hsps(ref_h)
            st[k] = (c / L * 100, w)
        covA, pidA = st.get("American", (0.0, 0.0))
        covD, pidD = st.get("Dutch", (0.0, 0.0))
        out.append(dict(contig=ctg, start1=ls, end1=le,
                        strand="+" if strands == {"+"} else
                        ("-" if strands == {"-"} else "mixed"),
                        locus_class=classify(covA, pidA, covD, pidD),
                        n_hsp=len(lh), covA=covA, pidA=pidA, covD=covD, pidD=pidD,
                        contig_len=fa_lens.get(ctg, 0)))
    return out, origin


def wilson(k, n, z=1.959963984540054):
    if n == 0:
        return None, None
    p = k / n
    denom = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / denom
    half = z * ((p * (1 - p) / n + z * z / (4 * n * n)) ** 0.5) / denom
    return centre - half, centre + half


def main():
    os.makedirs(OUT, exist_ok=True)
    roster = common.read_roster()
    meta = {r["strain"]: r for r in common.load_metadata()}
    group = {a: meta[a]["group_for_testing"] for a in roster}
    display = {r["accession"]: r for r in common.read_tsv(
        os.path.join(common.PACKET, "00_input/cohort/main152_display_metadata.tsv"))}
    coverage = {r["accession"]: r for r in common.read_tsv(os.path.join(
        common.PACKET,
        "00_input/reference_sources/per_tool_coverage_historical_main152.tsv"))}
    fa_of = {r["accession"]: r["resolved_path"] for r in common.read_tsv(
        os.path.join(common.RESULTS, "00_preflight/resolved_paths.tsv"))
        if r["category"] == "fasta"}

    hist_hashes = {os.path.basename(f): common.sha256_file(f) for f in HIST_FILES}
    v7_status = {r["accession"]: r for r in common.read_tsv(
        os.path.join(AMR_V7, "amr_per_genome_status_v7.tsv"))}
    v7_loci = [r for r in common.read_tsv(
        os.path.join(AMR_V7, "locus_evidence_grading_v7.tsv"))
        if r["accession"] in meta]
    v7_genome = {r["accession"]: r for r in common.read_tsv(
        os.path.join(AMR_V7, "tetm_genome_classification_v7.tsv"))}

    # ---- step 1: parse RGI raw outputs ------------------------------------
    per_genome = {}
    rgi_hash_rows = []
    for acc in roster:
        ev_dir, note = common.resolve_path(coverage[acc]["rgi_evidence"])
        txt = os.path.join(ev_dir, f"{acc}.txt") if ev_dir else None
        js = os.path.join(ev_dir, f"{acc}.json") if ev_dir else None
        rec = parse_rgi(acc, txt)
        if rec["run_status"] == "ok" and (not js or not os.path.isfile(js)
                                          or os.path.getsize(js) == 0):
            rec["run_status"] = "ok_txt_only"
            rec["detail"] += "|json_missing"
        rec["rgi_dir"] = ev_dir or "NA"
        rec["rgi_resolve"] = note
        per_genome[acc] = rec
        rgi_hash_rows.append({
            "accession": acc,
            "rgi_txt_path": txt or "NA",
            "rgi_txt_sha256": common.sha256_file(txt) if txt and os.path.isfile(txt) else "NA",
            "rgi_json_sha256": common.sha256_file(js) if js and os.path.isfile(js) else "NA",
            "resolve_note": note,
            "rgi_final_index": coverage[acc]["rgi_final"],
        })

    # cross-check vs v7 per-genome status
    status_mismatch = []
    for acc in roster:
        v7 = v7_status.get(acc)
        if v7 is None:
            status_mismatch.append((acc, "absent_in_v7_table"))
            continue
        mine = per_genome[acc]
        if mine["tetM_status"] != v7["tetM_status"]:
            status_mismatch.append((acc, f"mine={mine['tetM_status']}/v7={v7['tetM_status']}"))
        elif mine["tetM_status"] == "positive_complete" and (
                str(mine["tetM_best_cov"]) != v7["tetM_best_coverage_pct"]):
            status_mismatch.append((acc, "coverage_mismatch"))

    # ---- step 2: recompute loci from raw blast, verify vs v7 --------------
    locus_verify_rows = []
    my_loci = {}  # acc -> list
    for acc in roster:
        lens = contig_lengths(fa_of[acc])
        loci, origin = recompute_loci(acc, lens)
        my_loci[acc] = loci or []
        v7l = sorted([r for r in v7_loci if r["accession"] == acc],
                     key=lambda r: (r["contig"], int(r["start"])))
        mine = sorted(my_loci[acc], key=lambda d: (d["contig"], d["start1"]))
        match = (len(mine) == len(v7l) and all(
            m["contig"] == v["contig"] and m["start1"] == int(v["start"])
            and m["end1"] == int(v["end"]) and m["locus_class"] == v["locus_class"]
            and m["contig_len"] == int(v["contig_len"])
            for m, v in zip(mine, v7l)))
        locus_verify_rows.append({
            "accession": acc, "blast_origin": origin,
            "n_loci_recomputed": len(mine), "n_loci_v7": len(v7l),
            "exact_match": "yes" if match else "NO",
        })
        if not match and (mine or v7l):
            for m, v in zip(mine, v7l):
                locus_verify_rows[-1]["detail"] = (
                    f"recomputed={m['contig']}:{m['start1']}-{m['end1']}:"
                    f"{m['locus_class']} vs v7={v['contig']}:{v['start']}-"
                    f"{v['end']}:{v['locus_class']}")

    # tet(M) positive genomes must have >=1 verified locus
    pos_no_locus = [a for a in roster
                    if per_genome[a]["tetM_status"] == "positive_complete"
                    and not my_loci[a]]

    # ---- step 3: MGE same-contig / interval overlap -----------------------
    mge_by_gc = defaultdict(list)  # (genome, annotation_contig) -> rows
    for r in common.read_tsv(MGE_TSV):
        mge_by_gc[(r["genome"], r["contig"])].append(r)
    contig_map = {}
    for acc in roster:
        contig_map[acc] = {r["frozen_contig_id"]: r["annotation_contig_id"]
                           for r in common.read_tsv(
                               os.path.join(CONTIG_MAPS, f"{acc}.contig_map.tsv"))}

    # ---- step 4: per-locus evidence table ---------------------------------
    locus_rows = []
    for acc in roster:
        v7l = sorted([r for r in v7_loci if r["accession"] == acc],
                     key=lambda r: (r["contig"], int(r["start"])))
        mine = sorted(my_loci[acc], key=lambda d: (d["contig"], d["start1"]))
        m = meta[acc]
        asmb_complete = (m["assembly_level"] == "Complete Genome"
                         and int(m["n_contigs"]) == 1)
        for i, (v, mm) in enumerate(zip(v7l, mine), 1):
            clen = mm["contig_len"]
            start0, end0 = mm["start1"] - 1, mm["end1"]  # 1-based closed -> 0-based half-open
            mob_t, mob_c = v["mob_molecule_type"], v["mob_circularity"]
            gn = v["genomad_class"]
            is_plasmid = (mob_t == "plasmid") or (gn == "plasmid")
            near_end = (start0 < BOUNDARY_SLACK) or (clen - end0 < BOUNDARY_SLACK)
            single_chr = (asmb_complete and clen == int(m["genome_size"]))
            if is_plasmid:
                grade = "predicted_plasmid_associated"
                evalu = "evaluable_predicted"
            elif single_chr and mob_t == "chromosome":
                grade = "chromosomal_verified"
                evalu = "evaluable"
            elif clen > 1_000_000:
                grade = "chromosome_associated_assembly_supported"
                evalu = "evaluable"
            elif near_end:
                grade = "undetermined_contig_boundary"
                evalu = "undetermined"
            else:
                grade = "undetermined_small_contig"
                evalu = "undetermined"
            ann_contig = contig_map[acc].get(mm["contig"])
            mges = mge_by_gc.get((acc, ann_contig), []) if ann_contig else []
            local = [g for g in mges if g["whole_contig"] == "0"]
            whole = [g for g in mges if g["whole_contig"] != "0"]
            overlap = [g["mge_id"] for g in local
                       if int(g["start"]) < end0 and start0 < int(g["end"])]
            replicon_ev = (f"assembly_level={m['assembly_level']};n_contigs={m['n_contigs']};"
                           f"contig_len={clen};mob_molecule_type={mob_t or 'NA'};"
                           f"mob_circularity={mob_c or 'NA'};genomad_class={gn or 'NA'};"
                           f"whole_contig_mge_calls={len(whole)}")
            locus_rows.append({
                "accession": acc, "locus_id": v["locus_id"],
                "contig_id": mm["contig"],
                "coordinate_system": "0-based_half-open_frozen_FASTA_contig",
                "start": start0, "end": end0, "strand": mm["strand"],
                "location_grade": grade, "location_evaluability": evalu,
                "replicon_evidence": replicon_ev,
                "MGE_same_contig": "yes" if local else "no",
                "MGE_interval_overlap": "yes" if overlap else "no",
                "MGE_overlap_ids": ";".join(overlap) if overlap else "NA",
                "contig_end_distance": min(start0, clen - end0),
                "contig_len": clen,
                "locus_class": v["locus_class"],
                "cov_american_pct": v["cov_american_pct"],
                "pid_american": v["pid_american"],
                "cov_dutch_pct": v["cov_dutch_pct"],
                "pid_dutch": v["pid_dutch"],
                "n_hsp": v["n_hsp"],
                "identical_interval_copies_in_genome": v["identical_interval_copies_in_genome"],
                "evidence_origin": v["evidence_origin"],
                "status": "OK",
            })

    # ---- step 5: per-genome output table ----------------------------------
    pg_rows = []
    for acc in roster:
        rec = per_genome[acc]
        fam_ps = rec["fam_ps"]
        loose_only = sorted(f for f in rec["fam_loose"] if f not in fam_ps)
        pg_rows.append({
            "accession": acc, "species_group": group[acc],
            "RGI_run_status": rec["run_status"],
            "tetM_status": rec["tetM_status"],
            "tetM_best_cutoff": rec["tetM_best_cutoff"],
            "tetM_reference_coverage": rec["tetM_best_cov"],
            "n_tetM_loci": len(my_loci[acc]),
            "tetM_best_identity_pct": rec["tetM_best_identity"],
            "n_tetM_hits": rec["n_tetM_hits"],
            "tetM_hit_detail": rec["detail"],
            "n_amr_hits_perfect_strict": rec["n_ps"],
            "n_amr_families_perfect_strict": len(fam_ps),
            "amr_family_hitcounts_perfect_strict": ";".join(
                f"{f}:{fam_ps[f]}" for f in sorted(fam_ps)) or "none",
            "amr_families_loose_only": ";".join(loose_only) or "none",
            "n_amr_hits_loose": rec["n_loose"],
            "n_amr_hits_all": rec["n_hits_all"],
            "genome_category_v7": v7_genome.get(acc, {}).get("genome_category", "NA"),
            "rgi_txt_sha256": next(r["rgi_txt_sha256"] for r in rgi_hash_rows
                                   if r["accession"] == acc),
            "rgi_evidence_dir": rec["rgi_dir"],
            "rgi_path_resolve": rec["rgi_resolve"],
            "status": "OK" if rec["run_status"].startswith("ok") else "BLOCKED",
        })

    # ---- step 6: carriage by species --------------------------------------
    fam_ps_carrier = {acc: set(per_genome[acc]["fam_ps"]) for acc in roster}
    fam_loose_carrier = {acc: set(per_genome[acc]["fam_loose"]) for acc in roster}
    all_fams = sorted({f for acc in roster for f in fam_ps_carrier[acc]}
                      | {f for acc in roster for f in fam_loose_carrier[acc]})
    features = [("gene:tet(M)", "tet(M) gene, Perfect/Strict with reference "
                              "length coverage >=95% (four-state rule)")]
    features += [(f"family:{f}", f"AMR Gene Family P/S carrier: {f}")
                 for f in all_fams]
    car_rows = []
    for feat, desc in features:
        for g in common.GROUPS:
            accs = [a for a in roster if group[a] == g]
            n_total = len(accs)
            n_pos = n_neg = n_ind = n_miss = 0
            for a in accs:
                rec = per_genome[a]
                if not rec["run_status"].startswith("ok"):
                    n_miss += 1
                    continue
                if feat == "gene:tet(M)":
                    s = rec["tetM_status"]
                    if s == "positive_complete":
                        n_pos += 1
                    elif s == "negative":
                        n_neg += 1
                    elif s == "indeterminate":
                        n_ind += 1
                    else:
                        n_miss += 1
                else:
                    fname = feat.split(":", 1)[1]
                    if fname in fam_ps_carrier[a]:
                        n_pos += 1
                    elif fname in fam_loose_carrier[a]:
                        n_ind += 1
                    else:
                        n_neg += 1
            n_eval = n_pos + n_neg
            lo, hi = wilson(n_pos, n_eval)
            car_rows.append({
                "feature": feat, "feature_definition": desc,
                "species_group": g, "n_total": n_total,
                "n_evaluable": n_eval, "n_positive": n_pos,
                "n_negative": n_neg, "n_indeterminate": n_ind,
                "n_missing": n_miss,
                "prevalence": f"{n_pos / n_eval:.6f}" if n_eval else "NA",
                "CI_low": f"{lo:.6f}" if lo is not None else "NA",
                "CI_high": f"{hi:.6f}" if hi is not None else "NA",
                "ci_method": "Wilson 95% (binomial denominator n_evaluable)",
                "status": "OK" if n_eval else "NOT_ESTIMABLE",
                "reason": "" if n_eval else "n_evaluable=0",
            })

    # ---- step 7: Fisher input (binary P/S carrier, carrier>=5) ------------
    fisher_feats = [("gene:tet(M)",
                     lambda a: per_genome[a]["tetM_status"] == "positive_complete")]
    fisher_feats += [(f"family:{f}",
                      (lambda f: lambda a: f in fam_ps_carrier[a])(f))
                     for f in all_fams]
    ledger_rows, fisher_rows = [], []
    for feat, fn in fisher_feats:
        carriers = sum(1 for a in roster if fn(a))
        testable = carriers >= 5
        ledger_rows.append({"feature": feat, "carriers_152": carriers,
                            "tested": "yes" if testable else "no",
                            "reason": "" if testable else "carrier<5_in_152"})
        for g in common.GROUPS:
            tgt = [a for a in roster if group[a] == g]
            rst = [a for a in roster if group[a] != g]
            a_ = sum(1 for x in tgt if fn(x))
            b_ = len(tgt) - a_
            c_ = sum(1 for x in rst if fn(x))
            d_ = len(rst) - c_
            if not testable:
                st, reason = "FILTERED", "carrier<5_in_152"
            elif a_ + b_ == 0 or c_ + d_ == 0 or (a_ == 0 and c_ == 0) \
                    or (b_ == 0 and d_ == 0):
                st, reason = "NOT_ESTIMABLE", "no_variation_or_empty_margin"
            else:
                st, reason = "VALID", ""
            fisher_rows.append({
                "family_id": "AMR_prevalence", "feature": feat,
                "target_group": g, "a": a_, "b": b_, "c": c_, "d": d_,
                "target_n": len(tgt), "rest_n": len(rst),
                "status": st, "reason": reason,
            })

    # ---- step 8: reconciliation vs display metadata -----------------------
    rec_rows = []
    n_conc = 0
    for acc in roster:
        mine = per_genome[acc]["tetM_status"]
        disp = display[acc]["tetM_state"]
        mine_bin = {"positive_complete": "positive", "negative": "negative"}.get(mine, mine)
        concord = "yes" if mine_bin == disp else "NO"
        if concord == "yes":
            n_conc += 1
        rec_rows.append({
            "accession": acc, "species_group": group[acc],
            "tetM_status_m6": mine, "tetM_state_display": disp,
            "concordant": concord,
            "reason": "" if concord == "yes" else
            f"m6={mine} vs display={disp}; see evidence",
            "evidence": (f"RGI {per_genome[acc]['detail']}; "
                         f"n_loci={len(my_loci[acc])}; "
                         f"rgi_txt_sha256={next(r['rgi_txt_sha256'] for r in rgi_hash_rows if r['accession'] == acc)[:16]}..."),
        })

    # ---- write outputs -----------------------------------------------------
    common.write_tsv(os.path.join(OUT, "amr_tetm_per_genome.tsv"),
                     list(pg_rows[0].keys()), pg_rows)
    common.write_tsv(os.path.join(OUT, "carriage_by_species.tsv"),
                     list(car_rows[0].keys()), car_rows)
    common.write_tsv(os.path.join(OUT, "tetm_reconciliation.tsv"),
                     list(rec_rows[0].keys()), rec_rows)
    common.write_tsv(os.path.join(OUT, "tetm_locus_evidence.tsv"),
                     list(locus_rows[0].keys()) if locus_rows else
                     ["accession", "locus_id", "contig_id", "coordinate_system",
                      "start", "end", "strand", "location_grade",
                      "location_evaluability", "replicon_evidence",
                      "MGE_same_contig", "MGE_interval_overlap",
                      "MGE_overlap_ids", "contig_end_distance", "contig_len",
                      "locus_class", "cov_american_pct", "pid_american",
                      "cov_dutch_pct", "pid_dutch", "n_hsp",
                      "identical_interval_copies_in_genome", "evidence_origin",
                      "status"], locus_rows)
    common.write_tsv(os.path.join(OUT, "m6_fisher_input.tsv"),
                     list(fisher_rows[0].keys()), fisher_rows)
    common.write_tsv(os.path.join(OUT, "amr_feature_filter_ledger.tsv"),
                     list(ledger_rows[0].keys()), ledger_rows)
    common.write_tsv(os.path.join(OUT, "rgi_input_hashes.tsv"),
                     list(rgi_hash_rows[0].keys()), rgi_hash_rows)
    common.write_tsv(os.path.join(OUT, "locus_blast_verification.tsv"),
                     list(locus_verify_rows[0].keys()), locus_verify_rows)
    with open(os.path.join(OUT, "historical_v7_table_hashes.tsv"), "w") as fh:
        fh.write("file\tsha256\n")
        for k, v in hist_hashes.items():
            fh.write(f"{k}\t{v}\n")

    summary = {
        "module": "M6_amr_tetm", "run_id": common.RUN_ID,
        "n_genomes": len(roster),
        "tetM_status_by_group": {
            g: dict(Counter(per_genome[a]["tetM_status"]
                            for a in roster if group[a] == g))
            for g in common.GROUPS},
        "rgi_run_status": dict(Counter(per_genome[a]["run_status"]
                                       for a in roster)),
        "v7_status_mismatches": status_mismatch,
        "positive_without_verified_locus": pos_no_locus,
        "locus_recompute_mismatches": [r["accession"] for r in locus_verify_rows
                                       if r["exact_match"] != "yes"],
        "n_loci_152": len(locus_rows),
        "location_grade_counts": dict(Counter(r["location_grade"]
                                              for r in locus_rows)),
        "display_concordant": n_conc,
        "display_total": len(roster),
        "amr_families_ps": {f: sum(1 for a in roster if f in fam_ps_carrier[a])
                           for f in all_fams},
    }
    with open(os.path.join(OUT, "m6_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2, ensure_ascii=False)
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
