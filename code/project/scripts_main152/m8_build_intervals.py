#!/usr/bin/env python3
"""M8 (main152): build unified interval tables from raw per-genome tool outputs.

Unified key: (accession, contig_id) where contig_id = frozen FASTA contig ID.
Unified coordinates: 0-based half-open [start, end). Original coordinates and
the conversion are kept in every output row.

Coordinate systems of sources (verified 2026-10-15):
  - occurrence_map.tsv (new PGAP2 dictionary): Bakta annotation contig names
    (contig_N); coordinates are 0-based half-open [start, end) from PGAP2
    total.involved_annot -- NOT 1-based closed as once assumed. Verified
    against Bakta GFF CDS lines for all 328,376 rows: occurrence start ==
    GFF_start - 1 and occurrence end == GFF_end, zero exceptions (the 55
    rows with start=0 all map to GFF start=1). The verification is re-run
    by this script on every execution (see build_qc.json).
  - ISEScan <acc>.fna.tsv / .fna.gff: frozen contig IDs, 1-based closed
    (tsv isBegin/isEnd == gff start/end, confirmed identical values).
  - IntegronFinder .integrons: frozen contig IDs, 1-based closed
    (pos_beg/pos_end). Zero hits across all historical genomes.
  - geNomad find_proviruses <acc>_provirus.tsv: source_seq = frozen contig ID,
    start/end 1-based closed (verified by exact subsequence match against the
    frozen FASTA for GCA_000173935.1 ACEN01000080.1|provirus_34_44232).
  - geNomad summary virus/plasmid: whole-contig predictions (coordinates==NA
    in virus_summary; all plasmid_summary rows). Provirus rows in
    virus_summary carry coordinates and duplicate find_proviruses output ->
    not used here (count cross-checked only).
  - MOB-suite contig_report.txt molecule_type==plasmid: whole-contig plasmid.
  - IslandPath-DIMOB islands/<acc>__<contig_N>.gff: Bakta annotation contig
    names (per-contig GBK derived from Bakta annotation; LOCUS length matches
    contig_id_maps length), 1-based closed GFF -> mapped to frozen IDs via
    01_annotation/contig_id_maps/<acc>.contig_map.tsv.

Circularity evidence: only explicit '[topology=circular]' in the Bakta
annotation <acc>.fna headers (per annotation contig name).

Outputs (results_revision_20261015_main152_species/09_interval_overlap/tables/):
  contigs.tsv, genes.tsv, mge_local.tsv, whole_contig.tsv, gi.tsv,
  tool_coverage_ledger.tsv, build_qc.json
"""
from pathlib import Path
import csv
import json
import os
import re
import sys
from collections import Counter, defaultdict

ROOT = str(Path(os.environ["NEISSERIA_PROJECT_ROOT"]))
RES = os.path.join(ROOT, "results_revision_20261015_main152_species")
OUT = os.path.join(RES, "09_interval_overlap", "tables")
SOURCE_ROOT = os.environ.get("NEISSERIA_SOURCE_ROOT", ROOT)

META = os.path.join(RES, "00_freeze/input_packet_snapshot/00_input/cohort/main152_metadata.tsv")
COV = os.path.join(RES, "00_freeze/input_packet_snapshot/00_input/reference_sources/per_tool_coverage_historical_main152.tsv")
ANN_MAN = os.path.join(RES, "01_annotation/annotation_manifest_main152.tsv")
CMAP_DIR = os.path.join(RES, "01_annotation/contig_id_maps")
OCC = os.path.join(RES, "02_pangenome/pgap2_main152/pav_tables/occurrence_map.tsv")
PART = os.path.join(RES, "02_pangenome/pgap2_main152/pav_tables/partition_table.tsv")
RERUN_IP = os.path.join(ROOT, "results_revision_20261010_phaseB_v7/06_annotation_summaries/per_genome_reruns/islandpath/islands")


def resolve(p):
    """Map a historical absolute path to its current (possibly trashed) location."""
    if not p:
        return None
    if os.path.exists(p):
        return p
    rel = p.split("Neisseria_subflava/")[-1]
    alt = os.path.join(SOURCE_ROOT, rel)
    if os.path.exists(alt):
        return alt
    return None


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def main():
    os.makedirs(OUT, exist_ok=True)
    qc = {"warnings": [], "per_tool": {}}

    meta = {r["strain"]: r for r in read_tsv(META)}
    accs = [r["strain"] for r in read_tsv(META)]
    assert len(accs) == 152 and len(set(accs)) == 152
    group_of = {a: meta[a]["group_for_testing"] for a in accs}

    cov = {r["accession"]: r for r in read_tsv(COV)}
    ann = {r["accession"]: r for r in read_tsv(ANN_MAN)}
    part = {}
    for r in read_tsv(PART):
        if r["check"] != "MATCH":
            qc["warnings"].append(f"partition mismatch {r['cluster']}")
        part[r["cluster"]] = r["partition_recomputed"]

    # ---- contig maps + topology ----
    contigs = {}          # (acc, frozen) -> dict
    ann2frozen = {}       # acc -> {ann_id: frozen}
    for acc in accs:
        mpath = os.path.join(CMAP_DIR, f"{acc}.contig_map.tsv")
        rows = read_tsv(mpath)
        a2f = {}
        for r in rows:
            key = (acc, r["frozen_contig_id"])
            contigs[key] = {"accession": acc, "contig_id": r["frozen_contig_id"],
                            "annotation_contig_id": r["annotation_contig_id"],
                            "length_bp": int(r["length_bp"]),
                            "topology": "unknown", "topology_evidence": ""}
            a2f[r["annotation_contig_id"]] = r["frozen_contig_id"]
        if len(set(a2f.values())) != len(a2f):
            qc["warnings"].append(f"{acc}: contig map not bijective")
        ann2frozen[acc] = a2f
        # topology from Bakta fna headers
        adir = resolve(ann[acc]["annotation_dir"])
        fna = os.path.join(adir, f"{acc}.fna")
        topo = {}
        for line in open(fna):
            if line.startswith(">"):
                name = line[1:].split()[0]
                m = re.search(r"\[topology=(\w+)\]", line)
                topo[name] = m.group(1) if m else "unknown"
                if "topology" not in line:
                    qc["warnings"].append(f"{acc}: no topology tag for {name}")
        for key, c in contigs.items():
            if key[0] != acc:
                continue
            t = topo.get(c["annotation_contig_id"], "unknown")
            c["topology"] = "circular" if t == "circular" else ("linear" if t == "linear" else "unknown")
            if t == "circular":
                c["topology_evidence"] = f"bakta_fna_header [topology=circular] {os.path.basename(fna)}"
            elif t == "linear":
                c["topology_evidence"] = f"bakta_fna_header [topology=linear] {os.path.basename(fna)}"
    n_circ = sum(1 for c in contigs.values() if c["topology"] == "circular")
    qc["n_contigs"] = len(contigs)
    qc["n_circular_contigs"] = n_circ

    # ---- genes from occurrence_map ----
    genes = []
    ginv = Counter()
    for r in read_tsv(OCC):
        acc = r["accession"]
        # NOTE: occurrence_map coordinates are empirically 0-based half-open
        # (PGAP2 total.involved_annot [start:end) convention), NOT 1-based
        # closed: verified 2026-10-15 against the Bakta GFF for all
        # 328,376 rows (start+1 == GFF start, end == GFF end, 0 exceptions;
        # the 55 rows with start=0 all correspond to GFF start=1).
        rec = {"accession": acc, "gene_id": r["gene_id"], "cluster": r["cluster"],
               "partition": part.get(r["cluster"], "NA"),
               "orig_contig_id": r["contig"],
               "orig_coordinate_system": "0-based_half-open_pgap2_involved_annot(verified_vs_gff)",
               "orig_start": r["start"], "orig_end": r["end"], "strand": r["strand"],
               "contig_id": "NA", "start": "NA", "end": "NA", "valid": "1", "invalid_reason": ""}
        if r["cluster"] not in part:
            rec["valid"] = "0"; rec["invalid_reason"] = "cluster_not_in_partition_table"
            ginv[rec["invalid_reason"]] += 1
            genes.append(rec); continue
        frz = ann2frozen.get(acc, {}).get(r["contig"])
        if frz is None:
            rec["valid"] = "0"; rec["invalid_reason"] = "contig_not_in_contig_map"
            ginv[rec["invalid_reason"]] += 1
            genes.append(rec); continue
        try:
            s1, e1 = int(r["start"]), int(r["end"])
        except ValueError:
            rec["valid"] = "0"; rec["invalid_reason"] = "coordinate_not_integer"
            ginv[rec["invalid_reason"]] += 1
            genes.append(rec); continue
        L = contigs[(acc, frz)]["length_bp"]
        if not (0 <= s1 < e1 <= L):
            rec["valid"] = "0"
            rec["invalid_reason"] = f"coordinate_out_of_range(0..{L})"
            ginv["coordinate_out_of_range"] += 1
            genes.append(rec); continue
        rec["contig_id"] = frz
        rec["start"] = s1
        rec["end"] = e1
        genes.append(rec)
    qc["genes_total"] = len(genes)
    qc["genes_invalid"] = dict(ginv)

    # ---- local MGE union sources ----
    mge = []   # local only
    wc = []    # whole-contig predictions
    minv = Counter()
    ledger = []
    mid = 0

    def add_mge(acc, contig_frozen, klass, tool, s0, e0, os1, oe1, osys, extra=""):
        nonlocal mid
        L = contigs.get((acc, contig_frozen), {}).get("length_bp")
        if L is None:
            minv[f"{tool}:contig_not_in_frozen_map"] += 1
            return
        if not (0 <= s0 < e0 <= L):
            minv[f"{tool}:coordinate_out_of_range"] += 1
            return
        mid += 1
        mge.append({"accession": acc, "contig_id": contig_frozen,
                    "mge_id": f"M8MGE{mid:07d}", "mge_class": klass,
                    "source_tool": tool, "start": s0, "end": e0,
                    "orig_start": os1, "orig_end": oe1, "orig_coordinate_system": osys,
                    "extra": extra})

    for acc in accs:
        c = cov[acc]
        # --- ISEScan ---
        p = resolve(c["isescan_evidence"])
        n = 0
        tsv = None
        if p:
            cand = os.path.join(p, "fasta", f"{acc}.fna.tsv")
            tsv = cand if os.path.exists(cand) else None
        if tsv:
            for r in read_tsv(tsv):
                try:
                    s1, e1 = int(r["isBegin"]), int(r["isEnd"])
                except (KeyError, ValueError):
                    minv["isescan:unparseable_row"] += 1
                    continue
                if e1 <= s1:
                    minv["isescan:end_le_start"] += 1
                    continue
                before = mid
                add_mge(acc, r["seqID"], f"IS:{r.get('family','IS')}", "isescan",
                        s1 - 1, e1, s1, e1, "1-based_closed_frozen",
                        f"type={r.get('type','NA')};ncopy4is={r.get('ncopy4is','NA')}")
                n += int(mid > before)
        ledger.append({"accession": acc, "tool": "isescan", "resolved_path": tsv or "NA",
                       "run_status": c["isescan_final"], "n_records_used": n})
        # --- IntegronFinder ---
        p = resolve(c["integron_finder_evidence"])
        n = 0
        integ = None
        if p:
            hits = []
            for dirpath, _, files in os.walk(p):
                for f in files:
                    if f.endswith(".integrons"):
                        hits.append(os.path.join(dirpath, f))
            for integ in sorted(hits):
                with open(integ) as fh:
                    for line in fh:
                        if line.startswith("#") or not line.strip():
                            continue
                        f = line.rstrip("\n").split("\t")
                        # integron_finder 2.x .integrons: ID_replicon strand pos_beg pos_end ... type
                        if len(f) >= 8 and f[2].lstrip("-").isdigit() and f[3].lstrip("-").isdigit():
                            s1, e1 = int(f[2]), int(f[3])
                            before = mid
                            add_mge(acc, f[0], f"integron:{f[-1]}", "integron_finder",
                                    s1 - 1, e1, s1, e1, "1-based_closed_frozen")
                            n += int(mid > before)
                        else:
                            minv["integron_finder:unparseable_row"] += 1
        ledger.append({"accession": acc, "tool": "integron_finder",
                       "resolved_path": p or "NA", "run_status": c["integron_finder_final"],
                       "n_records_used": n})
        # --- geNomad provirus (local) ---
        p = resolve(c["genomad_evidence"])
        nprov = 0
        pv_tsv = os.path.join(p, f"{acc}_find_proviruses", f"{acc}_provirus.tsv") if p else None
        if pv_tsv and os.path.exists(pv_tsv):
            for r in read_tsv(pv_tsv):
                src = r.get("source_seq", "")
                try:
                    s1, e1 = int(r["start"]), int(r["end"])
                except (KeyError, ValueError):
                    minv["genomad_provirus:unparseable_row"] += 1
                    continue
                before = mid
                add_mge(acc, src, "provirus:genomad", "genomad_find_proviruses",
                        s1 - 1, e1, s1, e1, "1-based_closed_frozen_seqname_offset",
                        f"seq_name={r.get('seq_name','NA')};n_genes={r.get('n_genes','NA')}")
                nprov += int(mid > before)
        # whole-contig: plasmid_summary + virus_summary (coordinates==NA)
        nwp = nwv = 0
        psf = os.path.join(p, f"{acc}_summary", f"{acc}_plasmid_summary.tsv") if p else None
        if psf and os.path.exists(psf):
            for r in read_tsv(psf):
                wc.append({"accession": acc, "contig_id": r["seq_name"], "tool": "genomad",
                           "prediction": "plasmid_contig",
                           "detail": f"plasmid_score={r.get('plasmid_score','NA')};fdr={r.get('fdr','NA')};length={r.get('length','NA')}"})
                nwp += 1
        vsf = os.path.join(p, f"{acc}_summary", f"{acc}_virus_summary.tsv") if p else None
        nprov_summ = 0
        if vsf and os.path.exists(vsf):
            for r in read_tsv(vsf):
                coord = r.get("coordinates", "NA")
                if coord and coord != "NA" and "-" in coord:
                    nprov_summ += 1  # duplicates find_proviruses; cross-check only
                else:
                    wc.append({"accession": acc, "contig_id": r["seq_name"], "tool": "genomad",
                               "prediction": "virus_contig",
                               "detail": f"virus_score={r.get('virus_score','NA')};fdr={r.get('fdr','NA')};length={r.get('length','NA')}"})
                    nwv += 1
        if nprov_summ != nprov:
            qc["warnings"].append(f"{acc}: provirus.tsv n={nprov} vs virus_summary provirus n={nprov_summ}")
        ledger.append({"accession": acc, "tool": "genomad", "resolved_path": p or "NA",
                       "run_status": c["genomad_final"],
                       "n_records_used": nprov,
                       "note": f"whole_contig_plasmid={nwp};whole_contig_virus={nwv};provirus_xcheck={nprov_summ}"})
        # --- MOB-suite whole-contig plasmids ---
        p = resolve(c["mob_suite_evidence"])
        nmob = 0
        cr = os.path.join(p, "contig_report.txt") if p else None
        if cr and os.path.exists(cr):
            for r in read_tsv(cr):
                if r.get("molecule_type") == "plasmid":
                    wc.append({"accession": acc, "contig_id": r["contig_id"].split()[0],
                               "tool": "mob_suite", "prediction": "plasmid_contig",
                               "detail": f"primary_cluster_id={r.get('primary_cluster_id','NA')};predicted_mobility={r.get('predicted_mobility','NA')};size={r.get('size','NA')}"})
                    nmob += 1
        ledger.append({"accession": acc, "tool": "mob_suite", "resolved_path": cr or "NA",
                       "run_status": c["mob_suite_final"], "n_records_used": nmob,
                       "note": "whole_contig only"})
    qc["mge_local_rows"] = len(mge)
    qc["mge_local_invalid"] = dict(minv)
    qc["whole_contig_rows"] = len(wc)

    # ---- IslandPath GIs (annotation contig names -> frozen) ----
    gi = []
    ginv2 = Counter()
    gid = 0
    for acc in accs:
        c = cov[acc]
        if c["islandpath_final"] == "rerun_success" or not resolve(c["islandpath_evidence"]):
            idir = RERUN_IP
            src = "phaseB_v7_rerun"
        else:
            idir = os.path.join(resolve(c["islandpath_evidence"]), "islands")
            src = "revision_20261001"
        files = sorted(f for f in os.listdir(idir) if f.startswith(f"{acc}__"))
        if not files:
            qc["warnings"].append(f"{acc}: no islandpath files in {idir}")
        n = 0
        for f in files:
            ann_ctg = f[len(acc) + 2:-4]
            frz = ann2frozen[acc].get(ann_ctg)
            for r in csv.reader(open(os.path.join(idir, f)), delimiter="\t"):
                if not r or r[0].startswith("#"):
                    continue
                if len(r) < 9 or r[2] != "genomic_island":
                    ginv2["unexpected_row"] += 1
                    continue
                if r[0] != ann_ctg:
                    ginv2["seqid_filename_mismatch"] += 1
                    continue
                if frz is None:
                    ginv2["contig_not_in_map"] += 1
                    continue
                s1, e1 = int(r[3]), int(r[4])
                L = contigs[(acc, frz)]["length_bp"]
                if not (1 <= s1 < e1 <= L):
                    ginv2["coordinate_out_of_range"] += 1
                    continue
                gid += 1
                n += 1
                gi.append({"accession": acc, "contig_id": frz, "gi_id": f"M8GI{gid:06d}",
                           "start": s1 - 1, "end": e1, "strand": r[6],
                           "orig_feature_id": r[8].replace("ID=", "").split(";")[0],
                           "orig_start": s1, "orig_end": e1,
                           "orig_coordinate_system": "1-based_closed_annotation_contig",
                           "annotation_contig_id": ann_ctg, "source": src})
        ledger.append({"accession": acc, "tool": "islandpath", "resolved_path": idir,
                       "run_status": c["islandpath_final"], "n_records_used": n,
                       "note": f"source={src};n_contig_files={len(files)}"})
    qc["gi_rows"] = len(gi)
    qc["gi_invalid"] = dict(ginv2)

    # ---- full verification: occurrence_map coords vs Bakta GFF CDS ----
    gff_check = Counter()
    occ_by_acc = defaultdict(dict)
    for r in read_tsv(OCC):
        occ_by_acc[r["accession"]][r["gene_id"]] = (r["contig"], int(r["start"]), int(r["end"]), r["strand"])
    for acc, genes_acc in occ_by_acc.items():
        gff = resolve(ann[acc]["gff_path_resolved"])
        g = {}
        for line in open(gff):
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "CDS":
                continue
            m = re.search(r"ID=([^;]+)", f[8])
            if m:
                g[m.group(1)] = (f[0], int(f[3]), int(f[4]), f[6])
        for gid, (ctg, s, e, strand) in genes_acc.items():
            gc, gs, ge, gstrand = g.get(gid, (None, None, None, None))
            if ctg == gc and s + 1 == gs and e == ge and strand == gstrand:
                gff_check["match_0based_halfopen"] += 1
            elif ctg == gc and s == gs and e == ge and strand == gstrand:
                gff_check["match_1based_closed"] += 1
            else:
                gff_check["mismatch"] += 1
    qc["occurrence_vs_gff_coordinate_check"] = dict(gff_check)

    # ---- write ----
    def dump(name, rows, fields):
        with open(os.path.join(OUT, name), "w", newline="") as fh:
            w = csv.DictWriter(fh, delimiter="\t", fieldnames=fields, extrasaction="ignore")
            w.writeheader()
            for r in rows:
                w.writerow(r)

    dump("contigs.tsv", sorted(contigs.values(), key=lambda c: (c["accession"], c["annotation_contig_id"])),
         ["accession", "contig_id", "annotation_contig_id", "length_bp", "topology", "topology_evidence"])
    dump("genes.tsv", genes,
         ["accession", "gene_id", "cluster", "partition", "contig_id", "start", "end", "strand",
          "orig_contig_id", "orig_start", "orig_end", "orig_coordinate_system", "valid", "invalid_reason"])
    dump("mge_local.tsv", mge,
         ["accession", "contig_id", "mge_id", "mge_class", "source_tool", "start", "end",
          "orig_start", "orig_end", "orig_coordinate_system", "extra"])
    dump("whole_contig.tsv", wc,
         ["accession", "contig_id", "tool", "prediction", "detail"])
    dump("gi.tsv", gi,
         ["accession", "contig_id", "gi_id", "start", "end", "strand", "orig_feature_id",
          "orig_start", "orig_end", "orig_coordinate_system", "annotation_contig_id", "source"])
    dump("tool_coverage_ledger.tsv", ledger,
         ["accession", "tool", "resolved_path", "run_status", "n_records_used", "note"])

    # ---- cross-check vs historical 169-cohort tables (same raw files) ----
    old_mge = defaultdict(int)
    oldm = os.path.join(SOURCE_ROOT, "results_revision_20261001/08_mge_gi/tables/mge.tsv")
    if os.path.exists(oldm):
        for r in read_tsv(oldm):
            if r["whole_contig"] != "1":
                old_mge[r["genome"]] += 1
    new_mge = Counter(r["accession"] for r in mge)
    diff = {a: (old_mge.get(a, "NA"), new_mge.get(a, 0)) for a in accs
            if old_mge.get(a) is not None and old_mge.get(a, 0) != new_mge.get(a, 0)}
    qc["xcheck_old169_mge_local_count_diffs"] = diff
    old_gi = defaultdict(int)
    oldg = os.path.join(SOURCE_ROOT, "results_revision_20261001/08_mge_gi/tables/gi.tsv")
    if os.path.exists(oldg):
        for r in read_tsv(oldg):
            old_gi[r["genome"]] += 1
    new_gi = Counter(r["accession"] for r in gi)
    gdiff = {a: (old_gi.get(a, "NA"), new_gi.get(a, 0)) for a in accs
             if old_gi.get(a) is not None and old_gi.get(a, 0) != new_gi.get(a, 0)}
    qc["xcheck_old169_gi_count_diffs"] = gdiff

    with open(os.path.join(OUT, "build_qc.json"), "w") as fh:
        json.dump(qc, fh, indent=2, ensure_ascii=False)
    print(json.dumps({k: v for k, v in qc.items() if not isinstance(v, dict) or len(str(v)) < 400},
                     indent=2, ensure_ascii=False, default=str))
    print("MGE local rows:", len(mge), "| whole-contig rows:", len(wc), "| GI rows:", len(gi),
          "| genes:", len(genes), "| contigs:", len(contigs), "| circular:", n_circ)
    print("xcheck old169 MGE diffs:", len(diff), "| GI diffs:", len(gdiff))


if __name__ == "__main__":
    main()
