#!/usr/bin/env python3
"""M8 (main152): fixed-rule GI representative case selection + coordinates.

Selection rule (fixed before looking at any overlap statistic; adapted from
the historical d3_select_gi_examples.py rule to the new 152 dictionary):
  1. GI length 10-60 kb (visible detail at figure scale);
  2. carrier contig length >= 200 kb OR contig has explicit circular
     topology evidence (Bakta [topology=circular]);
  3. flanked on BOTH sides, within the same contig, by >=2 Strict_core
     anchor genes (global 152 partition) within 20 kb;
  4. at most one case per genome; genomes from different species groups are
     preferred (one per group before a second from any group);
  5. rank by (circular contig first, then longer contig, then accession,
     then GI start); take top 3.
No selection by MGE content, overlap significance, or visual salience.

For each selected GI: one feature_type=genomic_island row plus one row per
CDS overlapping the GI (coordinates 0-based half-open, unified frozen
contig_id), with product from the Bakta GFF and evidence strings.

Usage: m8_gi_cases.py <tables_dir> <outdir>
"""
from pathlib import Path
import csv
import os
import re
import sys
from collections import defaultdict

ROOT = str(Path(os.environ["NEISSERIA_PROJECT_ROOT"]))
RES = os.path.join(ROOT, "results_revision_20261015_main152_species")
SOURCE_ROOT = os.environ.get("NEISSERIA_SOURCE_ROOT", ROOT)
META = os.path.join(RES, "00_freeze/input_packet_snapshot/00_input/cohort/main152_metadata.tsv")
ANN_MAN = os.path.join(RES, "01_annotation/annotation_manifest_main152.tsv")

MIN_GI, MAX_GI = 10_000, 60_000
MIN_CONTIG = 200_000
FLANK = 20_000
MIN_ANCHORS = 2
N_CASES = 3


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def resolve(p):
    if os.path.exists(p):
        return p
    return os.path.join(SOURCE_ROOT, p.split("Neisseria_subflava/")[-1])


def main():
    tables, outdir = sys.argv[1], sys.argv[2]
    meta = {r["strain"]: r for r in read_tsv(META)}
    contigs = {(r["accession"], r["contig_id"]): r
               for r in read_tsv(os.path.join(tables, "contigs.tsv"))}
    genes = [r for r in read_tsv(os.path.join(tables, "genes.tsv")) if r["valid"] == "1"]
    gis = read_tsv(os.path.join(tables, "gi.tsv"))
    gi_mge = {r["gi_id"]: r for r in read_tsv(os.path.join(os.path.dirname(tables), "gi_mge_overlap.tsv"))}

    core_by_ctg = defaultdict(list)
    for r in genes:
        if r["partition"] == "Strict_core":
            core_by_ctg[(r["accession"], r["contig_id"])].append((int(r["start"]), int(r["end"])))
    for k in core_by_ctg:
        core_by_ctg[k].sort()

    cands = []
    for g in gis:
        key = (g["accession"], g["contig_id"])
        c = contigs.get(key)
        if c is None:
            continue
        L = int(c["length_bp"])
        circ = c["topology"] == "circular"
        gs, ge = int(g["start"]), int(g["end"])
        glen = ge - gs
        if not (MIN_GI <= glen <= MAX_GI):
            continue
        if not circ and L < MIN_CONTIG:
            continue
        anchors = core_by_ctg.get(key, [])
        left = [a for a in anchors if a[1] <= gs and gs - a[1] <= FLANK]
        right = [a for a in anchors if a[0] >= ge and a[0] - ge <= FLANK]
        if len(left) < MIN_ANCHORS or len(right) < MIN_ANCHORS:
            continue
        cands.append({"gi": g, "contig_len": L, "circular": circ,
                      "n_left": len(left), "n_right": len(right),
                      "group": meta[g["accession"]]["group_for_testing"]})
    cands.sort(key=lambda x: (not x["circular"], -x["contig_len"],
                              x["gi"]["accession"], int(x["gi"]["start"])))
    picked, used_acc, used_grp = [], set(), set()
    for cand in cands:  # group-diversity pass first
        a, gr = cand["gi"]["accession"], cand["group"]
        if a in used_acc or gr in used_grp:
            continue
        picked.append(cand); used_acc.add(a); used_grp.add(gr)
        if len(picked) >= N_CASES:
            break
    if len(picked) < N_CASES:
        for cand in cands:
            a = cand["gi"]["accession"]
            if a in used_acc:
                continue
            picked.append(cand); used_acc.add(a)
            if len(picked) >= N_CASES:
                break

    genes_by_acc_ctg = defaultdict(list)
    for r in genes:
        genes_by_acc_ctg[(r["accession"], r["contig_id"])].append(r)
    ann = {r["accession"]: r for r in read_tsv(ANN_MAN)}

    def products(acc, wanted):
        gff = resolve(ann[acc]["gff_path_resolved"])
        out = {}
        for line in open(gff):
            line = line.rstrip("\n")
            if line.startswith("#") or "\tCDS\t" not in line:
                continue
            m = re.search(r"ID=([^;]+)", line)
            if m and m.group(1) in wanted:
                pm = re.search(r"product=([^;]+)", line)
                gm = re.search(r"(?:^|;)gene=([^;]+)", line)
                out[m.group(1)] = (pm.group(1) if pm else "NA",
                                   gm.group(1) if gm else "NA")
        return out

    rows = []
    for cand in picked:
        g = cand["gi"]
        acc, ctg = g["accession"], g["contig_id"]
        gs, ge = int(g["start"]), int(g["end"])
        gm = gi_mge.get(g["gi_id"], {})
        rows.append({
            "accession": acc, "feature_id": g["gi_id"],
            "feature_type": "genomic_island", "contig_id": ctg,
            "coordinate_system": "0-based_half-open",
            "start": gs, "end": ge, "strand": g["strand"],
            "product": "genomic_island",
            "evidence": (f"IslandPath-DIMOB 1.0.6 ({g['source']});orig_id={g['orig_feature_id']};"
                         f"contig_topology={contigs[(acc, ctg)]['topology']};"
                         f"contig_len={cand['contig_len']};gi_len={ge - gs};"
                         f"core_anchors_left={cand['n_left']};right={cand['n_right']};"
                         f"mge_overlap_union_bp={gm.get('mge_overlap_union_bp', 'NA')};"
                         f"selection_rule=fixed(len10-60kb,contig>=200kb_or_circular,"
                         f">=2 Strict_core anchors each flank<=20kb)")})
        inside = [r for r in genes_by_acc_ctg[(acc, ctg)]
                  if int(r["start"]) < ge and int(r["end"]) > gs]
        prods = products(acc, {r["gene_id"] for r in inside})
        for r in sorted(inside, key=lambda x: int(x["start"])):
            p, gn = prods.get(r["gene_id"], ("NA", "NA"))
            rows.append({
                "accession": acc, "feature_id": r["gene_id"],
                "feature_type": "CDS", "contig_id": ctg,
                "coordinate_system": "0-based_half-open",
                "start": r["start"], "end": r["end"], "strand": r["strand"],
                "product": p,
                "evidence": (f"Bakta 1.12.0 annotation;gene={gn};cluster={r['cluster']};"
                             f"partition={r['partition']};within_gi={g['gi_id']}")})

    out = os.path.join(outdir, "gi_case_coordinates.tsv")
    fields = ["accession", "feature_id", "feature_type", "contig_id",
              "coordinate_system", "start", "end", "strand", "product", "evidence"]
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    print(f"candidates={len(cands)} picked={len(picked)} rows={len(rows)} -> {out}")
    for cand in picked:
        g = cand["gi"]
        print(" ", g["accession"], cand["group"], g["contig_id"], g["start"], g["end"],
              "circ" if cand["circular"] else f"len={cand['contig_len']}")


if __name__ == "__main__":
    main()
