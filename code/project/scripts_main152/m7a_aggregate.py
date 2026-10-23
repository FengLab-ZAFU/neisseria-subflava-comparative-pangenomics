#!/usr/bin/env python3
"""M7a (main152, 2026-10-15): per-genome aggregation of DefenseFinder, geNomad,
IslandPath-DIMOB, MOB-suite, IntegronFinder and ISEScan evidence for the frozen
152-genome cohort.

Inputs (all read-only; historical paths resolved via common.resolve_path,
2026-10-15 cleanup moved results_revision_2026100{1,6,8}* into
NEISSERIA_SOURCE_ROOT):
  - roster/groups: 00_input/cohort/main152_accessions.txt + main152_metadata.tsv
  - per-tool evidence index: 00_input/reference_sources/per_tool_coverage_historical_main152.tsv
  - assembly metrics: results_revision_20261015_main152_species/00_preflight/assembly_metrics.tsv

Counting rules (identical to scripts/revision_20261001/d2_defense_summary.py +
scripts/revision_20261010_phaseB_v7/p5c_01_build_feature_tables.py, minus the
contig-namespace remapping which is M8's job):
  DefenseFinder: systems rows in <acc>_defense_finder_systems.tsv; family =
    column "type"; partial = per-system min(sys_wholeness) < 1.0 from the genes
    file; complete = total - partial; complete families = distinct type among
    complete systems.
  ISEScan: rows of fasta/<acc>.fna.tsv with isEnd > isBegin (one row = one IS).
  IntegronFinder: non-comment data lines of Results_Integron_Finder_<acc>/<acc>.integrons.
  geNomad: <acc>_summary/<acc>_virus_summary.tsv rows with coordinates != NA =
    provirus (local); other virus rows = whole-contig virus;
    <acc>_plasmid_summary.tsv rows = whole-contig plasmid.
  MOB-suite: contig_report.txt rows with molecule_type == "plasmid" (whole-contig).
  IslandPath: islands/<acc>__*.gff records (1-based inclusive -> 0-based
    half-open); GI total length = union of intervals per genome; fraction =
    union_bp / total_length_bp.

Run status vocabulary (never infer): success_hits / success_empty when the key
evidence file(s) parse cleanly; failed when files exist but parsing errors;
missing when no evidence file is found. Failures are never written as 0.

Outputs (results_revision_20261015_main152_species/08_defense_mge_vfdb/):
  tool_coverage.tsv                 152 x 6 run ledger with resolved path + hash
  m7a_evidence_files.tsv            every parsed evidence file + sha256
  defense_per_genome.tsv            152 rows (F6_DEFENSE contract + raw label)
  intermediate/defense_type_counts.tsv  accession x system_type complete/partial
  per_genome_features.tsv           152 rows of MGE/GI counts + prevalence flags
  m7a_coordinate_sources.tsv        coordinate field provenance ledger for M8
Usage: python3 scripts_main152/m7a_aggregate.py
"""
import csv
import glob
import hashlib
import json
import os
import sys
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import (ROOT, RESULTS, PACKET, read_roster, read_tsv, write_tsv,
                    sha256_file, resolve_path, RUN_ID)

OUT = os.path.join(RESULTS, "08_defense_mge_vfdb")
INTER = os.path.join(OUT, "intermediate")
os.makedirs(INTER, exist_ok=True)

TOOLS = ["defense_finder", "genomad", "islandpath", "mob_suite",
         "integron_finder", "isescan"]

roster = read_roster()
meta = {r["strain"]: r for r in read_tsv(os.path.join(PACKET, "00_input/cohort/main152_metadata.tsv"))}
amet = {r["accession"]: r for r in read_tsv(os.path.join(RESULTS, "00_preflight/assembly_metrics.tsv"))}
cov = {r["accession"]: r for r in read_tsv(os.path.join(PACKET, "00_input/reference_sources/per_tool_coverage_historical_main152.tsv"))}
assert set(roster) == set(meta) == set(cov) == set(amet)

GROUP_ORDER = ["subflava", "perflava", "flavescens"]


def assembly_class(level):
    return "Complete Genome" if level == "Complete Genome" else "Draft"


def combined_digest(file_hashes):
    """file_hashes: list of (relpath, sha256) -> single digest string."""
    blob = "".join(f"{rp}:{h}\n" for rp, h in sorted(file_hashes))
    return hashlib.sha256(blob.encode("utf-8")).hexdigest()


evidence_files = []   # accession, tool, path(relative to ROOT), sha256, n_bytes
def register_files(acc, tool, paths):
    fh = []
    for p in paths:
        h = sha256_file(p)
        rel = os.path.relpath(p, ROOT)
        evidence_files.append({"accession": acc, "tool": tool,
                               "evidence_file": rel, "sha256": h,
                               "n_bytes": os.path.getsize(p)})
        fh.append((rel, h))
    return fh


coverage = {}       # (acc, tool) -> dict
defense = {}        # acc -> dict
features = defaultdict(dict)  # acc -> counts
type_counts = defaultdict(Counter)        # acc -> Counter of COMPLETE types
type_counts_partial = defaultdict(Counter)

qc_notes = []

for acc in roster:
    crow = cov[acc]
    # ---------------- DefenseFinder ----------------
    ev, note = resolve_path(crow["defense_finder_evidence"])
    sf = os.path.join(ev, f"{acc}_defense_finder_systems.tsv") if ev else None
    gf = os.path.join(ev, f"{acc}_defense_finder_genes.tsv") if ev else None
    run_status, n_sys, n_partial, n_complete = "missing", None, None, None
    n_complete_fam = None
    if sf and os.path.exists(sf):
        try:
            fh = register_files(acc, "defense_finder",
                                [p for p in (sf, gf) if p and os.path.exists(p)])
            whole = {}
            if os.path.exists(gf):
                for r in read_tsv(gf):
                    sid = r.get("sys_id", "")
                    if not sid:
                        continue
                    try:
                        w = float(r.get("sys_wholeness", "") or 1.0)
                    except ValueError:
                        w = 1.0
                    whole[sid] = min(w, whole.get(sid, w))
            srows = read_tsv(sf)
            n_sys = len(srows)
            n_partial = 0
            complete_cnt = Counter()
            missing_wholeness = 0
            for r in srows:
                sid = r.get("sys_id", "")
                w = whole.get(sid)
                if w is None:
                    missing_wholeness += 1
                    w = 1.0
                if w < 1.0:
                    n_partial += 1
                    type_counts_partial[acc][r.get("type", "")] += 1
                else:
                    complete_cnt[r.get("type", "")] += 1
            if missing_wholeness:
                qc_notes.append({"accession": acc, "tool": "defense_finder",
                                 "note": f"{missing_wholeness} sys_id in systems file absent from genes file; treated as wholeness=1.0"})
            n_complete = n_sys - n_partial
            n_complete_fam = len(complete_cnt)
            type_counts[acc] = complete_cnt
            run_status = "success_hits" if n_sys > 0 else "success_empty"
        except Exception as e:
            run_status = "failed"
            qc_notes.append({"accession": acc, "tool": "defense_finder",
                             "note": f"parse exception: {e!r}"})
            fh = []
    else:
        fh = []
    coverage[(acc, "defense_finder")] = {
        "run_status": run_status, "evidence_path": os.path.relpath(ev, ROOT) if ev else "NA",
        "path_resolution": note, "historical_status": crow["defense_finder_final"],
        "evidence_files_sha256": combined_digest(fh) if fh else "NA",
        "n_records": n_sys if n_sys is not None else "NA"}
    defense[acc] = {"run_status": run_status, "n_complete_systems": n_complete,
                    "n_complete_families": n_complete_fam, "n_partial_systems": n_partial,
                    "n_systems_total": n_sys}

    # ---------------- geNomad ----------------
    ev, note = resolve_path(crow["genomad_evidence"])
    pf = os.path.join(ev, f"{acc}_summary", f"{acc}_plasmid_summary.tsv") if ev else None
    vf = os.path.join(ev, f"{acc}_summary", f"{acc}_virus_summary.tsv") if ev else None
    n_plasmid = n_provirus = n_virus_contig = None
    provirus_suffix_fail = 0
    run_status = "missing"
    if pf and os.path.exists(pf) and os.path.exists(vf):
        try:
            fh = register_files(acc, "genomad", [pf, vf])
            n_plasmid = len(read_tsv(pf))
            n_provirus = 0
            n_virus_contig = 0
            for r in read_tsv(vf):
                coord = r.get("coordinates", "NA")
                if coord and coord != "NA" and "-" in coord:
                    n_provirus += 1
                    sn = r.get("seq_name", "")
                    try:
                        S, E = sn.split("|provirus_")[1].split("_")
                        a, b = coord.split("-")[:2]
                        if not (int(S) == int(a) and int(E) == int(b)):
                            provirus_suffix_fail += 1
                    except (IndexError, ValueError):
                        provirus_suffix_fail += 1
                else:
                    n_virus_contig += 1
            if provirus_suffix_fail:
                qc_notes.append({"accession": acc, "tool": "genomad",
                                 "note": f"{provirus_suffix_fail} provirus rows with seq_name suffix != coordinates (counted as provirus; coordinate QC left to M8)"})
            n = n_plasmid + n_provirus + n_virus_contig
            run_status = "success_hits" if n > 0 else "success_empty"
        except Exception as e:
            run_status = "failed"
            qc_notes.append({"accession": acc, "tool": "genomad", "note": f"parse exception: {e!r}"})
            fh = []
    else:
        fh = []
    coverage[(acc, "genomad")] = {
        "run_status": run_status, "evidence_path": os.path.relpath(ev, ROOT) if ev else "NA",
        "path_resolution": note, "historical_status": crow["genomad_final"],
        "evidence_files_sha256": combined_digest(fh) if fh else "NA",
        "n_records": (n_plasmid + n_provirus + n_virus_contig) if n_plasmid is not None else "NA"}
    features[acc].update(genomad_plasmid_contigs=n_plasmid,
                         genomad_provirus=n_provirus,
                         genomad_virus_contigs=n_virus_contig,
                         genomad_provirus_suffix_mismatch=provirus_suffix_fail if n_provirus is not None else "NA")

    # ---------------- IslandPath ----------------
    ev_raw = crow["islandpath_evidence"]
    ev, note = resolve_path(ev_raw)
    files = []
    resolved_dir = ev
    if ev:
        isl_sub = os.path.join(ev, "islands")
        files = sorted(glob.glob(os.path.join(isl_sub, f"{acc}__*.gff")))
        if files:
            resolved_dir = isl_sub
        else:
            files = sorted(glob.glob(os.path.join(ev, f"{acc}__*.gff")))
    if not files:
        # rerun index rows point to <rerun_root>/<acc> which does not exist;
        # real outputs live in <rerun_root>/islands/
        parent, note2 = resolve_path(os.path.dirname(ev_raw))
        if parent:
            files = sorted(glob.glob(os.path.join(parent, "islands", f"{acc}__*.gff")))
            if files:
                resolved_dir = os.path.join(parent, "islands")
                note = note + ";evidence_dir_corrected_to_islands"
    n_gi = None
    gi_union_bp = None
    run_status = "missing"
    if files:
        try:
            fh = register_files(acc, "islandpath", files)
            intervals = defaultdict(list)
            n_gi = 0
            for f in files:
                for line in open(f, encoding="utf-8", errors="replace"):
                    if line.startswith("#"):
                        continue
                    p = line.rstrip("\n").split("\t")
                    if len(p) < 5:
                        continue
                    s, e = int(p[3]) - 1, int(p[4])  # 1-based inclusive -> 0-based half-open
                    if e > s >= 0:
                        intervals[p[0]].append((s, e))
                        n_gi += 1
            u = 0
            for c, ivs in intervals.items():
                ivs.sort()
                cs, ce = ivs[0]
                for s, e in ivs[1:]:
                    if s <= ce:
                        ce = max(ce, e)
                    else:
                        u += ce - cs
                        cs, ce = s, e
                u += ce - cs
            gi_union_bp = u
            run_status = "success_hits" if n_gi > 0 else "success_empty"
        except Exception as e:
            run_status = "failed"
            qc_notes.append({"accession": acc, "tool": "islandpath", "note": f"parse exception: {e!r}"})
            fh = []
    else:
        fh = []
    coverage[(acc, "islandpath")] = {
        "run_status": run_status,
        "evidence_path": os.path.relpath(resolved_dir, ROOT) if resolved_dir and files else "NA",
        "path_resolution": note, "historical_status": crow["islandpath_final"],
        "evidence_files_sha256": combined_digest(fh) if fh else "NA",
        "n_records": n_gi if n_gi is not None else "NA"}
    features[acc].update(islandpath_gi=n_gi, islandpath_gi_union_bp=gi_union_bp)

    # ---------------- MOB-suite ----------------
    ev, note = resolve_path(crow["mob_suite_evidence"])
    f = os.path.join(ev, "contig_report.txt") if ev else None
    n_mob = None
    run_status = "missing"
    if f and os.path.exists(f):
        try:
            fh = register_files(acc, "mob_suite", [f])
            n_mob = sum(1 for r in read_tsv(f) if r.get("molecule_type") == "plasmid")
            run_status = "success_hits" if n_mob > 0 else "success_empty"
        except Exception as e:
            run_status = "failed"
            qc_notes.append({"accession": acc, "tool": "mob_suite", "note": f"parse exception: {e!r}"})
            fh = []
    else:
        fh = []
    coverage[(acc, "mob_suite")] = {
        "run_status": run_status, "evidence_path": os.path.relpath(ev, ROOT) if ev else "NA",
        "path_resolution": note, "historical_status": crow["mob_suite_final"],
        "evidence_files_sha256": combined_digest(fh) if fh else "NA",
        "n_records": n_mob if n_mob is not None else "NA"}
    features[acc].update(mob_suite_plasmid_contigs=n_mob)

    # ---------------- IntegronFinder ----------------
    ev, note = resolve_path(crow["integron_finder_evidence"])
    f = os.path.join(ev, f"Results_Integron_Finder_{acc}", f"{acc}.integrons") if ev else None
    n_int = None
    run_status = "missing"
    if f and os.path.exists(f):
        try:
            fh = register_files(acc, "integron_finder", [f])
            n_int = 0
            for line in open(f, encoding="utf-8", errors="replace"):
                if line.startswith("#"):
                    continue
                p = line.rstrip("\n").split("\t")
                if len(p) >= 8 and p[0] and p[2].isdigit():
                    n_int += 1
            run_status = "success_hits" if n_int > 0 else "success_empty"
        except Exception as e:
            run_status = "failed"
            qc_notes.append({"accession": acc, "tool": "integron_finder", "note": f"parse exception: {e!r}"})
            fh = []
    else:
        fh = []
    coverage[(acc, "integron_finder")] = {
        "run_status": run_status, "evidence_path": os.path.relpath(ev, ROOT) if ev else "NA",
        "path_resolution": note, "historical_status": crow["integron_finder_final"],
        "evidence_files_sha256": combined_digest(fh) if fh else "NA",
        "n_records": n_int if n_int is not None else "NA"}
    features[acc].update(integron_finder_integrons=n_int)

    # ---------------- ISEScan ----------------
    ev, note = resolve_path(crow["isescan_evidence"])
    tsves = sorted(glob.glob(os.path.join(ev, "fasta", "*.tsv"))) if ev else []
    n_is = None
    run_status = "missing"
    if tsves:
        try:
            fh = register_files(acc, "isescan", tsves)
            n_is = 0
            for tsv in tsves:
                for r in read_tsv(tsv):
                    try:
                        s, e = int(r["isBegin"]), int(r["isEnd"])
                    except (KeyError, ValueError, TypeError):
                        continue
                    if e > s:
                        n_is += 1
            run_status = "success_hits" if n_is > 0 else "success_empty"
        except Exception as e:
            run_status = "failed"
            qc_notes.append({"accession": acc, "tool": "isescan", "note": f"parse exception: {e!r}"})
            fh = []
    else:
        fh = []
    coverage[(acc, "isescan")] = {
        "run_status": run_status, "evidence_path": os.path.relpath(ev, ROOT) if ev else "NA",
        "path_resolution": note, "historical_status": crow["isescan_final"],
        "evidence_files_sha256": combined_digest(fh) if fh else "NA",
        "n_records": n_is if n_is is not None else "NA"}
    features[acc].update(isescan_is=n_is)

# ---------------- tool_coverage.tsv ----------------
rows = []
for acc in roster:
    for t in TOOLS:
        c = coverage[(acc, t)]
        rows.append({"accession": acc, "species_group": meta[acc]["group_for_testing"],
                     "tool": t, **c})
write_tsv(os.path.join(OUT, "tool_coverage.tsv"),
          ["accession", "species_group", "tool", "run_status", "historical_status",
           "evidence_path", "path_resolution", "evidence_files_sha256", "n_records"], rows)

write_tsv(os.path.join(OUT, "m7a_evidence_files.tsv"),
          ["accession", "tool", "evidence_file", "sha256", "n_bytes"], evidence_files)

# ---------------- defense_per_genome.tsv (F6_DEFENSE) ----------------
rows = []
for acc in roster:
    d = defense[acc]
    m = meta[acc]
    a = amet[acc]
    ok = d["run_status"] in ("success_hits", "success_empty")
    rows.append({"accession": acc, "species_group": m["group_for_testing"],
                 "run_status": d["run_status"],
                 "n_complete_systems": d["n_complete_systems"] if ok else "NA",
                 "n_complete_families": d["n_complete_families"] if ok else "NA",
                 "n_partial_systems": d["n_partial_systems"] if ok else "NA",
                 "assembly_class": assembly_class(m["assembly_level"]),
                 "assembly_level_raw": m["assembly_level"],
                 "n_contigs": a["n_contigs"], "genome_size_bp": a["total_length_bp"]})
write_tsv(os.path.join(OUT, "defense_per_genome.tsv"),
          ["accession", "species_group", "run_status", "n_complete_systems",
           "n_complete_families", "n_partial_systems", "assembly_class",
           "assembly_level_raw", "n_contigs", "genome_size_bp"], rows)

# ---------------- intermediate defense_type_counts.tsv ----------------
all_types = sorted({t for acc in roster for t in type_counts[acc]} |
                   {t for acc in roster for t in type_counts_partial[acc]})
rows = []
for acc in roster:
    d = defense[acc]
    ok = d["run_status"] in ("success_hits", "success_empty")
    for t in all_types:
        rows.append({"accession": acc, "system_type": t,
                     "run_status": d["run_status"],
                     "n_complete_systems": type_counts[acc].get(t, 0) if ok else "NA",
                     "n_partial_systems": type_counts_partial[acc].get(t, 0) if ok else "NA"})
write_tsv(os.path.join(INTER, "defense_type_counts.tsv"),
          ["accession", "system_type", "run_status", "n_complete_systems", "n_partial_systems"], rows)

# ---------------- per_genome_features.tsv ----------------
def flag(v):
    return "NA" if v in (None, "NA") else (1 if v >= 1 else 0)

rows = []
for acc in roster:
    f = features[acc]
    m = meta[acc]
    a = amet[acc]
    size = int(a["total_length_bp"])
    gi_bp = f["islandpath_gi_union_bp"]
    row = {"accession": acc, "species_group": m["group_for_testing"],
           "isescan_run_status": coverage[(acc, "isescan")]["run_status"],
           "isescan_is_count": f["isescan_is"] if f["isescan_is"] is not None else "NA",
           "isescan_is_present": flag(f["isescan_is"]),
           "integron_finder_run_status": coverage[(acc, "integron_finder")]["run_status"],
           "integron_count": f["integron_finder_integrons"] if f["integron_finder_integrons"] is not None else "NA",
           "integron_present": flag(f["integron_finder_integrons"]),
           "genomad_run_status": coverage[(acc, "genomad")]["run_status"],
           "genomad_provirus_count": f["genomad_provirus"] if f["genomad_provirus"] is not None else "NA",
           "genomad_provirus_present": flag(f["genomad_provirus"]),
           "genomad_provirus_suffix_mismatch": f["genomad_provirus_suffix_mismatch"],
           "genomad_plasmid_contig_count": f["genomad_plasmid_contigs"] if f["genomad_plasmid_contigs"] is not None else "NA",
           "genomad_plasmid_contig_present": flag(f["genomad_plasmid_contigs"]),
           "genomad_virus_contig_count": f["genomad_virus_contigs"] if f["genomad_virus_contigs"] is not None else "NA",
           "genomad_virus_contig_present": flag(f["genomad_virus_contigs"]),
           "mob_suite_run_status": coverage[(acc, "mob_suite")]["run_status"],
           "mob_suite_plasmid_contig_count": f["mob_suite_plasmid_contigs"] if f["mob_suite_plasmid_contigs"] is not None else "NA",
           "mob_suite_plasmid_contig_present": flag(f["mob_suite_plasmid_contigs"]),
           "islandpath_run_status": coverage[(acc, "islandpath")]["run_status"],
           "islandpath_gi_count": f["islandpath_gi"] if f["islandpath_gi"] is not None else "NA",
           "islandpath_gi_present": flag(f["islandpath_gi"]),
           "gi_total_length_bp": gi_bp if gi_bp is not None else "NA",
           "gi_genome_fraction": ("%.8f" % (gi_bp / size)) if gi_bp is not None else "NA",
           "genome_size_bp": size}
    rows.append(row)
write_tsv(os.path.join(OUT, "per_genome_features.tsv"), list(rows[0].keys()), rows)

# ---------------- coordinate source ledger (for M8) ----------------
coord_rows = [
    {"tool": "isescan", "evidence_file_pattern": "<evidence>/fasta/<acc>.fna.tsv",
     "coordinate_fields": "seqID,isBegin,isEnd", "coordinate_system": "1-based closed on frozen-FASTA contig ids",
     "note": "one row per IS copy; convert to 0-based half-open as (isBegin-1, isEnd); frozen contig namespace (e.g. ACEN...), remap via contig_id_maps for M8"},
    {"tool": "integron_finder", "evidence_file_pattern": "<evidence>/Results_Integron_Finder_<acc>/<acc>.integrons",
     "coordinate_fields": "column1=replicon, pos_beg(col3), pos_end(col4), last column=integron type",
     "coordinate_system": "1-based closed on frozen-FASTA contig ids",
     "note": "convert (pos_beg-1, pos_end); frozen contig namespace"},
    {"tool": "genomad_provirus", "evidence_file_pattern": "<evidence>/<acc>_summary/<acc>_virus_summary.tsv",
     "coordinate_fields": "seq_name (parent_contig|provirus_S_E), coordinates (S-E on parent contig)",
     "coordinate_system": "1-based closed interval on the PARENT contig (nested provirus sequence must be mapped back)",
     "note": "convert (S-1, E); suffix check S/E vs coordinates done at count stage (mismatch count in per_genome_features genomad_provirus_suffix_mismatch); parent contig is frozen namespace"},
    {"tool": "genomad_virus_contig", "evidence_file_pattern": "<evidence>/<acc>_summary/<acc>_virus_summary.tsv",
     "coordinate_fields": "seq_name, length (coordinates == NA)", "coordinate_system": "whole contig, frozen namespace",
     "note": "whole-contig prediction only; never mix into local interval sets"},
    {"tool": "genomad_plasmid_contig", "evidence_file_pattern": "<evidence>/<acc>_summary/<acc>_plasmid_summary.tsv",
     "coordinate_fields": "seq_name, length", "coordinate_system": "whole contig, frozen namespace",
     "note": "whole-contig prediction only"},
    {"tool": "mob_suite", "evidence_file_pattern": "<evidence>/contig_report.txt",
     "coordinate_fields": "contig_id (first whitespace-delimited token), size; molecule_type==plasmid",
     "coordinate_system": "whole contig, frozen namespace",
     "note": "whole-contig prediction only"},
    {"tool": "islandpath", "evidence_file_pattern": "<evidence>/islands/<acc>__<contig>.gff (rerun set: <rerun_root>/islands/)",
     "coordinate_fields": "GFF columns 1 (seqid), 4 (start), 5 (end)",
     "coordinate_system": "GFF 1-based inclusive; seqid uses annotation namespace (contig_N)",
     "note": "convert (start-1, end); union per contig before length summaries"},
]
write_tsv(os.path.join(OUT, "m7a_coordinate_sources.tsv"),
          ["tool", "evidence_file_pattern", "coordinate_fields", "coordinate_system", "note"], coord_rows)

if qc_notes:
    write_tsv(os.path.join(OUT, "m7a_qc_notes.tsv"), ["accession", "tool", "note"], qc_notes)
else:
    write_tsv(os.path.join(OUT, "m7a_qc_notes.tsv"), ["accession", "tool", "note"], [])

# ---------------- console summary ----------------
print("run_id:", RUN_ID)
for t in TOOLS:
    print(t, dict(Counter(coverage[(a, t)]["run_status"] for a in roster)))
print("defense complete systems: min/med/max",
      min(d["n_complete_systems"] for d in defense.values()),
      sorted(d["n_complete_systems"] for d in defense.values())[len(roster)//2],
      max(d["n_complete_systems"] for d in defense.values()))
print("distinct complete types:", len(all_types))
print("qc_notes:", len(qc_notes))
print("evidence files registered:", len(evidence_files))
