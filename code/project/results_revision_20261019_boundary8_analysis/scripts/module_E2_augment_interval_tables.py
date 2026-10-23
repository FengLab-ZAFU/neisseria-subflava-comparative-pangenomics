#!/usr/bin/env python3
"""Module E2: augment interval instance tables with boundary8 scope columns.

Regenerates the delivered layouts of:
  S/09_interval_overlap/occurrence_overlap.tsv  (from module E output + metadata)
  S/09_interval_overlap/gi_mge_overlap.tsv      (from module E output + metadata)

NOTE: module E writes plain versions of both tables and runs after module D2,
so this augmentation must run last (pipeline order: D2 -> E -> E2).

Column semantics (verified byte-identical against the accepted 2026-10-19 delivery):
  display_group   metadata display_group (species group or unresolved_species)
  group           display_group for comparison144 genomes, empty for boundary8
  scope_pooled    TRUE for all 152 (pooled152 scope)
  scope_<group>   primary144 group membership flags
  included_primary144 / included_sensitivity152 / included_descriptive8 /
  included_participant_subset_from144   0/1 membership flags

Run after module_E_intervals.py and module_D2_defense.py.
"""
from pathlib import Path
import csv, os, sys

ROOT = str(Path(os.environ["NEISSERIA_PROJECT_ROOT"]))
P = os.path.join(ROOT, "boundary8_split_handoff_20261019")
S = os.path.join(ROOT, "results_revision_20261019_boundary8_analysis")
OUT09 = os.path.join(S, "09_interval_overlap")

GROUPS = ["subflava", "perflava", "flavescens"]


def read_tsv(path):
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(path, fieldnames, rows):
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow(r)


def main():
    meta_overall = {r["strain"]: r for r in read_tsv(os.path.join(P, "00_cohort/metadata_overall152.tsv"))}
    meta_comp = {r["strain"] for r in read_tsv(os.path.join(P, "00_cohort/metadata_comparison144.tsv"))}
    psm = {r["accession"]: r for r in read_tsv(os.path.join(S, "08_defense_mge_vfdb/participant_subset_manifest.tsv"))}

    occ = read_tsv(os.path.join(OUT09, "occurrence_overlap.tsv"))
    occ_fields = ["accession", "display_group", "group",
                  "scope_pooled", "scope_subflava", "scope_perflava", "scope_flavescens",
                  "included_primary144", "included_sensitivity152", "included_descriptive8",
                  "included_participant_subset_from144",
                  "gene_id", "cluster", "partition", "contig_id", "coordinate_system",
                  "start", "end", "strand", "gene_length",
                  "mge_overlap_union_bp", "mge_any_bp", "mge_half_gene",
                  "gi_overlap_union_bp", "gi_half_gene",
                  "eligible_for_null", "contig_layout_status", "unmapped_invalid_reason",
                  "orig_contig_id", "orig_start", "orig_end", "orig_coordinate_system"]
    out = []
    for r in occ:
        acc = r["accession"]
        dg = meta_overall.get(acc, {}).get("display_group", "")
        in_comp = acc in meta_comp
        row = {k: r.get(k, "") for k in occ_fields}
        row["display_group"] = dg
        row["group"] = dg if in_comp else ""
        row["scope_pooled"] = "TRUE"
        for g in GROUPS:
            row[f"scope_{g}"] = "TRUE" if (in_comp and dg == g) else "FALSE"
        row["included_primary144"] = "1" if in_comp else "0"
        row["included_sensitivity152"] = "1" if acc in meta_overall else "0"
        row["included_descriptive8"] = "1" if (acc in meta_overall and not in_comp) else "0"
        row["included_participant_subset_from144"] = "1" if int(psm.get(acc, {}).get("selected", 0)) == 1 else "0"
        out.append(row)
    write_tsv(os.path.join(OUT09, "occurrence_overlap.tsv"), occ_fields, out)
    print(f"[E2] occurrence_overlap: {len(out)} rows")

    gi = read_tsv(os.path.join(OUT09, "gi_mge_overlap.tsv"))
    gi_fields = ["accession", "display_group", "group", "contig_id", "gi_id", "start", "end",
                 "gi_length", "mge_overlap_union_bp", "overlap_over_gi_length",
                 "gi_half_covered_by_mge", "included_primary144", "included_sensitivity152",
                 "included_descriptive8"]
    gi_out = []
    for r in gi:
        acc = r["accession"]
        dg = meta_overall.get(acc, {}).get("display_group", "")
        in_comp = acc in meta_comp
        row = {k: r.get(k, "") for k in gi_fields}
        row["display_group"] = dg
        row["group"] = dg if in_comp else ""
        row["included_primary144"] = "1" if in_comp else "0"
        row["included_sensitivity152"] = "1" if acc in meta_overall else "0"
        row["included_descriptive8"] = "1" if (acc in meta_overall and not in_comp) else "0"
        gi_out.append(row)
    write_tsv(os.path.join(OUT09, "gi_mge_overlap.tsv"), gi_fields, gi_out)
    print(f"[E2] gi_mge_overlap: {len(gi_out)} rows")


if __name__ == "__main__":
    main()
