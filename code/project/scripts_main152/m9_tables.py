#!/usr/bin/env python3
"""M9 (tables part): cohort flow, reference ANI interface, TableS1 consistency,
TableS2 species x source composition.
"""
import csv
import os
import shutil
import sys
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import RESULTS, PACKET, read_roster, read_tsv, write_tsv

OUT_T = os.path.join(RESULTS, "11_tables")
PRE = os.path.join(RESULTS, "00_preflight")


def main():
    roster = read_roster()
    meta = {r["strain"]: r for r in read_tsv(os.path.join(PACKET, "00_input/cohort/main152_metadata.tsv"))}
    tetm = {r["accession"]: r["tetM_status"]
            for r in read_tsv(os.path.join(RESULTS, "07_amr_tetm/amr_tetm_per_genome.tsv"))}
    am = {r["accession"]: r for r in read_tsv(os.path.join(PRE, "assembly_metrics.tsv"))}

    # F1_FLOW
    flow = read_tsv(os.path.join(PACKET, "00_input/cohort/screening_flow_counts.tsv"))
    rules = {
        "candidate_assemblies": "genus-wide freeze snapshot 2026-10-01, 1303 records",
        "type_reference_membership_pass": "ANI>=95% and min bidirectional AF>=0.70 to a type-anchored reference channel with margin>=0.10 pp",
        "after_provenance_and_QC": "provenance corrections + CheckM2 completeness>=90%, contamination<=5%, N50>=10 kb",
        "dereplicated_representatives_before_final_reference_assignment_filter": "independent-unit dereplication representatives",
        "excluded_non_unique_reference_assignment": "17 records without a unique reference assignment (ledger retained for traceability only)",
        "final_main_cohort": "152 genomes: subflava 31, perflava 114, flavescens 7",
    }
    write_tsv(os.path.join(PRE, "cohort_flow.tsv"), ["stage", "n", "rule"],
              [dict(stage=r["stage"], n=r["n"], rule=rules.get(r["stage"], "")) for r in flow])

    # F1_REFERENCE: freeze copy
    shutil.copy(os.path.join(PACKET, "00_input/cohort/main152_reference_ANI_evidence.tsv"),
                os.path.join(RESULTS, "00_freeze/main152_reference_ANI_evidence.tsv"))

    # TableS1 consistency: compare packet S1 against frozen metadata per key fields
    s1 = read_tsv(os.path.join(PACKET, "00_input/TableS1/Supplementary_Table_S1_main152_metadata.tsv"))
    diffs = []
    s1_by_acc = {}
    for r in s1:
        acc = r.get("assembly_accession") or r.get("accession") or r.get("strain") or ""
        s1_by_acc[acc] = r
    if len(s1_by_acc) != 152:
        diffs.append(f"TableS1 accession count {len(s1_by_acc)} != 152 (key column detection)")
    for acc in roster:
        if acc not in s1_by_acc:
            diffs.append(f"{acc} missing from TableS1")
    write_tsv(os.path.join(OUT_T, "TableS1_consistency_check.tsv"),
              ["item", "status"],
              [dict(item="row count == 152", status=len(s1_by_acc) == 152),
               dict(item="all roster accessions present", status=all(a in s1_by_acc for a in roster)),
               *[dict(item=d, status="DIFF") for d in diffs]])

    # TableS2: species x source
    groups = Counter(meta[a]["group_for_testing"] for a in roster)
    rows = []
    by_pair = {}
    for a in roster:
        m = meta[a]
        key = (m["group_for_testing"], m["study_namespace"], m["true_source_bioproject"])
        by_pair.setdefault(key, []).append(a)
    for (grp, ns, bp), accs in sorted(by_pair.items()):
        tet_c = Counter(tetm[a] for a in accs)
        part_ev = Counter(meta[a]["participant_status"] for a in accs)
        asm = Counter(meta[a]["assembly_level"] for a in accs)
        rows.append(dict(
            species_group=grp, study_namespace=ns, source_BioProject=bp,
            n_genomes=len(accs), species_total=groups[grp],
            within_species_proportion=f"{len(accs) / groups[grp]:.6f}",
            participant_status_counts="; ".join(f"{k}={v}" for k, v in sorted(part_ev.items())),
            assembly_level_counts="; ".join(f"{k}={v}" for k, v in sorted(asm.items())),
            tetM_positive=tet_c.get("positive_complete", 0),
            tetM_negative=tet_c.get("negative", 0),
            tetM_indeterminate=tet_c.get("indeterminate", 0),
            tetM_missing=tet_c.get("missing", 0)))
    write_tsv(os.path.join(OUT_T, "TableS2_species_source.tsv"),
              ["species_group", "study_namespace", "source_BioProject", "n_genomes",
               "species_total", "within_species_proportion", "participant_status_counts",
               "assembly_level_counts", "tetM_positive", "tetM_negative",
               "tetM_indeterminate", "tetM_missing"], rows)
    print(f"TableS2 rows: {len(rows)}; flow stages: {len(flow)}; S1 diffs: {len(diffs)}")


if __name__ == "__main__":
    main()
