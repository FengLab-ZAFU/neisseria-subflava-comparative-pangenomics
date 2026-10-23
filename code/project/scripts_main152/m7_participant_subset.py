#!/usr/bin/env python3
"""Participant-representative sensitivity subset .

Group key: (study_namespace, participant_id) only when participant_status ==
'identifiable' and participant_id non-missing; every other record is its own
accession-level unit. Representative ranking: completeness desc, contamination
asc, N50 desc, n_contigs asc, accession lexicographic. Compare against the
packet's expected selection table; any mismatch is reported, not silently
accepted.
"""
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import RESULTS, PACKET, read_tsv, write_tsv

OUT = os.path.join(RESULTS, "08_defense_mge_vfdb")


def main():
    meta = read_tsv(os.path.join(PACKET, "00_input/cohort/main152_metadata.tsv"))
    c2 = {r["accession"]: r for r in read_tsv(os.path.join(PACKET, "00_input/reference_sources/checkm2_main152.tsv"))}
    am = {r["accession"]: r for r in read_tsv(os.path.join(RESULTS, "00_preflight/assembly_metrics.tsv"))}

    groups = {}
    for r in meta:
        acc = r["strain"]
        if r["participant_status"] == "identifiable" and r["participant_id"] not in ("", "unknown", "not applicable", "NA"):
            key = json.dumps(["participant", r["study_namespace"], r["participant_id"]])
        else:
            key = json.dumps(["accession", acc])
        groups.setdefault(key, []).append(acc)

    rows = []
    for key, accs in sorted(groups.items()):
        def rank(a):
            return (-float(c2[a]["Completeness"]), float(c2[a]["Contamination"]),
                    -int(am[a]["n50_bp"]), int(am[a]["n_contigs"]), a)
        ordered = sorted(accs, key=rank)
        for i, a in enumerate(ordered):
            m = next(r for r in meta if r["strain"] == a)
            rows.append(dict(
                accession=a, participant_key=key, species_group=m["group_for_testing"],
                rank_in_participant=i + 1, selected=1 if i == 0 else 0,
                reason="highest prespecified quality rank" if i == 0 else
                       "lower prespecified quality rank than selected representative",
                n_in_participant=len(accs)))
    write_tsv(os.path.join(OUT, "participant_subset_manifest.tsv"),
              ["accession", "participant_key", "species_group", "rank_in_participant",
               "selected", "reason", "n_in_participant"], rows)

    sel = [r for r in rows if r["selected"] == 1]
    gc = {}
    for r in sel:
        gc[r["species_group"]] = gc.get(r["species_group"], 0) + 1

    exp = {r["accession"]: r for r in read_tsv(os.path.join(PACKET, "00_input/reference_sources/participant_sensitivity_selection_expected.tsv"))}
    mism = []
    for r in rows:
        e = exp.get(r["accession"])
        if e is None or str(r["selected"]) != str(e["selected"]):
            mism.append((r["accession"], r["selected"], e["selected"] if e else "ABSENT"))
    print(f"selected n={len(sel)} groups={gc}")
    print(f"expected-table mismatches: {len(mism)}")
    for m in mism[:20]:
        print("  MISMATCH", m)
    with open(os.path.join(OUT, "participant_subset_selection_check.json"), "w") as fh:
        json.dump(dict(n_selected=len(sel), group_counts=gc,
                       expected_n=121, expected_groups={"subflava": 24, "perflava": 90, "flavescens": 7},
                       mismatches=[list(m) for m in mism],
                       status="PASS" if not mism and len(sel) == 121 else "REVIEW"), fh, indent=1)
    sys.exit(0 if not mism and len(sel) == 121 else 1)


if __name__ == "__main__":
    main()
