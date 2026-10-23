#!/usr/bin/env python3
"""M9 (iTOL part): assemble the Figure2_iTOL directory.

- copies the four packet tip annotation templates after verifying each has
  exactly the 152 roster keys
- 04_tetM_carriage.txt is REGENERATED from the M6 verified per-genome status
  (4-state: positive/indeterminate/missing/negative; indeterminate=#E7B95F,
  missing=#B8B8B8) — never white-for-nonnegative
- writes local_install_map.tsv and settings note
"""
import json
import os
import shutil
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import RESULTS, PACKET, read_roster, read_tsv

FIG = os.path.join(RESULTS, "10_figures_source/Figure2_iTOL")
TPL = os.path.join(PACKET, "00_input/iTOL_template/annotations")
ANN = os.path.join(FIG, "annotations")

COLOR = {"positive_complete": "#000000", "negative": "#FFFFFF",
         "indeterminate": "#E7B95F", "missing": "#B8B8B8"}
LABEL = {"positive_complete": "tet(M)+", "negative": "tet(M)-",
         "indeterminate": "tet(M) indeterminate", "missing": "tet(M) missing"}


def tpl_keys(fn):
    keys = []
    for ln in open(os.path.join(TPL, fn)):
        if ln.startswith("GCA_"):
            keys.append(ln.split("\t")[0])
    return keys


def main():
    roster = read_roster()
    os.makedirs(ANN, exist_ok=True)
    os.makedirs(os.path.join(FIG, "reference"), exist_ok=True)
    problems = []

    for fn in ["01_tip_labels.txt", "02_operational_ANI_assignment.txt", "03_collection_region.txt"]:
        keys = tpl_keys(fn)
        if not (len(keys) == 152 and len(set(keys)) == 152 and set(keys) == set(roster)):
            problems.append(f"{fn}: key set invalid ({len(keys)})")
            continue
        shutil.copy(os.path.join(TPL, fn), os.path.join(ANN, fn))

    # tetM ring from verified M6 status
    m6 = {r["accession"]: r for r in read_tsv(os.path.join(RESULTS, "07_amr_tetm/amr_tetm_per_genome.tsv"))}
    with open(os.path.join(ANN, "04_tetM_carriage.txt"), "w") as fh:
        fh.write("DATASET_COLORSTRIP\nSEPARATOR TAB\nDATASET_LABEL\ttet(M) carriage\n")
        fh.write("COLOR\t#555555\nSTRIP_WIDTH\t12\n")
        fh.write("SHOW_INTERNAL\t0\n")
        fh.write("LEGEND_TITLE\ttet(M) status\n")
        fh.write("LEGEND_SHAPES\t1\t1\t1\t1\n")
        fh.write("LEGEND_COLORS\t#000000\t#FFFFFF\t#E7B95F\t#B8B8B8\n")
        fh.write("LEGEND_LABELS\ttet(M)+\ttet(M)-\ttet(M) indeterminate\ttet(M) missing\n")
        fh.write("DATA\n")
        for acc in roster:
            st = m6[acc]["tetM_status"]
            if st not in COLOR:
                problems.append(f"{acc}: unexpected tetM_status {st}")
                continue
            fh.write(f"{acc}\t{COLOR[st]}\t{LABEL[st]}\n")

    # settings / readme
    with open(os.path.join(FIG, "SETTINGS_NOTES.md"), "w") as fh:
        fh.write("""# Figure2 iTOL settings (server-delivered text inputs; drawing is local Windows)

- Tree: single-copy core-gene concatenation, maximum likelihood (IQ-TREE 2.4.0,
  GTR+F+I+R4, UFBoot 1000 + SH-aLRT 1000), midpoint-rooted for visualization only.
  No recombination correction; no outgroup; root position carries no ancestral
  or directional interpretation.
- Upload tree: Figure2_main152_midpoint_iTOL.nwk (root of this directory; unique
  internal node IDs N0001...; tip IDs = assembly accessions; branch lengths kept).
- reference/: science record trees (unrooted ML, midpoint display with raw
  SH-aLRT/UFBoot labels), branch_support.tsv, itol_node_id_map.tsv.
- annotations/: 01 tip labels, 02 reference-ANI operational assignment ring,
  03 collection region ring, 04 tet(M) 4-state ring, 05 branch support symbols
  (DATASET_SYMBOL; only edges with SH-aLRT>=80 AND UFBoot>=95; internal node
  name display must be OFF so N IDs are not shown).
- iTOL view settings to reproduce historical layout: circular mode, branch
  lengths visible, rotation 350 degrees, existing ring colors from the templates.
- Ring order inside-out: 02 operational ANI assignment, 03 collection region,
  04 tet(M).
""")

    # local install map
    with open(os.path.join(FIG, "local_install_map.tsv"), "w", newline="") as fh:
        fh.write("server_relative_path\tlocal_target\n")
        fh.write("10_figures_source/Figure2_iTOL\t"
                 "D:\\Work\\Neisseria_subflava\\manuscript\\Figures\\new\\Figure2_iTOL\n")
        fh.write("# note\tLinux cannot create a literal D:\\ path; local install must archive "
                 "previous preview files, verify hashes, then update this single directory.\n")

    print("iTOL assembly:", "PASS" if not problems else problems)
    sys.exit(0 if not problems else 1)


if __name__ == "__main__":
    main()
