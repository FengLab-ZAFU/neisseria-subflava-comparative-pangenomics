#!/usr/bin/env python3
"""Write 10_figures_source/Figure2_iTOL/reference/tip_metadata.tsv (contract F2_TIPS).

Columns: accession;display_label;species_group;collection_region;tetM_status;unit_id
Grain: 152 tip rows, accession-ordered by the frozen roster order.
Sources: 00_input/cohort/main152_display_metadata.tsv + 07_amr_tetm/amr_tetm_per_genome.tsv
"""
import csv, os, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import RESULTS, read_roster

META = os.path.join(RESULTS, "00_input/cohort/main152_display_metadata.tsv")
TETM = os.path.join(RESULTS, "07_amr_tetm/amr_tetm_per_genome.tsv")
OUT = os.path.join(RESULTS, "10_figures_source/Figure2_iTOL/reference/tip_metadata.tsv")

roster = read_roster()
meta = {}
with open(META, newline="") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        meta[r["accession"]] = r
tetm = {}
with open(TETM, newline="") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        tetm[r["accession"]] = r["tetM_status"]

COLS = ["accession", "display_label", "species_group", "collection_region",
        "tetM_status", "unit_id"]
with open(OUT, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=COLS, delimiter="\t", lineterminator="\n")
    w.writeheader()
    for acc in roster:
        m = meta[acc]
        w.writerow(dict(accession=acc, display_label=m["display_label"],
                        species_group=m["operational_label"],
                        collection_region=m["collection_region"],
                        tetM_status=tetm[acc], unit_id=m["unit_id"]))

# cross-check against iTOL annotation colors already installed
print(f"wrote {OUT} rows={len(roster)}")
from collections import Counter
print("tetM:", Counter(tetm[a] for a in roster))
print("species:", Counter(meta[a]["operational_label"] for a in roster))
