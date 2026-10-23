#!/usr/bin/env python3
"""A2: parse bidirectional FastANI results and apply associated-group membership rules.

Output: 01_retrieval_taxonomy/ani_screen_results.tsv
Columns per genome x reference:
  per-direction ANI, matched/total frags (both directions), bidirectional mean ANI,
  min bidirectional AF, membership pass (mean ANI>=95 & min AF>=0.70, both dirs present)
Then per genome: best reference, margin (vs 2nd best valid bidirectional estimate),
operational_group, group_for_testing (margin>=0.10 & unique member ref else ambiguous).
"""
import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent.parent
OUT = ROOT / "results_revision_20261001" / "01_retrieval_taxonomy"
RAW = OUT / "fastani_raw"
REFS = ["subflava", "perflava", "flavescens"]

ANI_MIN = 95.0
AF_MIN = 0.70
MARGIN_TEST = 0.10


def parse_fastaani(path):
    """Return (ani, matched, total) or None if file missing/empty."""
    p = Path(path)
    if not p.exists() or p.stat().st_size == 0:
        return None
    line = p.read_text().strip()
    if not line:
        return None
    f = line.split("\t")
    # fastANI: query reference ANI fragments_matching total_query_fragments
    if len(f) < 5:
        return None
    try:
        return float(f[2]), int(f[3]), int(f[4])
    except ValueError:
        return None


genomes = sorted(p.stem for p in (OUT / "fasta").glob("*.fna"))
rows = []
for acc in genomes:
    rec = {"accession": acc}
    per_ref = []
    for ref in REFS:
        q2r = parse_fastaani(RAW / f"{acc}_vs_{ref}.txt")     # genome query -> ref
        r2q = parse_fastaani(RAW / f"{ref}_vs_{acc}.txt")     # ref query -> genome
        rec[f"ani_{ref}_genome2ref"], rec[f"match_{ref}_g2r"], rec[f"total_{ref}_g2r"] = q2r if q2r else ("", "", "")
        rec[f"ani_{ref}_ref2genome"], rec[f"match_{ref}_r2q"], rec[f"total_{ref}_r2q"] = r2q if r2q else ("", "", "")
        if q2r and r2q and q2r[2] > 0 and r2q[2] > 0:
            mean_ani = (q2r[0] + r2q[0]) / 2
            af_g2r = q2r[1] / q2r[2]
            af_r2q = r2q[1] / r2q[2]
            min_af = min(af_g2r, af_r2q)
            member = mean_ani >= ANI_MIN and min_af >= AF_MIN
            per_ref.append({"ref": ref, "mean_ani": mean_ani, "min_af": min_af, "member": member})
            rec[f"mean_ani_{ref}"] = round(mean_ani, 4)
            rec[f"af_{ref}_genome2ref"] = round(af_g2r, 4)
            rec[f"af_{ref}_ref2genome"] = round(af_r2q, 4)
            rec[f"min_af_{ref}"] = round(min_af, 4)
            rec[f"member_{ref}"] = member
        else:
            rec[f"mean_ani_{ref}"] = ""
            rec[f"af_{ref}_genome2ref"] = ""
            rec[f"af_{ref}_ref2genome"] = ""
            rec[f"min_af_{ref}"] = ""
            rec[f"member_{ref}"] = "ANI_UNRESOLVED"
    valid = [p for p in per_ref]
    members = [p for p in valid if p["member"]]
    valid_est = [p for p in valid]  # all with bidirectional estimates
    if valid_est:
        ranked = sorted(valid_est, key=lambda x: -x["mean_ani"])
        best = ranked[0]
        rec["best_ref"] = best["ref"]
        rec["best_mean_ani"] = round(best["mean_ani"], 4)
        if len(ranked) > 1:
            rec["second_ref"] = ranked[1]["ref"]
            rec["second_mean_ani"] = round(ranked[1]["mean_ani"], 4)
            rec["margin"] = round(best["mean_ani"] - ranked[1]["mean_ani"], 4)
        else:
            rec["second_ref"] = ""; rec["second_mean_ani"] = ""; rec["margin"] = "NA"
        rec["operational_group"] = best["ref"]
        # group_for_testing: margin>=0.10 AND not two refs both passing membership
        n_members = len(members)
        if rec["margin"] != "NA" and rec["margin"] >= MARGIN_TEST and n_members <= 1:
            rec["group_for_testing"] = best["ref"]
        else:
            rec["group_for_testing"] = "ambiguous"
        rec["associated_group_member"] = n_members >= 1
    else:
        for k in ["best_ref", "best_mean_ani", "second_ref", "second_mean_ani", "margin",
                  "operational_group", "group_for_testing"]:
            rec[k] = "ANI_UNRESOLVED" if k != "margin" else "NA"
        rec["associated_group_member"] = "ANI_UNRESOLVED"
    rows.append(rec)

fields = ["accession"]
for ref in REFS:
    fields += [f"ani_{ref}_genome2ref", f"match_{ref}_g2r", f"total_{ref}_g2r",
               f"ani_{ref}_ref2genome", f"match_{ref}_r2q", f"total_{ref}_r2q",
               f"mean_ani_{ref}", f"af_{ref}_genome2ref", f"af_{ref}_ref2genome",
               f"min_af_{ref}", f"member_{ref}"]
fields += ["best_ref", "best_mean_ani", "second_ref", "second_mean_ani", "margin",
           "operational_group", "group_for_testing", "associated_group_member"]

with open(OUT / "ani_screen_results.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, delimiter="\t", fieldnames=fields, extrasaction="ignore")
    w.writeheader()
    w.writerows(rows)

# summary
n = len(rows)
n_member = sum(1 for r in rows if r["associated_group_member"] is True)
n_unres = sum(1 for r in rows if r["associated_group_member"] == "ANI_UNRESOLVED")
margins = [r["margin"] for r in rows if isinstance(r["margin"], float)]
import statistics
print(f"genomes screened: {n}")
print(f"associated-group members: {n_member} | non-members: {n - n_member - n_unres} | ANI_UNRESOLVED: {n_unres}")
print(f"margin <0.10: {sum(1 for m in margins if m < 0.10)} | margin <0.50: {sum(1 for m in margins if m < 0.50)} | margin NA: {n_unres}")
