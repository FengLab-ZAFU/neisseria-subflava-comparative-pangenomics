#!/usr/bin/env python3
"""fix3-H2 (Phase A4 reference calibration): type-material-anchored re-screen
of ALL 1,303 screened records.

v4 anchor rule (task 2026-10-09 sec 2.5): membership decided ONLY by verified
type material -
  subflava    ATCC 49275      GCA_005221305.1 (complete chromosome, root ref)
  flavescens  ATCC 13120      GCA_005221285.1 (complete chromosome, root ref)
  perflava    CCUG 17915      two assemblies of the SAME type culture
                              (GCA_023472875.1 + GCA_042649725.1); per-species
                              channel = max mean ANI over the two - one culture
                              is NOT two independent species supports.
The historical non-type perflava reference (27098_8_142, GCA_041433205.1) is
a CONTROL only: it defines the expanded_reference_sensitivity channel and may
aid candidate discovery, but never decides main-cohort membership or labels.

Membership (unchanged project rule): bidirectional mean ANI >= 95 AND
min bidirectional AF >= 0.70 for at least one anchored species channel.
Missing results stay UNRESOLVED - never zero-filled.

Outputs (results_revision_20261009_fix3/01_reference_v4/):
  rescreen_v4.tsv                 per-ref bidirectional ANI/matched/total/AF
                                  for all 5 refs + channels + v4 membership/labels
  reference_rule_impact_v4.tsv    v4 vs v3 membership/label changes, incl.
                                  disposition for v3-primary records
  ani_unresolved_disposition.tsv  all 18 ANI_UNRESOLVED records, closed out
  reference_identity_table_v4.tsv ref identity + sha256 (incl. type_dup)
  ref_vs_ref_ani_v4.tsv           all ref pairs incl. type_dup
"""
import os
import csv
import hashlib
from collections import Counter
from pathlib import Path

ROOT = Path(os.environ["NEISSERIA_PROJECT_ROOT"])
OLD = ROOT / "results_revision_20261001" / "01_retrieval_taxonomy"
FIX = ROOT / "results_revision_20261006_fix"
F2 = ROOT / "results_revision_20261008_fix2"
OUT = ROOT / "results_revision_20261009_fix3" / "01_reference_v4"
RAW_T = FIX / "02_reference_check/fastani_raw"          # perflava_type
RAW_D = OUT / "fastani_raw_dup"                          # perflava_type_dup

ANI_MIN, AF_MIN, MARGIN_TEST = 95.0, 0.70, 0.10

def parse_fa(p):
    p = Path(p)
    if not p.exists() or p.stat().st_size == 0:
        return None
    f = p.read_text().strip().split("\t")
    if len(f) < 5:
        return None
    try:
        return float(f[2]), int(f[3]), int(f[4])
    except ValueError:
        return None

def num(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None

# ---------- reference identity + hashes ----------
REFS = [
    ("subflava", "N_subflava_ATCC49275.fna", "Neisseria subflava", "ATCC 49275", "GCA_005221305.1",
     "TYPE_MATERIAL (NCBI); type strain of N. subflava", "v4 MAIN ANCHOR"),
    ("flavescens", "N_flavescens_ATCC13120.fna", "Neisseria flavescens", "ATCC 13120", "GCA_005221285.1",
     "TYPE_MATERIAL (NCBI); type strain of N. flavescens", "v4 MAIN ANCHOR"),
    ("perflava_type", "GCA_023472875.1.fna", "Neisseria perflava", "CCUG 17915", "GCA_023472875.1",
     "TYPE_MATERIAL (NCBI); type strain of N. perflava", "v4 MAIN ANCHOR"),
    ("perflava_type_dup", "GCA_042649725.1.fna", "Neisseria perflava", "CCUG 17915", "GCA_042649725.1",
     "TYPE_MATERIAL (NCBI); same type culture as perflava_type (second assembly)", "v4 MAIN ANCHOR (same culture, not independent support)"),
    ("perflava_historical", "N_perflava_27098_8_142_historical.fna", "Neisseria perflava", "27098_8_142", "GCA_041433205.1",
     "NOT type material; historical reference", "CONTROL ONLY - excluded from main anchoring; defines expanded sensitivity"),
]
REFDIR = FIX / "02_reference_check/refs"
with open(OUT / "reference_identity_table_v4.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["ref_id", "file", "organism", "strain", "assembly_accession",
                "type_status_evidence", "v4_role", "fasta_sha256"])
    for rid, fn, org, strain, acc, tm, role in REFS:
        h = hashlib.sha256((REFDIR / fn).read_bytes()).hexdigest()
        w.writerow([rid, fn, org, strain, acc, tm, role, h])

# ---------- ref-vs-ref ANI (v2 table + new type_dup pairs) ----------
rvr = list(csv.DictReader(open(FIX / "02_reference_check/ref_vs_ref_ani.tsv"), delimiter="\t"))
have_pairs = {(r["query_ref"], r["target_ref"]) for r in rvr}
for a in ("perflava_type_dup",):
    for b in ("subflava", "flavescens", "perflava_type", "perflava_historical"):
        for q, t in ((a, b), (b, a)):
            if (q, t) in have_pairs:
                continue
            r = parse_fa(OUT / "ref_vs_ref" / f"{q}_vs_{t}.txt")
            rvr.append({"query_ref": q, "target_ref": t,
                        "ani": r[0] if r else "", "matched_frags": r[1] if r else "",
                        "total_frags": r[2] if r else "",
                        "af": round(r[1] / r[2], 4) if r else "NO_RESULT"})
with open(OUT / "ref_vs_ref_ani_v4.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, delimiter="\t", fieldnames=list(rvr[0].keys()))
    w.writeheader()
    w.writerows(rvr)

# ---------- per-record channels ----------
old = {r["accession"]: r for r in csv.DictReader(open(OLD / "ani_screen_results.tsv"), delimiter="\t")}

def channel_from_old(o, sp):
    """subflava / flavescens / perflava(historical) channels from the frozen screen."""
    g2r = (num(o.get(f"ani_{sp}_genome2ref")), num(o.get(f"match_{sp}_g2r")), num(o.get(f"total_{sp}_g2r")))
    r2q = (num(o.get(f"ani_{sp}_ref2genome")), num(o.get(f"match_{sp}_r2q")), num(o.get(f"total_{sp}_r2q")))
    if None in (*g2r, *r2q) or g2r[2] == 0 or r2q[2] == 0:
        return None
    return {"ani_g2r": g2r[0], "match_g2r": int(g2r[1]), "total_g2r": int(g2r[2]),
            "af_g2r": g2r[1] / g2r[2], "ani_r2q": r2q[0], "match_r2q": int(r2q[1]),
            "total_r2q": int(r2q[2]), "af_r2q": r2q[1] / r2q[2],
            "mean_ani": (g2r[0] + r2q[0]) / 2, "min_af": min(g2r[1] / g2r[2], r2q[1] / r2q[2])}

def channel_from_raw(acc, refname, rawdir):
    g2r = parse_fa(rawdir / f"{acc}_vs_{refname}.txt")
    r2q = parse_fa(rawdir / f"{refname}_vs_{acc}.txt")
    if not g2r or not r2q or g2r[2] == 0 or r2q[2] == 0:
        return None
    return {"ani_g2r": g2r[0], "match_g2r": g2r[1], "total_g2r": g2r[2], "af_g2r": g2r[1] / g2r[2],
            "ani_r2q": r2q[0], "match_r2q": r2q[1], "total_r2q": r2q[2], "af_r2q": r2q[1] / r2q[2],
            "mean_ani": (g2r[0] + r2q[0]) / 2, "min_af": min(g2r[1] / g2r[2], r2q[1] / r2q[2])}

def passes(c):
    return c is not None and c["mean_ani"] >= ANI_MIN and c["min_af"] >= AF_MIN

rows = []
for acc, o in sorted(old.items()):
    refs = {"subflava": channel_from_old(o, "subflava"),
            "flavescens": channel_from_old(o, "flavescens"),
            "perflava_type": channel_from_raw(acc, "perflava_type", RAW_T),
            "perflava_type_dup": channel_from_raw(acc, "perflava_type_dup", RAW_D),
            "perflava_historical": channel_from_old(o, "perflava")}
    # v4 anchored species channels (perflava = best of the two same-culture type assemblies)
    chan = {}
    for sp in ("subflava", "flavescens"):
        if refs[sp]:
            chan[sp] = refs[sp]
    pcands = [c for c in (refs["perflava_type"], refs["perflava_type_dup"]) if c]
    if pcands:
        chan["perflava"] = max(pcands, key=lambda c: c["mean_ani"])
    members_v4 = [s for s, c in chan.items() if passes(c)]
    # expanded control channel (historical non-type reference)
    exp_members = list(members_v4)
    if passes(refs["perflava_historical"]) and "perflava" not in exp_members:
        exp_members.append("perflava")
    member_v4 = len(members_v4) >= 1
    member_exp = len(exp_members) >= 1

    def label(channels, members):
        if not channels:
            return "ANI_UNRESOLVED", "", "", "NA", "ANI_UNRESOLVED"
        ranked = sorted(channels.items(), key=lambda kv: -kv[1]["mean_ani"])
        best_sp, best = ranked[0]
        second = ranked[1] if len(ranked) > 1 else (None, None)
        margin = round(best["mean_ani"] - second[1]["mean_ani"], 4) if second[1] else "NA"
        gft = best_sp if (margin != "NA" and margin >= MARGIN_TEST and len(members) <= 1) else "ambiguous"
        return best_sp, round(best["mean_ani"], 4), (second[0] or ""), margin, gft

    b, ba, s, mg, gft = label(chan, members_v4)
    chan_exp = dict(chan)
    if refs["perflava_historical"] and "perflava" not in chan_exp:
        chan_exp["perflava"] = refs["perflava_historical"]
    elif refs["perflava_historical"]:
        chan_exp["perflava"] = max([chan_exp["perflava"], refs["perflava_historical"]],
                                   key=lambda c: c["mean_ani"])
    eb, eba, es, emg, egft = label(chan_exp, exp_members)

    row = {"accession": acc}
    for rid in ("subflava", "flavescens", "perflava_type", "perflava_type_dup", "perflava_historical"):
        c = refs[rid]
        for k in ("ani_g2r", "match_g2r", "total_g2r", "af_g2r", "ani_r2q", "match_r2q",
                  "total_r2q", "af_r2q", "mean_ani", "min_af"):
            v = c[k] if c else ""
            row[f"{rid}_{k}"] = round(v, 4) if isinstance(v, float) else v
    row.update({
        "member_species_v4": ";".join(members_v4),
        "associated_group_member_v4": member_v4,
        "best_species_v4": b, "best_mean_ani_v4": ba, "second_species_v4": s,
        "margin_v4": mg, "group_for_testing_v4": gft,
        "member_species_expanded": ";".join(exp_members),
        "associated_group_member_expanded": member_exp,
        "best_species_expanded": eb, "group_for_testing_expanded": egft,
        "historical_only_membership": (not member_v4) and member_exp,
    })
    rows.append(row)

with open(OUT / "rescreen_v4.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, delimiter="\t", fieldnames=list(rows[0].keys()))
    w.writeheader()
    w.writerows(rows)

# ---------- impact vs v3 ----------
v3 = {r["accession"]: r for r in csv.DictReader(open(FIX / "02_reference_check/rescreen_4ref.tsv"), delimiter="\t")}
v3prim = {r["accession"] for r in csv.DictReader(open(F2 / "03_cohorts/cohort_revision_primary_v3.tsv"), delimiter="\t")}
imp = []
for r in rows:
    acc = r["accession"]
    o = v3.get(acc, {})
    old_m = o.get("associated_group_member_new") == "True"
    new_m = r["associated_group_member_v4"] is True
    old_g, new_g = o.get("group_for_testing_new"), r["group_for_testing_v4"]
    if old_m != new_m or (new_m and old_g != new_g):
        imp.append({"accession": acc,
                    "v3_member": old_m, "v4_member": new_m,
                    "v3_group": old_g, "v4_group": new_g,
                    "change": ("membership_loss" if old_m and not new_m else
                               "membership_gain" if new_m and not old_m else
                               f"label_change:{old_g}->{new_g}"),
                    "in_v3_primary_199": acc in v3prim,
                    "historical_only_membership": r["historical_only_membership"],
                    "subflava_mean_ani": r["subflava_mean_ani"], "subflava_min_af": r["subflava_min_af"],
                    "flavescens_mean_ani": r["flavescens_mean_ani"], "flavescens_min_af": r["flavescens_min_af"],
                    "perflava_type_mean_ani": r["perflava_type_mean_ani"], "perflava_type_min_af": r["perflava_type_min_af"],
                    "perflava_type_dup_mean_ani": r["perflava_type_dup_mean_ani"], "perflava_type_dup_min_af": r["perflava_type_dup_min_af"],
                    "perflava_historical_mean_ani": r["perflava_historical_mean_ani"],
                    "perflava_historical_min_af": r["perflava_historical_min_af"]})
with open(OUT / "reference_rule_impact_v4.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, delimiter="\t", fieldnames=list(imp[0].keys()) if imp else ["accession"])
    w.writeheader()
    w.writerows(imp)

# ---------- ANI_UNRESOLVED close-out ----------
mag = {r["accession"]: r for r in csv.DictReader(open(FIX / "01_mag_audit/mag_source_audit.tsv"), delimiter="\t")}
chk = {r["accession"] for r in csv.DictReader(open(OLD / "checkm2_merged.tsv"), delimiter="\t")}
disp = []
for r in rows:
    if r["best_species_v4"] != "ANI_UNRESOLVED":
        continue
    acc = r["accession"]
    m = mag.get(acc, {})
    if m.get("decision") == "EXCLUDE_confirmed_MAG":
        d = ("excluded as confirmed MAG (independent of ANI): " + m.get("decision_reason", ""))
        cat = "MAG_excluded"
    elif acc == "GCA_051132615.1":
        d = ("NOT a MAG; FASTA present (3.24 MB, 'Complete Genome' Ec49-e6-T10). fastANI 1.34 produced "
             "EMPTY output in BOTH directions vs all 5 references (ran to completion, no error); skani 0.2.2 "
             "second method: no detectable ANI vs any of the 5 references (positive control on a known member "
             "returns 94.94). ANI is below the ~80% detection floor of both methods, i.e. far under the 95% "
             "membership threshold -> confirmed NON-member. CheckM2 row absent -> QC would also be unevaluable. "
             "No zero-fill: result recorded as evidence-based exclusion, not missing data.")
        cat = "confirmed_non_member_two_methods"
    else:
        d = "no MAG basis and no two-method confirmation -> REVIEW required"
        cat = "REVIEW"
    disp.append({"accession": acc, "category": cat, "mag_decision": m.get("decision", ""),
                 "checkm2_row_present": acc in chk, "disposition": d})
with open(OUT / "ani_unresolved_disposition.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, delimiter="\t", fieldnames=list(disp[0].keys()))
    w.writeheader()
    w.writerows(disp)

n_v4 = sum(1 for r in rows if r["associated_group_member_v4"] is True)
n_exp = sum(1 for r in rows if r["associated_group_member_expanded"] is True)
print(f"records: {len(rows)} | v4 anchored members: {n_v4} | expanded (incl historical): {n_exp}")
print(f"historical-only memberships: {sum(1 for r in rows if r['historical_only_membership'])}")
print(f"ANI_UNRESOLVED: {sum(1 for r in rows if r['best_species_v4']=='ANI_UNRESOLVED')} "
      f"(categories: {dict(Counter(d['category'] for d in disp))})")
print(f"impact rows: {len(imp)} | membership_loss: {sum(1 for i in imp if i['change']=='membership_loss')} "
      f"| of which in v3 primary 199: {sum(1 for i in imp if i['change']=='membership_loss' and i['in_v3_primary_199'])}")
print(f"label changes among v4 members: {sum(1 for i in imp if i['change'].startswith('label_change'))} "
      f"| in v3 primary: {sum(1 for i in imp if i['change'].startswith('label_change') and i['in_v3_primary_199'])}")
