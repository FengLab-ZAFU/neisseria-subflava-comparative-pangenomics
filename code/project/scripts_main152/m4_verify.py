#!/usr/bin/env python3
"""M4 independent verification: recompute key M4 results from the delivered
files with a code path separate from the producers, and check the plan's
acceptance criteria. Exits non-zero on any failed check."""
import csv
import os
import sys
from collections import Counter

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import RESULTS, read_roster, read_tsv

S04 = os.path.join(RESULTS, "04_species_pangenome")
S05 = os.path.join(RESULTS, "05_structure")
PAV = os.path.join(RESULTS, "02_pangenome/pgap2_main152/pav_tables")

fails = []


def check(name, ok, detail=""):
    print(f"{'PASS' if ok else 'FAIL'}\t{name}\t{detail}")
    if not ok:
        fails.append(name)


roster = read_roster()
meta = {r["strain"]: r for r in read_tsv(
    os.path.join(RESULTS, "00_freeze/input_packet_snapshot/00_input/cohort/main152_metadata.tsv"))}
group_of = {a: meta[a]["group_for_testing"] for a in roster}

# raw PAV (independent reload)
with open(os.path.join(PAV, "pav_binary.tsv"), newline="", encoding="utf-8") as fh:
    r = csv.reader(fh, delimiter="\t")
    header = next(r)
    cols = header[1:]
    clusters, rows = [], []
    for row in r:
        clusters.append(row[0])
        rows.append([int(x) for x in row[1:]])
M = np.asarray(rows, dtype=np.uint8)
check("pav shape 13215x152, cols == roster", M.shape == (13215, 152) and cols == roster)

# 1. prevalence table
prev = read_tsv(os.path.join(S04, "cluster_species_prevalence.tsv"))
check("prevalence rows == 13215", len(prev) == 13215, str(len(prev)))
check("prevalence cluster set == dictionary", {r["cluster"] for r in prev} == set(clusters))
ok = all(int(r["subflava_total"]) == 31 and int(r["perflava_total"]) == 114
         and int(r["flavescens_total"]) == 7 for r in prev)
check("prevalence totals 31/114/7 on all rows", ok)
gi = {g: np.array([i for i, a in enumerate(roster) if group_of[a] == g]) for g in
      ("subflava", "perflava", "flavescens")}
pres = {g: M[:, gi[g]].sum(axis=1) for g in gi}
cidx = {c: i for i, c in enumerate(clusters)}
bad = sum(1 for r in prev for g in gi
          if int(r[f"{g}_present"]) != int(pres[g][cidx[r["cluster"]]]))
check("prevalence counts match independent PAV recount", bad == 0, f"{bad} mismatches")

# 2. sharing closure recomputed from the prevalence file
pat = Counter()
for r in prev:
    key = tuple(g for g in ("subflava", "perflava", "flavescens") if int(r[f"{g}_present"]) > 0)
    pat[key] += 1
name_map = {("subflava",): "subflava_only", ("perflava",): "perflava_only",
            ("flavescens",): "flavescens_only",
            ("subflava", "perflava"): "subflava_perflava",
            ("subflava", "flavescens"): "subflava_flavescens",
            ("perflava", "flavescens"): "perflava_flavescens",
            ("subflava", "perflava", "flavescens"): "all_three"}
sh = {r["pattern"]: int(r["n_clusters"]) for r in read_tsv(os.path.join(S04, "sharing_patterns.tsv"))}
check("sharing has 7 patterns, sum 13215", len(sh) == 7 and sum(sh.values()) == 13215)
rec = {name_map[k]: v for k, v in pat.items()}
check("sharing counts match prevalence-file recount", all(sh.get(k, 0) == rec.get(k, 0) for k in sh),
      str(sh))

# 3. group_conserved recomputed
gc = {r["group"]: r for r in read_tsv(os.path.join(S04, "group_conserved.tsv"))}
for g in gi:
    n = len(gi[g])
    others = [o for o in gi if o != g]
    allp = int((pres[g] == n).sum())
    ge95 = int((100 * pres[g] >= 95 * n).sum())
    lim = int(((pres[g] > 0) & (pres[others[0]] == 0) & (pres[others[1]] == 0)).sum())
    r = gc[g]
    check(f"group_conserved {g}", int(r["n_clusters_all_present"]) == allp
          and int(r["n_clusters_ge95pct"]) == ge95
          and int(r["n_clusters_group_limited"]) == lim
          and int(r["ge95pct_min_present"]) * 100 >= 95 * n
          and (int(r["ge95pct_min_present"]) - 1) * 100 < 95 * n,
          f"all={allp} ge95={ge95} min_ge95={r['ge95pct_min_present']} limited={lim}")
check("flavescens ge95 denominator is 7/7", gc["flavescens"]["ge95pct_min_present"] == "7")

# 4. burden-by-group recomputed
burden = {r["accession"]: r for r in read_tsv(os.path.join(PAV, "per_genome_gene_burden.tsv"))}
keymap = {"n_CDS_copies": "n_CDS_copies", "n_clusters": "n_clusters_present",
          "n_accessory_clusters": "n_accessory_clusters", "n_singletons": "n_singleton_clusters"}
bb = read_tsv(os.path.join(S04, "per_genome_burden_by_group.tsv"))
check("burden-by-group rows == 12", len(bb) == 12)
bad = 0
for r in bb:
    vals = np.array([float(burden[roster[i]][keymap[r["metric"]]]) for i in gi[r["group"]]])
    if not (abs(float(r["mean"]) - vals.mean()) < 1e-5
            and abs(float(r["sd"]) - vals.std(ddof=1)) < 1e-5
            and abs(float(r["median"]) - np.median(vals)) < 1e-6
            and int(r["min"]) == vals.min() and int(r["max"]) == vals.max()
            and int(r["n"]) == len(vals)):
        bad += 1
check("burden-by-group stats match recomputation", bad == 0, f"{bad} bad rows")

# contract per_genome_gene_burden
pgb = read_tsv(os.path.join(S04, "per_genome_gene_burden.tsv"))
check("F3_GENOMES 152 rows, keys = roster", len(pgb) == 152 and [r["accession"] for r in pgb] == roster)

# partition summary
psum = read_tsv(os.path.join(S04, "partition_summary.tsv"))
check("partition summary 5 rows, total 13215", len(psum) == 5
      and sum(int(r["n_clusters"]) for r in psum) == 13215)
n_pres = M.sum(axis=1)
check("singletons in Cloud == recount",
      int([r for r in psum if r["partition"] == "Cloud"][0]["n_singletons"])
      == int((n_pres == 1).sum()), str(int((n_pres == 1).sum())))

# 5. accumulation raw/summary
raw = read_tsv(os.path.join(S04, "accumulation_raw.tsv"))
check("accumulation raw rows == 304000", len(raw) == 304000, str(len(raw)))
sizes = {"pooled152": 152, "subflava": 31, "perflava": 114, "flavescens": 7}
cnt = Counter((r["group"], r["iteration"]) for r in raw)
check("raw complete grid", all(cnt[(g, str(i))] == sizes[g] for g in sizes for i in range(1, 1001)))
seeds = {(r["group"], r["iteration"]): int(r["seed"]) for r in raw if r["n"] == "1"}
expect = {"pooled152": 0, "subflava": 1, "perflava": 2, "flavescens": 3}
check("raw seeds match deterministic scheme",
      all(v == 20261015000000 + expect[g] * 100000 + int(i) for (g, i), v in seeds.items()))
summ = {(r["group"], int(r["n"]), k): r for r in read_tsv(os.path.join(S04, "accumulation_summary.tsv"))
        for k in ("pan", "core")}
store = {}
for r in raw:
    store.setdefault((r["group"], int(r["n"])), []).append((int(r["pan_clusters"]), int(r["core_clusters"])))
bad = 0
for (g, n), vals in store.items():
    for k, j in (("pan", 0), ("core", 1)):
        v = np.array([x[j] for x in vals], dtype=float)
        r = summ[(g, n, k)]
        if not (abs(float(r[f"{k}_mean"]) - v.mean()) < 1e-5
                and abs(float(r[f"{k}_median"]) - np.median(v)) < 1e-6
                and abs(float(r[f"{k}_q025"]) - np.percentile(v, 2.5)) < 1e-2
                and abs(float(r[f"{k}_q975"]) - np.percentile(v, 97.5)) < 1e-2):
            bad += 1
check("summary recomputed from raw matches file", bad == 0, f"{bad} bad cells")
# endpoints from raw vs independent matrix
for g in sizes:
    idx = gi[g] if g in gi else np.arange(152)
    union = int(M[:, idx].any(axis=1).sum()); inter = int(M[:, idx].all(axis=1).sum())
    endp = {x[0] for x in store[(g, sizes[g])]}
    endc = {x[1] for x in store[(g, sizes[g])]}
    n1 = np.array([x[0] for x in store[(g, 1)]], dtype=float)
    pg = M.sum(axis=0)[idx].astype(float)
    check(f"accum endpoints {g}", endp == {union} and endc == {inter}
          and abs(n1.mean() - pg.mean()) <= 5 * pg.std(ddof=1) / 10
          and set(n1.astype(int)) <= set(M.sum(axis=0)[idx].tolist()),
          f"union={union} inter={inter} n1mean={n1.mean():.3f} expect={pg.mean():.3f}")

# matched sampling
ms = read_tsv(os.path.join(S04, "matched_sampling_summary.tsv"))
n37 = sum(1 for r in ms if r["comparison"] == "three_groups_n1_7")
n231 = sum(1 for r in ms if r["comparison"] == "subflava_perflava_n1_31")
check("matched rows 42 + 124", n37 == 42 and n231 == 124, f"{n37}/{n231}")
bad = 0
for r in ms:
    v = np.array([x[0] if r["metric"] == "pan_clusters" else x[1]
                  for x in store[(r["group"], int(r["n"]))]], dtype=float)
    if abs(float(r["mean"]) - v.mean()) > 1e-5 or abs(float(r["median"]) - np.median(v)) > 1e-6:
        bad += 1
check("matched summary matches raw store", bad == 0)
eff = read_tsv(os.path.join(S04, "matched_sampling_effects.tsv"))
check("effects rows == 8 (3 pairs x2 metrics at n7 + 1 pair x2 at n31)", len(eff) == 8, str(len(eff)))

# 6. distance matrices + PCoA
for fs, feat_sel in (("shell_cloud", None), ("all_variable", None)):
    dm = {}
    with open(os.path.join(S05, f"jaccard_distance_{fs}.tsv"), newline="", encoding="utf-8") as fh:
        r = csv.reader(fh, delimiter="\t")
        hd = next(r)
        assert hd[0] == "accession" and hd[1:] == roster
        rows = []
        for row in r:
            rows.append([float(x) for x in row[1:]])
    D = np.array(rows)
    check(f"{fs} D symmetric, diag 0, shape", D.shape == (152, 152)
          and np.allclose(D, D.T, atol=1e-9) and np.all(np.diag(D) == 0))
    # independent Jaccard recomputation for a sample of pairs
    part = {r["cluster"]: r["partition_recomputed"] for r in read_tsv(os.path.join(PAV, "partition_table.tsv"))}
    if fs == "shell_cloud":
        sel = np.array([part[c] in ("Shell", "Cloud") for c in clusters])
    else:
        sel = (n_pres >= 1) & (n_pres <= 151)
    X = M[sel, :]
    rng = np.random.default_rng(1)
    ok = True
    for _ in range(500):
        i, j = rng.integers(0, 152, 2)
        a, b = X[:, i].astype(bool), X[:, j].astype(bool)
        u = (a | b).sum()
        d = 1.0 - (a & b).sum() / u if u > 0 else 0.0
        if abs(D[i, j] - d) > 5e-7:
            ok = False
            break
    check(f"{fs} 500 random Jaccard entries recomputed", ok)
    # independent PCoA recompute
    n = 152
    Jm = np.eye(n) - 1.0 / n
    B = -0.5 * Jm @ (D ** 2) @ Jm
    w = np.linalg.eigvalsh((B + B.T) / 2)[::-1]
    scale = np.abs(w).max()
    coords = read_tsv(os.path.join(S05, "pcoa_coordinates.tsv"))
    cc = [r for r in coords if r["feature_set"] == fs]
    check(f"{fs} pcoa rows 152, order matches roster", len(cc) == 152
          and [r["accession"] for r in cc] == roster)
    ling = w.min() < -1e-8 * scale and scale > 0
    if ling:
        D2 = D ** 2 - w.min()
        np.fill_diagonal(D2, 0.0)
        wr = np.linalg.eigvalsh((-0.5 * Jm @ D2 @ Jm))[::-1]
    else:
        wr = w
    pos = wr[wr > 0]
    fr1, fr2 = pos[0] / pos.sum(), pos[1] / pos.sum()
    check(f"{fs} axis fractions match independent recompute",
          abs(float(cc[0]["axis1_fraction"]) - fr1) < 1e-6
          and abs(float(cc[0]["axis2_fraction"]) - fr2) < 1e-6
          and cc[0]["correction"] == ("lingoes" if ling else "none"),
          f"axis1={fr1:.6f} axis2={fr2:.6f} correction={'lingoes' if ling else 'none'}")

# 7. distance statistics / permutation audit
ds = read_tsv(os.path.join(S05, "distance_statistics.tsv"))
check("distance_statistics 8 rows, all OK", len(ds) == 8 and all(r["status"] == "OK" for r in ds))
check("distance statistics B=999, p in (0,1]",
      all(int(r["B_used"]) == 999 and 0 < float(r["p"]) <= 1 for r in ds))
au = read_tsv(os.path.join(S05, "permutation_design_audit.tsv"))
tot = {r["permutation_design"]: r for r in au if r["block_id"] == "__TOTAL__"}
check("audit totals: both designs estimable, MC not enumeration",
      all(tot[d]["estimable"] == "yes" and tot[d]["enumeration_used"] == "no"
          and tot[d]["B_used"] == "999" for d in tot), str(len(tot)))
blocks = [r for r in au if r["permutation_design"] == "restricted_within_study_namespace"
          and r["block_id"] != "__TOTAL__"]
check("audit restricted blocks sum to 152", sum(int(r["block_size"]) for r in blocks) == 152,
      f"{len(blocks)} blocks")
# permutation label files: 999 unique, non-identity, block-consistent
for nm, design in (("free", None), ("restricted", None)):
    pl = read_tsv(os.path.join(S05, f"permutation_labels_{nm}.tsv"))
    cfgs = {tuple(r[a] for a in roster) for r in pl}
    ident = tuple(group_of[a] for a in roster)
    ok = len(pl) == 999 and len(cfgs) == 999 and ident not in cfgs
    if nm == "restricted":
        # block-consistency: multiset of labels within each namespace preserved
        ns = {a: meta[a]["study_namespace"] for a in roster}
        obs = Counter((ns[a], group_of[a]) for a in roster)
        ok = ok and all(Counter((ns[a], r[a]) for a in roster) == obs for r in pl)
    check(f"permutation_labels_{nm}: 999 unique non-identity" +
          (", within-block label multisets preserved" if nm == "restricted" else ""), ok)

# 8. encoding: UTF-8 no BOM, no CRLF
for root, _, files in os.walk(S04):
    for f in files:
        if f.endswith(".tsv"):
            b = open(os.path.join(root, f), "rb").read()
            if b.startswith(b"\xef\xbb\xbf") or b"\r" in b:
                fails.append(f"encoding {f}")
for root, _, files in os.walk(S05):
    for f in files:
        if f.endswith(".tsv"):
            b = open(os.path.join(root, f), "rb").read()
            if b.startswith(b"\xef\xbb\xbf") or b"\r" in b:
                fails.append(f"encoding {f}")
check("all TSV UTF-8 no BOM, LF endings", not any(x.startswith("encoding") for x in fails))

print()
if fails:
    print("FAILED CHECKS:", fails)
    sys.exit(1)
print("ALL M4 VERIFICATION CHECKS PASSED")
