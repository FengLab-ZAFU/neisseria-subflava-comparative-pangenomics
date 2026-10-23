#!/usr/bin/env python3
"""M3c: midpoint display rooting + iTOL text interface for the main152 ML tree.

Steps:
1. read IQ-TREE .treefile (internal labels "SH-aLRT/UFBoot"); save unrooted copy.
2. midpoint-root a copy for display only.
3. split_key binding: for each internal edge, tip sets on both sides sorted
   lexicographically; take the smaller side (ties: lexicographically smaller
   joined list); key = sha256(LF-joined tips + trailing LF).
4. invariance audit: same 152 tips; unrooted RF=0; pairwise patristic
   distances and total length within 1e-8; support mapping consistent.
5. derived iTOL upload tree: internal nodes renamed to unique N0001... ordered
   by sorted rooted-descendant-tip lists (root included); tips/lengths
   unchanged; itol_node_id_map.tsv maps node_id -> split_key/supports/role.
   Root gets no support; a supported nontrivial original edge split by the root
   is counted once (display side = the side containing the globally smallest
   tip; other side marked root_split_duplicate).

Run with: conda run -n gubbins python3 m3c_root_itol.py  (dendropy 5.0.11)
"""
import hashlib
import json
import os
import sys

import dendropy

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import RESULTS, read_roster, read_tsv, write_tsv

CORE = os.path.join(RESULTS, "03_phylogeny/strict_singlecopy_core")
IQ = os.path.join(CORE, "iqtree")
FIG = os.path.join(RESULTS, "10_figures_source/Figure2_iTOL")
TOL = 1e-8


def parse_label(node):
    """IQ-TREE internal label 'SHaLRT/UFBoot' -> (float, float) or None."""
    lb = node.label
    if not lb:
        return None
    parts = lb.split("/")
    try:
        if len(parts) == 2:
            return float(parts[0]), float(parts[1])
        return None, float(parts[0])
    except ValueError:
        return None


def split_key(tip_subset, all_tips):
    a = sorted(tip_subset)
    b = sorted(set(all_tips) - set(tip_subset))
    if not a or not b:
        return None
    small = a if (len(a), a) < (len(b), b) else b
    if len(a) == len(b):
        small = min(a, b)
    return hashlib.sha256(("\n".join(small) + "\n").encode()).hexdigest()


def edge_split_keys(tree, all_tips):
    """edge -> (split_key, side_tips_smaller) for internal edges; nontrivial only."""
    tree.update_bipartitions()
    out = {}
    taxa = tree.taxon_namespace
    for e in tree.postorder_edge_iter():
        if e.is_leaf():
            continue
        bp = e.bipartition
        side1 = sorted(t.label for t in bp.leafset_taxa(taxa))
        side2 = sorted(set(all_tips) - set(side1))
        if len(side1) < 2 or len(side2) < 2:
            continue  # trivial split (pendant edge): no internal support stats
        small = side1 if (len(side1), side1) <= (len(side2), side2) else side2
        key = hashlib.sha256(("\n".join(small) + "\n").encode()).hexdigest()
        out[key] = dict(edge=e, n_side1=len(side1), n_side2=len(side2))
    return out


def tip_labels(tree):
    return sorted(t.label for t in tree.taxon_namespace)


def patristic(tree):
    dm = tree.phylogenetic_distance_matrix()
    taxa = sorted(tree.taxon_namespace, key=lambda t: t.label)
    return {t.label: {t2.label: dm.distance(t, t2) for t2 in taxa} for t in taxa}


def main():
    roster = read_roster()
    treefile = os.path.join(IQ, "main152.treefile")
    if not os.path.isfile(treefile):
        sys.exit("FATAL: iqtree treefile missing")

    taxa = dendropy.TaxonNamespace()
    ml = dendropy.Tree.get(path=treefile, schema="newick", taxon_namespace=taxa,
                           preserve_underscores=True)
    tips = tip_labels(ml)
    assert set(tips) == set(roster), "tip set != roster"
    assert len(tips) == 152

    # 1. unrooted copy (as-is, with support labels)
    ml.write(path=os.path.join(CORE, "main152_ML_unrooted.nwk"), schema="newick",
             suppress_rooting=True, real_value_format_specifier=".17g")

    # support per edge via split keys (unrooted view)
    all_tips = tips
    edge_supports = {}  # split_key -> (sh, uf)
    edge_info = edge_split_keys(ml, all_tips)
    for key, info in edge_info.items():
        lab = parse_label(info["edge"].head_node)
        edge_supports[key] = lab

    total_len_unrooted = sum(e.length or 0.0 for e in ml.postorder_edge_iter())

    # 2. midpoint display tree
    mid = dendropy.Tree.get(path=treefile, schema="newick", taxon_namespace=taxa,
                            preserve_underscores=True)
    mid.reroot_at_midpoint(update_bipartitions=True, suppress_unifurcations=False)
    # re-bind internal labels from the UNROOTED split_key->support map:
    # reroot_at_midpoint exchanges the labels of two adjacent internal nodes
    # (found by independent audit G3.10), so labels traveling with nodes are
    # not trustworthy; the unrooted binding is authoritative. Both sides of a
    # root-split edge share the same original split and get the same label.
    for e in mid.postorder_edge_iter():
        if e.is_leaf():
            continue
        bp = e.bipartition
        side = sorted(t.label for t in bp.leafset_taxa(taxa))
        comp = sorted(set(all_tips) - set(side))
        if len(side) < 2 or len(comp) < 2:
            continue
        small = side if (len(side), side) <= (len(comp), comp) else comp
        key = hashlib.sha256(("\n".join(small) + "\n").encode()).hexdigest()
        lab = edge_supports.get(key)
        e.head_node.label = f"{lab[0]}/{lab[1]}" if lab else None
    mid.write(path=os.path.join(CORE, "main152_ML_midpoint_display.nwk"), schema="newick",
              suppress_rooting=True, real_value_format_specifier=".17g")

    # 3. invariance audit
    dm1 = ml.phylogenetic_distance_matrix()
    dm2 = mid.phylogenetic_distance_matrix()
    taxa_sorted = sorted(ml.taxon_namespace, key=lambda t: t.label)
    max_d = 0.0
    for i, t1 in enumerate(taxa_sorted):
        for t2 in taxa_sorted[i + 1:]:
            d = abs(dm1.distance(t1, t2) - dm2.distance(t1, t2))
            if d > max_d:
                max_d = d
    total_len_mid = sum(e.length or 0.0 for e in mid.postorder_edge_iter())
    # unrooted RF: reroot_at_midpoint leaves the copy flagged rooted, whose
    # bipartition encoding double-stores the root split; compare as unrooted
    ml_u = dendropy.Tree(ml)
    mid_u = dendropy.Tree(mid)
    ml_u.is_rooted = False
    mid_u.is_rooted = False
    ml_u.update_bipartitions()
    mid_u.update_bipartitions()
    rf = dendropy.calculate.treecompare.symmetric_difference(ml_u, mid_u)
    # support mapping: every nontrivial split key present in both trees
    keys_mid = set(edge_split_keys(mid, all_tips))
    keys_ml = set(edge_split_keys(ml, all_tips))
    # root-split edge in midpoint tree may duplicate one ml key; sets must be equal
    audit = dict(
        n_tips_unrooted=len(tips), n_tips_midpoint=len(tip_labels(mid)),
        tip_sets_identical=set(tip_labels(mid)) == set(roster),
        unrooted_RF=rf,
        max_abs_patristic_diff=max_d,
        total_length_unrooted=total_len_unrooted,
        total_length_midpoint=total_len_mid,
        total_length_abs_diff=abs(total_len_unrooted - total_len_mid),
        n_nontrivial_splits_unrooted=len(keys_ml),
        n_nontrivial_splits_midpoint=len(keys_mid),
        split_key_sets_equal=keys_ml == keys_mid,
        tolerance=TOL,
        status="PASS" if (rf == 0 and max_d <= TOL
                          and abs(total_len_unrooted - total_len_mid) <= TOL
                          and keys_ml == keys_mid
                          and set(tip_labels(mid)) == set(roster)) else "FAILED")
    with open(os.path.join(CORE, "midpoint_invariance_audit.json"), "w") as fh:
        json.dump(audit, fh, indent=1)
    print(json.dumps(audit, indent=1))
    if audit["status"] != "PASS":
        sys.exit(1)

    # 4. branch_support.tsv (unrooted, split_key-keyed)
    bs_rows = []
    for key, lab in sorted(edge_supports.items()):
        bs_rows.append(dict(split_key=key,
                            SH_aLRT=lab[0] if lab else "NA",
                            UFBoot=lab[1] if lab else "NA",
                            n_tips_smaller_side=min(edge_info[key]["n_side1"], edge_info[key]["n_side2"])))
    write_tsv(os.path.join(CORE, "branch_support.tsv"),
              ["split_key", "SH_aLRT", "UFBoot", "n_tips_smaller_side"], bs_rows)

    # 5. iTOL upload tree with unique internal IDs
    itol = dendropy.Tree.get(path=treefile, schema="newick",
                             taxon_namespace=dendropy.TaxonNamespace(),
                             preserve_underscores=True)
    itol.reroot_at_midpoint(update_bipartitions=True, suppress_unifurcations=False)
    all_tips_i = tip_labels(itol)

    # collect internal nodes with descendant tip lists
    node_desc = []
    for node in itol.postorder_node_iter():
        if node.is_leaf():
            continue
        desc = sorted(t.taxon.label for t in node.leaf_nodes())
        node_desc.append((desc, node))
    node_desc.sort(key=lambda x: x[0])
    id_of = {}
    for i, (desc, node) in enumerate(node_desc):
        nid = f"N{i + 1:04d}"
        id_of[id(node)] = nid

    # root split duplicate detection: the two root-child edges correspond to the
    # same original split (if nontrivial and supported)
    map_rows = []
    root = itol.seed_node
    # find which split_key each internal node corresponds to
    for desc, node in node_desc:
        nid = id_of[id(node)]
        if node is root:
            map_rows.append(dict(node_id=nid, split_key="NA", SH_aLRT="NA", UFBoot="NA",
                                 role="root", display_symbol="no",
                                 note="midpoint root; no support fabricated"))
            continue
        key = split_key(desc, all_tips_i)
        lab = edge_supports.get(key) if key else None
        map_rows.append(dict(node_id=nid, split_key=key or "NA_trivial",
                             SH_aLRT=lab[0] if lab else "NA",
                             UFBoot=lab[1] if lab else "NA",
                             role="internal", display_symbol="", note=""))

    # mark root-split duplicates: two child edges of root sharing one split key
    key_count = {}
    for r in map_rows:
        if r["split_key"] not in ("NA", "NA_trivial"):
            key_count[r["split_key"]] = key_count.get(r["split_key"], 0) + 1
    dup_keys = {k for k, c in key_count.items() if c > 1}
    for r in map_rows:
        if r["split_key"] in dup_keys:
            # decide display side: the rooted clade containing the global smallest tip
            # find the node's descendant tips via id order — we stored them in node_desc
            pass
    # second pass with descendant lists
    desc_of_nid = {}
    for desc, node in node_desc:
        desc_of_nid[id_of[id(node)]] = desc
    for r in map_rows:
        if r["split_key"] in dup_keys and r["role"] == "internal":
            desc = desc_of_nid[r["node_id"]]
            # the split's two sides: desc and complement; display on side with smallest tip
            comp = sorted(set(all_tips_i) - set(desc))
            side_with_min = desc if min(desc) < min(comp) else comp
            # this node displays only if its desc side IS the side with min tip
            # (the duplicate is the complementary rooting artifact)
            # compare with the other row carrying same key
            r["note"] = "root_split_candidate"
    # finalize duplicates: for each dup key pair, choose display once
    for key in dup_keys:
        rows_k = [r for r in map_rows if r["split_key"] == key and r["role"] == "internal"]
        sides = [desc_of_nid[r["node_id"]] for r in rows_k]
        # choose the row whose descendant side contains the global smallest tip
        winner = None
        for r, desc in zip(rows_k, sides):
            if min(desc) == min(all_tips_i):
                winner = r
        if winner is None:
            winner = rows_k[0]
        for r in rows_k:
            if r is winner:
                r["note"] = "root_split_display_once"
            else:
                r["note"] = "root_split_duplicate_not_counted"
                r["display_symbol"] = "no"

    for r in map_rows:
        if r["role"] == "internal" and r["display_symbol"] != "no":
            try:
                ok = (r["SH_aLRT"] != "NA" and r["UFBoot"] != "NA"
                      and float(r["SH_aLRT"]) >= 80 and float(r["UFBoot"]) >= 95)
            except (TypeError, ValueError):
                ok = False
            r["display_symbol"] = "yes" if ok else "no"

    # write iTOL upload tree: rename internal labels to N ids
    for desc, node in node_desc:
        node.label = id_of[id(node)]
    os.makedirs(FIG, exist_ok=True)
    os.makedirs(os.path.join(FIG, "reference"), exist_ok=True)
    os.makedirs(os.path.join(FIG, "annotations"), exist_ok=True)
    itol.write(path=os.path.join(FIG, "Figure2_main152_midpoint_iTOL.nwk"), schema="newick",
               suppress_rooting=True, real_value_format_specifier=".17g")

    write_tsv(os.path.join(FIG, "reference/itol_node_id_map.tsv"),
              ["node_id", "split_key", "SH_aLRT", "UFBoot", "role", "display_symbol", "note"],
              map_rows)

    # copy science record trees + branch support into reference/
    import shutil
    shutil.copy(os.path.join(CORE, "main152_ML_unrooted.nwk"),
                os.path.join(FIG, "reference/main152_ML_unrooted.nwk"))
    shutil.copy(os.path.join(CORE, "main152_ML_midpoint_display.nwk"),
                os.path.join(FIG, "reference/main152_ML_midpoint_display.nwk"))
    shutil.copy(os.path.join(CORE, "branch_support.tsv"),
                os.path.join(FIG, "reference/branch_support.tsv"))

    # 6. 05_branch_support.txt (iTOL DATASET_SYMBOL) — only edges passing 80/95
    sym_rows = [r for r in map_rows if r["display_symbol"] == "yes"]
    with open(os.path.join(FIG, "annotations/05_branch_support.txt"), "w") as fh:
        fh.write("DATASET_SYMBOL\nSEPARATOR TAB\nDATASET_LABEL\tSH-aLRT>=80 & UFBoot>=95\n")
        fh.write("COLOR\t#2166ac\nDATA\n")
        for r in sym_rows:
            fh.write(f"{r['node_id']}\t2\t3\t#2166ac\t1\t1\t"
                     f"SH-aLRT {r['SH_aLRT']}/UFBoot {r['UFBoot']}\n")
    n_sym = len(sym_rows)
    n_pass = sum(1 for r in map_rows if r["role"] == "internal")
    print(f"iTOL: {len(map_rows)} internal nodes, {n_sym} symbol edges (80/95), dup root splits: {len(dup_keys)}")
    with open(os.path.join(FIG, "validation_summary.json"), "w") as fh:
        json.dump(dict(n_internal_nodes=len(map_rows), n_symbol_edges=n_sym,
                       n_root_split_duplicates=len(dup_keys),
                       unique_internal_ids=len({r['node_id'] for r in map_rows}) == len(map_rows),
                       note="server delivers text inputs only; no image rendering"), fh, indent=1)


if __name__ == "__main__":
    main()
