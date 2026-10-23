#!/usr/bin/env python3
"""M5 step 1: protein-sequence SHA-256 bridge from historical eggNOG annotations
(emapper-3.0.0-beta6, revision_20261001 PGAP2 proteins) to the new main152 PGAP2
dictionary. No old cluster IDs are used; bridging is by normalized protein
sequence hash only. md5 is computed in parallel to cross-check against the
historical hash_bridge_db.tsv.

Normalization rule (both old and new): uppercase, strip all whitespace,
remove exactly one trailing '*'; internal '*' retained.

Outputs (06_function/):
  protein_hash_bridge.tsv   unique new-dictionary protein seqs -> COG/state
  copy_cog_map.tsv          per CDS copy: cluster/accession/gene_id/gene_index/hash/COG/state
  cluster_cog_map.tsv       per main cluster: representative protein annotation + conflicts
  bridge_qc.json            coverage + cross-check summary
"""
from pathlib import Path
import os
import csv, hashlib, json, re, sys
from collections import defaultdict

ROOT = str(Path(os.environ["NEISSERIA_PROJECT_ROOT"]))
PGAP2 = f"{ROOT}/results_revision_20261015_main152_species/02_pangenome/pgap2_main152"
OLD_PROT = f"{ROOT}/manuscript/Figures/new/Figure4_GO_KEGG_20261019/reference_inputs/total.involved_prot.fa"
OLD_ANNOT = f"{ROOT}/manuscript/Figures/new/Figure4_GO_KEGG_20261019/reference_inputs/pgap2_revision_primary.emapper.annotations"
HIST_BRIDGE = f"{ROOT}/results_revision_20261010_phaseB_v7/07_defense_function/cog/mapping/hash_bridge_db.tsv"
OUT = f"{ROOT}/results_revision_20261015_main152_species/06_function"

GI_RE = re.compile(r"^\d+:\d+:\d+$")

def norm_seq(s):
    s = re.sub(r"\s+", "", s).upper()
    if s.endswith("*"):
        s = s[:-1]
    return s

def hashes(seq):
    n = norm_seq(seq).encode()
    return hashlib.sha256(n).hexdigest(), hashlib.md5(n).hexdigest()

def read_fasta_indexed(path):
    """gene_index -> raw sequence (headers are '<gi> [..](.)')."""
    d = {}
    name, chunks = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    d[name] = "".join(chunks)
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
    if name is not None:
        d[name] = "".join(chunks)
    return d

def main():
    import os
    os.makedirs(OUT, exist_ok=True)

    # ---- 0. COG ID -> category letters dictionary
    cogid2letters = {}
    with open(f"{ROOT}/db/cog-24.def.tab") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 2 and f[0].startswith("COG"):
                cogid2letters[f[0]] = f[1]

    # ---- 1. old emapper annotations: gene_index -> COG letters
    # emapper-3.0.0-beta6 writes either category letters (e.g. "S", "KL") or a
    # COG ID (e.g. "COG2878") in the COG_category column; IDs are translated via
    # db/cog-24.def.tab (same rule as the historical hash bridge).
    old_annot = {}
    cogid_missing = set()
    with open(OLD_ANNOT) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            raw = f[7]
            if re.fullmatch(r"COG\d+", raw):
                letters = cogid2letters.get(raw)
                if letters is None:
                    cogid_missing.add(raw)
                    letters = ""
            elif raw == "-":
                letters = ""
            else:
                letters = raw
            old_annot[f[0]] = (letters, f[8])
    sys.stderr.write(f"old emapper rows: {len(old_annot)}; COG IDs missing from dict: {len(cogid_missing)}\n")

    # ---- 2. old proteins -> hash -> annotation state
    old_prot = read_fasta_indexed(OLD_PROT)
    sys.stderr.write(f"old proteins: {len(old_prot)}\n")
    # hash -> {"md5":..., "letters_set": set, "states": set, "n_gi": int}
    hmap = {}
    for gi, seq in old_prot.items():
        s256, m5 = hashes(seq)
        e = hmap.setdefault(s256, {"md5": m5, "letters": set(), "states": set(), "n_gi": 0})
        e["n_gi"] += 1
        if gi in old_annot:
            cog, _pname = old_annot[gi]
            if cog and cog != "-":
                e["letters"].update(cog)
                e["states"].add("annotated")
            else:
                e["states"].add("matched_no_cog")
        else:
            e["states"].add("matched_no_emapper_row")

    def bridge_state_letters(e):
        letters = "".join(sorted(e["letters"]))
        if letters:
            state = "annotated"
        elif "matched_no_cog" in e["states"]:
            state = "matched_no_cog"
        else:
            state = "matched_no_emapper_row"
        # conflict = mixed annotated/non-annotated gene_index for identical sequence
        conflict = 1 if (e["letters"] and ("matched_no_cog" in e["states"] or "matched_no_emapper_row" in e["states"])) else 0
        return letters, state, conflict

    # ---- 3. cross-check md5 vs historical bridge
    hist = {}
    with open(HIST_BRIDGE) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            hist[r["seq_md5"]] = (r["cog_letters"], r["annotation_state"])
    md5_to_sha = {}
    for s256, e in hmap.items():
        md5_to_sha[e["md5"]] = s256
    xcheck = {"hist_rows": len(hist), "md5_matched": 0, "md5_missing_in_new_bridge": 0,
              "state_concordant": 0, "state_discordant": 0, "letters_concordant": 0,
              "letters_discordant": 0, "discordant_examples": []}
    for m5, (hletters, hstate) in hist.items():
        s256 = md5_to_sha.get(m5)
        if s256 is None:
            xcheck["md5_missing_in_new_bridge"] += 1
            continue
        xcheck["md5_matched"] += 1
        letters, state, _ = bridge_state_letters(hmap[s256])
        if state == hstate:
            xcheck["state_concordant"] += 1
        else:
            xcheck["state_discordant"] += 1
            if len(xcheck["discordant_examples"]) < 10:
                xcheck["discordant_examples"].append([m5, hstate, state, hletters, letters])
        if letters == hletters:
            xcheck["letters_concordant"] += 1
        else:
            xcheck["letters_discordant"] += 1
    sys.stderr.write(f"md5 cross-check: {xcheck['md5_matched']}/{xcheck['hist_rows']} matched, "
                     f"state discordant {xcheck['state_discordant']}, letters discordant {xcheck['letters_discordant']}\n")

    # ---- 3b. optional supplementary emapper run (same version/DB snapshot) for
    # sequences absent from the historical bridge; queries are sha256 keys.
    supp = {}
    supp_path = None
    if "--supp" in sys.argv:
        supp_path = sys.argv[sys.argv.index("--supp") + 1]
        supp_cogid_missing = set()
        with open(supp_path) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 9:
                    continue
                raw = f[7]
                if re.fullmatch(r"COG\d+", raw):
                    letters = cogid2letters.get(raw)
                    if letters is None:
                        supp_cogid_missing.add(raw)
                        letters = ""
                elif raw == "-":
                    letters = ""
                else:
                    letters = raw
                state = "annotated" if letters else "matched_no_cog"
                supp[f[0]] = ("".join(sorted(set(letters))), state)
        sys.stderr.write(f"supplement annotation rows: {len(supp)}; COG IDs missing from dict: {len(supp_cogid_missing)}\n")

    # ---- 4. new run: gene_index -> (accession, gene_id); sequences
    gi2acc = {}
    with open(f"{PGAP2}/total.involved_annot.tsv") as fh:
        rdr = csv.reader(fh, delimiter="\t")
        header = next(rdr)
        for f in rdr:
            gi2acc[f[0]] = (f[1], f[5])  # Strain, Gene_ID
    new_prot = read_fasta_indexed(f"{PGAP2}/total.involved_prot.fa")
    sys.stderr.write(f"new annot rows: {len(gi2acc)}, new proteins: {len(new_prot)}\n")

    # sanity: annot strain+gene_id vs occurrence_map
    occ = []  # (cluster, accession, gene_id)
    with open(f"{PGAP2}/pav_tables/occurrence_map.tsv") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            occ.append((r["cluster"], r["accession"], r["gene_id"]))
    occ_set = {(c, a, g) for c, a, g in occ}
    sys.stderr.write(f"occurrence rows: {len(occ)}\n")

    # cross-validate cluster membership: detail table last column (','/';' split) vs occurrence_map
    detail_members = defaultdict(set)
    cluster_product = {}
    cluster_repre = {}
    with open(f"{PGAP2}/pgap2.partition.gene_content.detail.tsv") as fh:
        header = fh.readline()
        for line in fh:
            f = line.rstrip("\n").split("\t")
            cl = f[0]
            cluster_product[cl] = f[2]
            cluster_repre[cl] = f[4]
            for tok in re.split(r"[;,]", f[13]):
                tok = tok.strip()
                if GI_RE.match(tok):
                    detail_members[cl].add(tok)
    occ_members = defaultdict(set)
    gi_lookup = {}
    for gi, (acc, gid) in gi2acc.items():
        gi_lookup[(acc, gid)] = gi
    missing_join = 0
    occ_rows = []
    for cl, acc, gid in occ:
        gi = gi_lookup.get((acc, gid))
        if gi is None:
            missing_join += 1
        else:
            occ_members[cl].add(gi)
        occ_rows.append((cl, acc, gid, gi))
    membership_mismatch = [cl for cl in set(detail_members) | set(occ_members)
                           if detail_members.get(cl, set()) != occ_members.get(cl, set())]
    sys.stderr.write(f"occ->annot missing join: {missing_join}; "
                     f"cluster membership mismatches (detail vs occurrence): {len(membership_mismatch)}\n")

    # ---- 5. hash new proteins, bridge
    gi_hash = {}
    new_hashes = {}
    for gi, seq in new_prot.items():
        s256, m5 = hashes(seq)
        gi_hash[gi] = (s256, len(norm_seq(seq)))
        new_hashes.setdefault(s256, {"md5": m5, "n_copies": 0})
        new_hashes[s256]["n_copies"] += 1

    # supplement queries that produced no emapper row: evaluated, no annotation.
    # The submitted query set is exactly the new-dictionary sequences absent from
    # the historical bridge (preserved as supplement_emapper/input_unmatched_proteins.fasta).
    n_supp_norow = 0
    if supp_path:
        for s256 in new_hashes:
            if s256 not in hmap and s256 not in supp:
                supp[s256] = ("", "matched_no_emapper_row")
                n_supp_norow += 1

    def lookup(s256):
        e = hmap.get(s256)
        if e is not None:
            letters, state, _ = bridge_state_letters(e)
            return letters, state, "historical_bridge"
        s = supp.get(s256)
        if s is not None:
            return s[0], s[1], "supplement_emapper"
        return "", "no_sequence_match", "none"

    # protein_hash_bridge.tsv (unique sequences in new dictionary)
    with open(f"{OUT}/protein_hash_bridge.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["protein_sha256", "protein_md5", "n_copies_main152",
                    "cog_letters", "annotation_state", "annotation_source",
                    "hist_md5_letters", "hist_md5_state", "hist_md5_consistent"])
        for s256 in sorted(new_hashes):
            e = new_hashes[s256]
            letters, state, src = lookup(s256)
            hm = hist.get(e["md5"])
            if hm is None:
                hl, hs, hc = "NA", "NA", "NA"
            else:
                hl, hs = hm[0] or "", hm[1]
                hc = "yes" if (hl == letters and hs == state) else "no"
                if hl == "":
                    hl = "-"
            w.writerow([s256, e["md5"], e["n_copies"], letters if letters else "-",
                        state, src, hl, hs, hc])

    # unmatched unique sequences -> FASTA for optional supplementary emapper
    n_unmatched = 0
    with open(f"{OUT}/unmatched_proteins.fasta", "w") as fh:
        seen = set()
        for gi, seq in new_prot.items():
            s256, _l = gi_hash[gi]
            if s256 in seen:
                continue
            if hmap.get(s256) is None and (supp.get(s256) is None):
                fh.write(f">{s256}\n{norm_seq(seq)}\n")
                seen.add(s256)
                n_unmatched += 1
    sys.stderr.write(f"unique unmatched sequences written: {n_unmatched}\n")

    # ---- 6. copy-level map
    copy_counts = defaultdict(int)
    with open(f"{OUT}/copy_cog_map.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["cluster", "accession", "gene_id", "gene_index",
                    "protein_sha256", "protein_len_aa", "cog_letters", "annotation_state"])
        for cl, acc, gid, gi in occ_rows:
            if gi is None or gi not in gi_hash:
                copy_counts["unresolved_copy"] += 1
                continue
            s256, plen = gi_hash[gi]
            letters, state, _src = lookup(s256)
            copy_counts[state] += 1
            w.writerow([cl, acc, gid, gi, s256, plen, letters if letters else "-", state])

    # ---- 7. cluster-level map (representative = longest valid protein;
    #         ties by accession then gene_id lexicographic)
    cluster_copies = defaultdict(list)
    for cl, acc, gid, gi in occ_rows:
        if gi is not None and gi in gi_hash:
            cluster_copies[cl].append((acc, gid, gi))
    state_rank = {"annotated": 0, "matched_no_cog": 1, "matched_no_emapper_row": 2,
                  "no_sequence_match": 3}
    with open(f"{OUT}/cluster_cog_map.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["cluster", "repre_gene_index", "repre_accession", "repre_gene_id",
                    "repre_protein_sha256", "repre_protein_len_aa", "cog_letters",
                    "n_copies", "n_seq_matched", "n_annotated_copies",
                    "n_member_conflicts", "annotation_state", "product"])
        for cl in sorted(cluster_copies):
            copies = cluster_copies[cl]
            # representative: longest protein; tie -> accession, gene_id lexicographic
            rep = min(copies, key=lambda t: (-gi_hash[t[2]][1], t[0], t[1]))
            acc, gid, gi = rep
            s256, plen = gi_hash[gi]
            rep_letters, rep_state, _src = lookup(s256)
            n_matched = 0
            n_annot = 0
            n_conf = 0
            for a2, g2, gi2 in copies:
                s2, _l = gi_hash[gi2]
                l2, st2, _s2 = lookup(s2)
                if st2 != "no_sequence_match":
                    n_matched += 1
                if st2 == "annotated":
                    n_annot += 1
                    if l2 != rep_letters:
                        n_conf += 1
            w.writerow([cl, gi, acc, gid, s256, plen, rep_letters if rep_letters else "-",
                        len(copies), n_matched, n_annot, n_conf, rep_state,
                        cluster_product.get(cl, "NA")])

    qc = {
        "old_emapper_rows": len(old_annot),
        "old_proteins": len(old_prot),
        "old_unique_seq_hashes": len(hmap),
        "new_proteins": len(new_prot),
        "new_unique_seq_hashes": len(new_hashes),
        "occurrence_rows": len(occ_rows),
        "occ_to_annot_missing_joins": missing_join,
        "cluster_membership_mismatches_detail_vs_occurrence": len(membership_mismatch),
        "membership_mismatch_examples": sorted(membership_mismatch)[:20],
        "copy_state_counts": dict(copy_counts),
        "hist_md5_crosscheck": xcheck,
        "supplement_emapper_file": supp_path,
        "supplement_rows": len(supp),
        "supplement_queries_without_emapper_row": n_supp_norow,
        "unique_unmatched_after_supplement": n_unmatched,
        "normalization": "uppercase, whitespace removed, exactly one trailing '*' removed; internal '*' retained",
    }
    with open(f"{OUT}/bridge_qc.json", "w") as fh:
        json.dump(qc, fh, indent=1, ensure_ascii=False)
    sys.stderr.write("done\n")

if __name__ == "__main__":
    main()
