#!/usr/bin/env python3
"""M1: per-genome annotation binding verification + GFF staging (main152).

Binding logic mirrors the historically verified Phase-B checker
(scripts/revision_20261010_phaseB_v7/p1_contig_binding_check.py):
  1. GFF seqids subset of annotation .fna; annotation contigs map bijectively
     onto frozen FASTA contigs by exact sequence (bakta renames to contig_N;
     per-accession map written to contig_id_maps/).
  2. All feature coords within bounds; circular origin-wrapping legal only on
     bakta [topology=circular] contigs with single wrap.
  3./4. Every FFN record resolves to a GFF feature (by ID, else exact sequence,
     else bakta JSON) and reproduces the extraction.
  5. FAA == table-11 translation of FFN (terminal stop / literal first residue
     for partial CDS tolerated).
Plus main152 additions: annotation file sha256 ledger, reuse ledger, and
PGAP2/Panaroo GFF staging (##FASTA appended from the *annotation* sibling fna
after bijective sequence verification, keeping contig_N seqids consistent).
"""
import csv
import json
import os
import sys
from multiprocessing import Pool
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import (RESULTS, PACKET, read_roster, read_tsv, write_tsv,
                    resolve_path, sha256_file, sha256_bytes)

OUT = os.path.join(RESULTS, "01_annotation")
STAGE = os.path.join(OUT, "gff_main152")
MAPS = os.path.join(OUT, "contig_id_maps")

BASES = "TCAG"
AA = ("FFLLSSSSYY**CC*W" "LLLLPPPPHHQQRRRR" "IIIMTTTTNNKKSSRR" "VVVVAAAADDEEGGGG")
CODON = {a + b + c: AA[i * 16 + j * 4 + k]
         for i, a in enumerate(BASES) for j, b in enumerate(BASES) for k, c in enumerate(BASES)}
STARTS = {"TTG", "CTG", "GTG", "ATT", "ATC", "ATA", "ATG"}
_IUPAC = str.maketrans("ACGTRYSWKMBDHVNacgtryswkmbdhvn", "TGCAYRSWMKVHDBNtgcayrswmkvhdbn")
FEATURE_TYPES = {"CDS", "tRNA", "rRNA", "tmRNA", "ncRNA", "ncRNA_region",
                 "regulatory", "regulatory_region", "rep_origin", "repeat_region", "sORF"}


def translate(seq, start_rule=True):
    seq = seq.upper()
    aa = []
    for i in range(0, len(seq) - 2, 3):
        c = seq[i:i + 3]
        aa.append("M" if (i == 0 and start_rule and c in STARTS) else CODON.get(c, "X"))
    return "".join(aa)


def revcomp(s):
    return s.translate(_IUPAC)[::-1]


def read_fasta(p):
    recs, name, cur = {}, None, []
    with open(p) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    recs[name] = "".join(cur).upper()
                name = line[1:].split()[0]
                cur = []
            else:
                cur.append(line.strip())
    if name is not None:
        recs[name] = "".join(cur).upper()
    return recs


def read_fasta_topology(p):
    topo = {}
    with open(p) as fh:
        for line in fh:
            if line.startswith(">"):
                topo[line[1:].split()[0]] = "[topology=circular]" in line
    return topo


def parse_gff(p):
    seqids, cds, feats = set(), {}, {}
    with open(p) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            seqids.add(f[0])
            if f[2] not in FEATURE_TYPES:
                continue
            attrs = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
            cid = attrs.get("ID") or attrs.get("locus_tag")
            feats.setdefault(cid, []).append((f[0], int(f[3]), int(f[4]), f[6]))
            if f[2] in ("CDS", "sORF"):
                cds.setdefault(cid, []).append((f[0], int(f[3]), int(f[4]), f[6]))
    return seqids, cds, feats


BINDING = None


def init():
    global BINDING
    BINDING = {r["accession"]: r for r in read_tsv(os.path.join(RESULTS, "00_preflight/input_check/input_binding.tsv"))}


def check_one(acc_row):
    acc = acc_row["accession"]
    res = {"accession": acc, "annotation_source": acc_row["annotation_source"]}
    hard_fail, notes = [], []

    paths = {}
    for col in ("gff_path", "gbff_path", "faa_path", "ffn_path", "log_path"):
        act, note = resolve_path(acc_row[col])
        paths[col] = act
        res[col + "_resolved"] = act or ""
        res[col + "_sha256"] = sha256_file(act) if act else ""
        if not act:
            hard_fail.append(f"missing_{col}")
    ann_dir = os.path.dirname(paths["gff_path"]) if paths["gff_path"] else None
    res["annotation_dir"] = ann_dir or ""
    if not ann_dir:
        res["status"] = "FAILED"
        res["detail"] = ";".join(hard_fail)
        res["notes"] = ""
        return res

    gff, faa, ffn = paths["gff_path"], paths["faa_path"], paths["ffn_path"]
    ann_fna = os.path.join(ann_dir, acc + ".fna")
    fasta_path = BINDING[acc]["path"]

    frozen = read_fasta(fasta_path)
    ann = read_fasta(ann_fna) if os.path.isfile(ann_fna) else {}
    circular = read_fasta_topology(ann_fna) if os.path.isfile(ann_fna) else {}
    seqids, cds, feats = parse_gff(gff)
    ffn_recs = read_fasta(ffn) if faa else {}
    faa_recs = read_fasta(faa) if faa else {}

    # 1. contig sets
    if not seqids <= set(ann):
        hard_fail.append("gff_seqids_not_in_annotation_fna")
    frozen_by_seq = {}
    for cid, s in frozen.items():
        frozen_by_seq.setdefault(s, []).append(cid)
    mapping, ambiguous, unmatched = [], 0, 0
    for cid, s in ann.items():
        hits = frozen_by_seq.get(s, [])
        if len(hits) == 1:
            mapping.append((cid, hits[0], len(s)))
        elif len(hits) > 1:
            ambiguous += 1
            mapping.append((cid, "|".join(hits), len(s)))
        else:
            unmatched += 1
    bijective = (unmatched == 0 and ambiguous == 0
                 and len({m[1] for m in mapping}) == len(frozen) == len(ann))
    os.makedirs(MAPS, exist_ok=True)
    with open(os.path.join(MAPS, f"{acc}.contig_map.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["annotation_contig_id", "frozen_contig_id", "length_bp"])
        w.writerows(sorted(mapping, key=lambda m: m[0]))
    if not bijective:
        hard_fail.append(f"contig_map_not_bijective(unmatched={unmatched},ambiguous={ambiguous},"
                         f"ann={len(ann)},frozen={len(frozen)})")
    renamed = any(a != b for a, b, _ in mapping)

    # 2. bounds
    bounds_bad, n_wrapped = 0, 0
    for cid, parts in feats.items():
        for seqid, start, end, strand in parts:
            L = len(ann.get(seqid, ""))
            if strand not in ("+", "-", "?", ".") or start < 1 or end < start:
                bounds_bad += 1
            elif end > L:
                if circular.get(seqid) and end <= 2 * L and (end - L) < start:
                    n_wrapped += 1
                else:
                    bounds_bad += 1
            elif start > L:
                bounds_bad += 1
    if bounds_bad:
        hard_fail.append(f"feature_out_of_bounds={bounds_bad}")
    if n_wrapped:
        notes.append(f"origin_spanning_wrapped_features={n_wrapped}")

    def extract_one(seqid, s, e):
        seq = ann[seqid]
        L = len(seq)
        if e <= L:
            return seq[s - 1:e]
        return seq[s - 1:] + seq[:e - L]

    def extract(parts):
        parts = sorted(parts, key=lambda p: p[1], reverse=(parts[0][3] == "-"))
        sub = "".join(extract_one(p0, p1, p2) for p0, p1, p2, _ in parts)
        return revcomp(sub) if parts[0][3] == "-" else sub

    # 3./4. FFN
    feat_extractions = {}
    for cid, parts in feats.items():
        feat_extractions.setdefault(extract(parts), []).append(cid)
    json_feats = {}
    jpath = os.path.join(ann_dir, f"{acc}.json")
    if os.path.isfile(jpath):
        try:
            jdata = json.loads(open(jpath).read())
            for jf in jdata.get("features", []):
                seqid, s, e, st = jf.get("sequence"), jf.get("start"), jf.get("stop"), jf.get("strand")
                if seqid in ann and isinstance(s, int) and isinstance(e, int):
                    json_feats.setdefault(jf.get("id"), (seqid, s, e, st))
        except Exception:
            pass
    ffn_id_missing = ffn_seq_bad = ffn_alias = ffn_json = 0
    for rid, seq in ffn_recs.items():
        parts = feats.get(rid)
        if parts:
            if extract(parts) != seq.upper():
                ffn_seq_bad += 1
        elif seq.upper() in feat_extractions:
            ffn_alias += 1
        elif rid in json_feats:
            seqid, s, e, st = json_feats[rid]
            sub = extract_one(seqid, s, e)
            if st == "-":
                sub = revcomp(sub)
            if sub == seq.upper():
                ffn_json += 1
            else:
                ffn_seq_bad += 1
        else:
            ffn_id_missing += 1
    if ffn_id_missing:
        hard_fail.append(f"ffn_records_unresolvable={ffn_id_missing}")
    if ffn_seq_bad:
        hard_fail.append(f"ffn_seq_mismatch={ffn_seq_bad}")
    if ffn_alias:
        notes.append(f"ffn_ids_resolved_by_sequence_not_gff_id={ffn_alias}")
    if ffn_json:
        notes.append(f"ffn_ids_resolved_via_bakta_json_not_in_gff={ffn_json}")

    # 5. FAA
    faa_id_missing = sum(1 for rid in faa_recs if rid not in cds)
    faa_bad = faa_partial = 0
    for rid, prot in faa_recs.items():
        nt = ffn_recs.get(rid)
        if nt is None:
            faa_bad += 1
            continue
        p = prot.rstrip("*")
        if translate(nt).rstrip("*") == p:
            continue
        if translate(nt, start_rule=False).rstrip("*") == p:
            faa_partial += 1
            continue
        faa_bad += 1
    if faa_id_missing:
        hard_fail.append(f"faa_ids_not_in_gff_cds={faa_id_missing}")
    if faa_bad:
        hard_fail.append(f"faa_translation_mismatch={faa_bad}")
    if faa_partial:
        notes.append(f"faa_first_residue_literal_not_M={faa_partial}")
    if len(ffn_recs) != len(feats):
        notes.append(f"n_ffn({len(ffn_recs)})!=n_ffn_features_gff({len(feats)})")
    if len(faa_recs) != len(cds):
        notes.append(f"n_faa({len(faa_recs)})!=n_cds_gff({len(cds)})")

    # staging: PGAP2/Panaroo parse the GFF ##FASTA section. If the historical
    # gff3 lacks it, append the *annotation* sibling fna (contig_N ids match the
    # GFF seqids; bijective sequence map to frozen FASTA verified above).
    staged_sha = ""
    fasta_appended = ""
    if not hard_fail:
        gff_bytes = open(gff, "rb").read()
        if b"##FASTA" not in gff_bytes:
            if not os.path.isfile(ann_fna):
                hard_fail.append("gff_lacks_FASTA_and_annotation_fna_missing")
            else:
                if not gff_bytes.endswith(b"\n"):
                    gff_bytes += b"\n"
                gff_bytes += b"##FASTA\n" + open(ann_fna, "rb").read()
                fasta_appended = "yes_annotation_fna"
        else:
            # embedded FASTA must equal annotation fna sequences
            idx = gff_bytes.find(b"##FASTA")
            emb = {}
            name, buf = None, []
            for ln in gff_bytes[idx:].decode("utf-8", "strict").splitlines()[1:]:
                if ln.startswith(">"):
                    if name is not None:
                        emb[name] = "".join(buf).upper()
                    name, buf = ln[1:].split()[0], []
                else:
                    buf.append(ln.strip())
            if name is not None:
                emb[name] = "".join(buf).upper()
            if set(emb) != set(ann) or any(emb[k] != ann[k] for k in emb):
                hard_fail.append("embedded_fasta_differs_from_annotation_fna")
            fasta_appended = "no_embedded_verified"
        if not hard_fail:
            dst = os.path.join(STAGE, acc + ".gff3")
            with open(dst, "wb") as fh:
                fh.write(gff_bytes)
            staged_sha = sha256_bytes(gff_bytes)

    res.update(dict(
        n_contigs_frozen=len(frozen), n_contigs_annotation=len(ann), n_contigs_gff=len(seqids),
        contig_renamed="yes" if renamed else "no", contig_map_bijective="yes" if bijective else "no",
        n_cds_gff=len(cds), n_ffn_features_gff=len(feats), n_faa=len(faa_recs), n_ffn=len(ffn_recs),
        fasta_appended=fasta_appended, staged_sha256=staged_sha,
        status="PASS" if not hard_fail else "FAILED",
        detail=";".join(hard_fail), notes=";".join(notes)))
    return res


def main():
    os.makedirs(STAGE, exist_ok=True)
    roster = read_roster()
    ann = {r["accession"]: r for r in read_tsv(os.path.join(PACKET, "00_input/reference_sources/annotation_manifest_main152.tsv"))}
    rows_in = [ann[a] for a in roster]
    with Pool(32, initializer=init) as pool:
        results = pool.map(check_one, rows_in)

    fields = ["accession", "annotation_source", "annotation_dir", "status", "detail", "notes",
              "n_contigs_frozen", "n_contigs_annotation", "n_contigs_gff",
              "contig_renamed", "contig_map_bijective", "n_cds_gff", "n_ffn_features_gff",
              "n_faa", "n_ffn", "fasta_appended", "staged_sha256",
              "gff_path_resolved", "gff_path_sha256", "gbff_path_resolved", "gbff_path_sha256",
              "faa_path_resolved", "faa_path_sha256", "ffn_path_resolved", "ffn_path_sha256",
              "log_path_resolved", "log_path_sha256"]
    write_tsv(os.path.join(OUT, "annotation_manifest_main152.tsv"), fields, results)
    write_tsv(os.path.join(OUT, "contig_binding_check.tsv"),
              ["accession", "n_contigs_frozen", "n_contigs_annotation", "n_contigs_gff",
               "contig_renamed", "contig_map_bijective", "status"],
              [dict(**{k: r.get(k, "") for k in ["accession", "n_contigs_frozen", "n_contigs_annotation",
                                                 "n_contigs_gff", "contig_renamed", "contig_map_bijective", "status"]})
               for r in results])
    write_tsv(os.path.join(OUT, "cds_binding_check.tsv"),
              ["accession", "n_cds_gff", "n_faa", "n_ffn", "status", "detail", "notes"],
              [{k: r.get(k, "") for k in ["accession", "n_cds_gff", "n_faa", "n_ffn", "status", "detail", "notes"]}
               for r in results])
    write_tsv(os.path.join(OUT, "annotation_reuse_ledger.tsv"),
              ["accession", "annotation_source", "decision", "evidence"],
              [dict(accession=r["accession"], annotation_source=r["annotation_source"],
                    decision="reuse" if r["status"] == "PASS" else "rerun_required",
                    evidence=(r["detail"] or "contig/CDS/FFN/FAA sequence-level binding PASS")
                             + (" | notes: " + r["notes"] if r["notes"] else "")) for r in results])
    write_tsv(os.path.join(OUT, "gff_staging_log.tsv"),
              ["accession", "staged_gff", "fasta_appended", "staged_sha256", "n_cds", "status"],
              [dict(accession=r["accession"], staged_gff=os.path.join(STAGE, r["accession"] + ".gff3"),
                    fasta_appended=r.get("fasta_appended", ""), staged_sha256=r.get("staged_sha256", ""),
                    n_cds=r.get("n_cds_gff", ""), status=r["status"]) for r in results])

    n_pass = sum(1 for r in results if r["status"] == "PASS")
    print(f"M1 binding: {n_pass}/152 PASS")
    for r in results:
        if r["status"] != "PASS":
            print("FAILED:", r["accession"], r["detail"])
    sys.exit(0 if n_pass == 152 else 1)


if __name__ == "__main__":
    main()
