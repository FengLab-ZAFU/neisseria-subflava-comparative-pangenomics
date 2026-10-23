#!/usr/bin/env python3
"""P1 step 4: contig boundary / coordinate binding check for all 212 annotations
(reused cache + newly annotated).

Per accession (annotation dir from annotation_sources.tsv, frozen FASTA from
fasta_binding_212.tsv):
  1. contig set: seqids used in GFF3 feature lines must equal the contig IDs of
     the annotation .fna; the annotation .fna contig sequences must map
     bijectively (by exact sequence) onto the frozen FASTA contigs
     (bakta --compliant renames contigs to contig_N; the mapping table written
     to contig_id_maps/<acc>.tsv makes the renaming explicit and reversible).
  2. every CDS (seqid, start, end) lies within its contig bounds.
  3. every FFN/FAA record ID resolves to a GFF3 CDS (grouped by ID attribute).
  4. FFN nucleotide sequence == sequence extracted from the contig by
     (start, end, strand) for 100% of records.
  5. FAA == translation(table 11, start codons -> M) of the FFN sequence,
     allowing terminal-stop differences.

Outputs:
  01_annotation/annotation_contig_binding.tsv   one row per accession
  01_annotation/contig_id_maps/<acc>.tsv        contig rename map
"""
import os
import csv, json, sys
from pathlib import Path

ROOT = Path(os.environ["NEISSERIA_PROJECT_ROOT"])
OUT = ROOT / "results_revision_20261010_phaseB_v7/01_annotation"
MAPS = OUT / "contig_id_maps"
MAPS.mkdir(parents=True, exist_ok=True)

# NCBI translation table 11 (bacterial): standard code, starts TTG/CTG/GTG/ATT/ATC/ATA -> M
BASES = "TCAG"
AA = ("FFLLSSSSYY**CC*W" "LLLLPPPPHHQQRRRR" "IIIMTTTTNNKKSSRR" "VVVVAAAADDEEGGGG")
CODON = {a + b + c: AA[i * 16 + j * 4 + k]
         for i, a in enumerate(BASES) for j, b in enumerate(BASES) for k, c in enumerate(BASES)}
STARTS = {"TTG", "CTG", "GTG", "ATT", "ATC", "ATA", "ATG"}


def translate(seq, start_rule=True):
    seq = seq.upper()
    aa = []
    for i in range(0, len(seq) - 2, 3):
        c = seq[i:i + 3]
        if i == 0 and start_rule and c in STARTS:
            aa.append("M")
        else:
            aa.append(CODON.get(c, "X"))
    return "".join(aa)


_IUPAC = str.maketrans("ACGTRYSWKMBDHVNacgtryswkmbdhvn",
                       "TGCAYRSWMKVHDBNtgcayrswmkvhdbn")


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
    """contig id -> True if bakta header marks [topology=circular]."""
    topo = {}
    with open(p) as fh:
        for line in fh:
            if line.startswith(">"):
                topo[line[1:].split()[0]] = "[topology=circular]" in line
    return topo


def read_tsv(p):
    with open(p, newline="", encoding="utf-8-sig") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


FEATURE_TYPES = {"CDS", "tRNA", "rRNA", "tmRNA", "ncRNA", "ncRNA_region",
                 "regulatory", "regulatory_region", "rep_origin", "repeat_region", "sORF"}


def parse_gff_cds(p):
    """Return (seqid_set, cds dict, feats dict).
    cds:  id -> list of (seqid,start,end,strand) for CDS features.
    feats: id -> list of (seqid,start,end,strand) for every FFN-exported
           feature type (CDS + RNA + regulatory/rep_origin/repeat_region).
           bakta writes non-locus-tag features (regulatory, rep_origin,
           repeat_region) with GFF IDs (featNNN) that differ from their FFN
           record IDs, so FFN records that miss the ID lookup are resolved by
           exact sequence match against all feature extractions."""
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


def check_one(acc, ann_dir, fasta_path):
    d = Path(ann_dir)
    gff, faa, ffn, ann_fna = (d / f"{acc}.{e}" for e in ("gff3", "faa", "ffn", "fna"))
    notes, hard_fail = [], []

    frozen = read_fasta(fasta_path)
    ann = read_fasta(ann_fna) if ann_fna.is_file() else {}
    circular = read_fasta_topology(ann_fna) if ann_fna.is_file() else {}
    seqids, cds, feats = parse_gff_cds(gff)
    ffn_recs = read_fasta(ffn)
    faa_recs = read_fasta(faa)

    # 1. contig sets
    gff_seqids_in_annfna = seqids <= set(ann)
    if not gff_seqids_in_annfna:
        hard_fail.append("gff_seqids_not_in_annotation_fna")
    # map annotation contigs -> frozen contigs by exact sequence
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
    with open(MAPS / f"{acc}.contig_map.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["annotation_contig_id", "frozen_contig_id", "length_bp"])
        w.writerows(sorted(mapping, key=lambda m: m[0]))
    contigs_match = (gff_seqids_in_annfna and bijective)
    if not bijective:
        hard_fail.append(f"contig_map_not_bijective(unmatched={unmatched},ambiguous={ambiguous},"
                         f"ann={len(ann)},frozen={len(frozen)})")
    renamed = any(a != b for a, b, _ in mapping)

    # 2. coordinate bounds (all FFN-exported feature types; rep_origin and
    # repeat_region features may carry strand "?" or "."). On circular contigs
    # bakta unwraps origin-spanning features past the contig end (end > L);
    # those are legal iff the contig is circular, end <= 2L and the wrapped
    # part does not overlap the start (single wrap).
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
    coord_ok = bounds_bad == 0
    if not coord_ok:
        hard_fail.append(f"feature_out_of_bounds={bounds_bad}")
    if n_wrapped:
        notes.append(f"origin_spanning_wrapped_features={n_wrapped}")

    def extract_one(seqid, s, e):
        seq = ann[seqid]
        L = len(seq)
        if e <= L:
            return seq[s - 1:e]
        return seq[s - 1:] + seq[:e - L]  # circular wrap

    def extract(parts):
        parts = sorted(parts, key=lambda p: p[1], reverse=(parts[0][3] == "-"))
        sub = "".join(extract_one(p0, p1, p2) for p0, p1, p2, _ in parts)
        return revcomp(sub) if parts[0][3] == "-" else sub

    # 3./4. FFN records resolve to a GFF feature (by ID, else by exact
    # extracted-sequence match) and reproduce the FASTA extraction (100%)
    feat_extractions = {}
    for cid, parts in feats.items():
        feat_extractions.setdefault(extract(parts), []).append(cid)
    # bakta exports oriT (and possibly other non-GFF3 feature types) to FFN but
    # not to GFF3; resolve those against the annotation JSON coordinates.
    json_feats = {}
    jpath = d / f"{acc}.json"
    if jpath.is_file():
        try:
            jdata = json.loads(jpath.read_text())
            for jf in jdata.get("features", []):
                seqid, s, e, st = jf.get("sequence"), jf.get("start"), jf.get("stop"), jf.get("strand")
                if seqid in ann and isinstance(s, int) and isinstance(e, int):
                    json_feats.setdefault(jf.get("id"), (seqid, s, e, st))
        except Exception:
            pass

    ffn_id_missing, ffn_seq_bad, ffn_alias, ffn_json = 0, 0, 0, 0
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

    # 5. FAA translation check (terminal stop tolerated; for partial CDS bakta
    # does not apply the alternative-start -> M rule, so the literal table-11
    # translation of the first codon is also accepted)
    faa_id_missing = sum(1 for rid in faa_recs if rid not in cds)
    faa_bad, faa_partial_start = 0, 0
    for rid, prot in faa_recs.items():
        nt = ffn_recs.get(rid)
        if nt is None:
            faa_bad += 1
            continue
        p = prot.rstrip("*")
        if translate(nt).rstrip("*") == p:
            continue
        if translate(nt, start_rule=False).rstrip("*") == p:
            faa_partial_start += 1
            continue
        faa_bad += 1
    if faa_id_missing:
        hard_fail.append(f"faa_ids_not_in_gff_cds={faa_id_missing}")
    if faa_bad:
        hard_fail.append(f"faa_translation_mismatch={faa_bad}")
    if faa_partial_start:
        notes.append(f"faa_first_residue_literal_not_M(partial_cds)={faa_partial_start}")

    if len(ffn_recs) != len(feats):
        notes.append(f"n_ffn({len(ffn_recs)})!=n_ffn_features_gff({len(feats)})")
    if len(faa_recs) != len(cds):
        notes.append(f"n_faa({len(faa_recs)})!=n_cds_gff({len(cds)})")

    status = "PASS" if not hard_fail else "FAIL"
    return {
        "accession": acc,
        "annotation_dir": str(d),
        "n_contigs_frozen": len(frozen), "n_contigs_annotation": len(ann),
        "n_contigs_gff": len(seqids), "n_contigs_match": "yes" if contigs_match else "no",
        "contig_renamed": "yes" if renamed else "no", "contig_map_bijective": "yes" if bijective else "no",
        "n_cds_gff": len(cds), "n_ffn_features_gff": len(feats),
        "n_faa": len(faa_recs), "n_ffn": len(ffn_recs),
        "seq_identity_check": "pass" if (ffn_id_missing == 0 and ffn_seq_bad == 0) else "fail",
        "coord_bounds_check": "pass" if coord_ok else "fail",
        "faa_translation_check": "pass" if (faa_id_missing == 0 and faa_bad == 0) else "fail",
        "status": status,
        "detail": ";".join(hard_fail),
        "notes": ";".join(notes),
    }


def main():
    sources = {r["accession"]: r for r in read_tsv(OUT / "annotation_sources.tsv")}
    binding = {r["accession"]: r for r in read_tsv(OUT / "fasta_binding_212.tsv")}
    rows = []
    for acc, s in sorted(sources.items()):
        b = binding[acc]
        if b["match"] != "PASS":
            rows.append({"accession": acc, "annotation_dir": s["annotation_dir"],
                         "status": "FAIL", "detail": "frozen_fasta_binding_not_pass",
                         "n_contigs_frozen": "", "n_contigs_annotation": "", "n_contigs_gff": "",
                         "n_contigs_match": "", "contig_renamed": "", "contig_map_bijective": "",
                         "n_cds_gff": "", "n_faa": "", "n_ffn": "", "seq_identity_check": "",
                         "coord_bounds_check": "", "faa_translation_check": "", "notes": ""})
            continue
        try:
            rows.append(check_one(acc, s["annotation_dir"], b["path"]))
        except Exception as e:
            rows.append({"accession": acc, "annotation_dir": s["annotation_dir"],
                         "status": "FAIL", "detail": f"exception:{type(e).__name__}:{e}",
                         "n_contigs_frozen": "", "n_contigs_annotation": "", "n_contigs_gff": "",
                         "n_contigs_match": "", "contig_renamed": "", "contig_map_bijective": "",
                         "n_cds_gff": "", "n_faa": "", "n_ffn": "", "seq_identity_check": "",
                         "coord_bounds_check": "", "faa_translation_check": "", "notes": ""})
    fields = ["accession", "annotation_dir", "n_contigs_frozen", "n_contigs_annotation",
              "n_contigs_gff", "n_contigs_match", "contig_renamed", "contig_map_bijective",
              "n_cds_gff", "n_ffn_features_gff", "n_faa", "n_ffn", "seq_identity_check",
              "coord_bounds_check", "faa_translation_check", "status", "detail", "notes"]
    with open(OUT / "annotation_contig_binding.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader(); w.writerows(rows)
    npass = sum(1 for r in rows if r["status"] == "PASS")
    print(f"contig binding: {npass}/{len(rows)} PASS")
    for r in rows:
        if r["status"] != "PASS":
            print("FAIL", r["accession"], r["detail"])
        elif r["notes"]:
            print("NOTE", r["accession"], r["notes"])


if __name__ == "__main__":
    main()
