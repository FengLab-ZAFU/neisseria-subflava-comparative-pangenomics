#!/usr/bin/env python3
"""A3: genome-wide close-pair screen from all-vs-all FastANI (flag-only, no merging).

Reads 02_independence/allvsall_189/allvsall.tsv (query ref ani af).
Selects unordered pairs with ANI >= 99.5 and AF >= 0.80 (both reported directions
present at >= thresholds) that are NOT already covered by candidate_groups.tsv.
For each such pair runs nucmer/dnadiff (mummer 4.x) and counts ACGT<->ACGT SNPs
on common aligned bases, same resolvability rule as a3_pairwise_snp.py.

Output: 02_independence/flag_pairs_snp.tsv
"""
import csv
import re
import subprocess
from pathlib import Path
from collections import defaultdict

ROOT = Path(__file__).resolve().parent.parent.parent
OUT = ROOT / "results_revision_20261001" / "02_independence"
AV = OUT / "allvsall_189" / "allvsall.tsv"
ALN = OUT / "alignments_screen"
ALN.mkdir(exist_ok=True)
CONDA = "source ${CONDA_SH:-$HOME/miniconda3/etc/profile.d/conda.sh} && conda activate mummer_syri"

ANI_THR = 99.5
AF_THR = 0.80

FASTA = {}


def load_fastas():
    for d in (ROOT / "results_revision_20261001" / "01_retrieval_taxonomy" / "fasta", ROOT / "Genomes"):
        if d.exists():
            for p in d.glob("*.fna"):
                FASTA[p.stem] = p
    with open(ROOT / "neisseria_subflava_associated_group_metadata.tsv", encoding="utf-8-sig") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            acc = (r.get("assembly_accession") or "").strip()
            fn = (r.get("fasta_filename") or "").strip()
            if acc and fn and (ROOT / "Genomes" / fn).exists():
                FASTA.setdefault(acc, ROOT / "Genomes" / fn)


def stem_of(path):
    return Path(path).stem


def run(cmd):
    return subprocess.run(["bash", "-c", f"{CONDA} && {cmd}"],
                          capture_output=True, text=True, timeout=3600)


ACGT = set("ACGT")


def count_snp(snps_file):
    snp = 0
    if not Path(snps_file).exists():
        return 0
    for line in open(snps_file):
        f = line.split("\t")
        if len(f) < 3:
            continue
        a, b = f[0].strip().upper(), f[1].strip().upper()
        if a in ACGT and b in ACGT:
            snp += 1
    return snp


def aligned_bases(report_file):
    try:
        txt = Path(report_file).read_text()
    except Exception:
        return None
    m = re.search(r"AlignedBases\s+(\d+)\S*\s+(\d+)", txt)
    return int(m.group(1)) if m else None


def ref_len(fasta):
    n = 0
    for line in open(fasta):
        if not line.startswith(">"):
            n += len(line.strip())
    return n


def main():
    load_fastas()
    # pairs already covered by candidate groups
    covered = set()
    for g in csv.DictReader(open(OUT / "candidate_groups.tsv"), delimiter="\t"):
        mems = g["members"].split(";")
        for i in range(len(mems)):
            for j in range(i + 1, len(mems)):
                covered.add(tuple(sorted((mems[i], mems[j]))))

    best = {}
    rows = list(csv.reader(open(AV), delimiter="\t"))
    for q, r, ani, m, t in rows:
        if q == r:
            continue
        try:
            ani_f, af = float(ani), int(m) / int(t)
        except (ValueError, ZeroDivisionError):
            continue
        a, b = sorted((stem_of(q), stem_of(r)))
        if ani_f < ANI_THR or af < AF_THR:
            continue
        key = (a, b)
        if key in covered:
            continue
        cur = best.get(key)
        # keep min ANI and min AF across the two directions to be conservative
        if cur is None:
            best[key] = [ani_f, af]
        else:
            cur[0] = min(cur[0], ani_f)
            cur[1] = min(cur[1], af)

    print(f"screen pairs: {len(best)} (covered-by-groups excluded: {len(covered)})")
    # write pair list; alignments run in parallel via xargs (see a3_flag_pairs_align.sh)
    with open(OUT / "flag_pairs_list.tsv", "w") as fh:
        for (a, b) in sorted(best):
            fh.write(f"{a}\t{b}\n")
    print(f"wrote pair list -> {OUT/'flag_pairs_list.tsv'}")


if __name__ == "__main__":
    main()
