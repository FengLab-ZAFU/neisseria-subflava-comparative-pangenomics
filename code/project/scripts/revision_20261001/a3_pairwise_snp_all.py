#!/usr/bin/env python3
"""A3: all-pairs nucmer/dnadiff SNP distances within each candidate group.

Complete-linkage clustering needs the full pairwise distance matrix, so every
C(n,2) pair in each candidate group is aligned (star-only would leave unknown
member-member distances at 0, wrongly merging clones).

Output: 02_independence/pairs_snp.tsv (full pairwise; group_id, pair ends as
reference/member, common length/ratio, ACGT-only SNP count, resolvability rule:
common >= 1 Mb AND >= 80% of left-genome length, else CLONE_UNRESOLVED).
"""
import csv
import re
import subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent.parent
OUT = ROOT / "results_revision_20261001" / "02_independence"
ALN = OUT / "alignments_groups"
ALN.mkdir(parents=True, exist_ok=True)
CONDA = "source ${CONDA_SH:-$HOME/miniconda3/etc/profile.d/conda.sh} && conda activate mummer_syri"
PAR = 32
ACGT = set("ACGT")

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


def run(cmd):
    return subprocess.run(["bash", "-c", f"{CONDA} && {cmd}"],
                          capture_output=True, text=True, timeout=3600)


def ref_len(fasta):
    n = 0
    for line in open(fasta):
        if not line.startswith(">"):
            n += len(line.strip())
    return n


def aligned_bases(rp):
    try:
        txt = Path(rp).read_text()
    except Exception:
        return None
    m = re.search(r"AlignedBases\s+(\d+)\S*\s+(\d+)", txt)
    return int(m.group(1)) if m else None


def count_snp(sp):
    n = 0
    if not Path(sp).exists():
        return 0
    for line in open(sp):
        f = line.split("\t")
        if len(f) < 3:
            continue
        a, b = f[0].strip().upper(), f[1].strip().upper()
        if a in ACGT and b in ACGT:
            n += 1
    return n


def align_pair(args):
    gid, gtype, a, b = args
    fa, fb = FASTA.get(a), FASTA.get(b)
    if not fa or not fb:
        return {"group_id": gid, "group_type": gtype, "reference": a, "member": b,
                "ref_length": "", "common_aligned_bases": "", "common_ratio": "",
                "pairwise_snp_ACGT": "", "resolvable_by_snp": "FASTA_MISSING"}
    d = ALN / f"{gid}__{a}__{b}"
    d.mkdir(exist_ok=True)
    prefix = d / "aln"
    if not Path(str(prefix) + ".report").exists():
        run(f"nucmer --maxmatch -p {prefix} {fa} {fb} 2>/dev/null")
        run(f"dnadiff -d {prefix}.delta -p {prefix} 2>/dev/null")
    rl = ref_len(fa)
    common = aligned_bases(str(prefix) + ".report") or 0
    snp = count_snp(str(prefix) + ".snps")
    ratio = common / rl if rl else 0
    resolvable = common >= 1_000_000 and ratio >= 0.80
    return {"group_id": gid, "group_type": gtype, "reference": a, "member": b,
            "ref_length": rl, "common_aligned_bases": common, "common_ratio": round(ratio, 4),
            "pairwise_snp_ACGT": snp,
            "resolvable_by_snp": "yes" if resolvable else "CLONE_UNRESOLVED"}


def main():
    load_fastas()
    groups = list(csv.DictReader(open(OUT / "candidate_groups.tsv"), delimiter="\t"))
    jobs = []
    for g in groups:
        mems = sorted(set(g["members"].split(";")))
        for i in range(len(mems)):
            for j in range(i + 1, len(mems)):
                jobs.append((g["group_id"], g["group_type"], mems[i], mems[j]))
    print(f"groups: {len(groups)} | pairs: {len(jobs)}")
    with ThreadPoolExecutor(max_workers=PAR) as ex:
        rows = list(ex.map(align_pair, jobs))
    fields = ["group_id", "group_type", "reference", "member", "ref_length",
              "common_aligned_bases", "common_ratio", "pairwise_snp_ACGT", "resolvable_by_snp"]
    with open(OUT / "pairs_snp.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    n_un = sum(1 for r in rows if r["resolvable_by_snp"] != "yes")
    print(f"wrote {len(rows)} pairs -> {OUT/'pairs_snp.tsv'} | unresolved/missing: {n_un}")


if __name__ == "__main__":
    main()
