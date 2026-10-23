#!/usr/bin/env python3
"""A3: pairwise whole-genome SNP distance via fixed-reference coordinate alignment.

For each candidate group: pick reference by quality (complete>scaffold>contig,
then N50 desc, contigs asc, accession lexicographic). Align each member to the
reference with nucmer (mummer 4.0.1), summarize with dnadiff, and count
ACGT<->ACGT substitutions from the .snps file.

Outputs per group into 02_independence/alignments/<group_id>/:
  ref.txt, <acc>.report, <acc>.snps.filtered, pairs_snp.tsv (merged later)

Operational rules (task doc A3):
  - pairwise SNP only on common A/C/G/T sites
  - report common length, reference length, ratio
  - merging allowed by SNP only if common >= 1 Mb AND >= 80% of reference length
    else CLONE_UNRESOLVED
"""
import csv
import re
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent.parent
OUT = ROOT / "results_revision_20261001" / "02_independence"
ALN = OUT / "alignments"
ALN.mkdir(parents=True, exist_ok=True)
CONDA = "source ${CONDA_SH:-$HOME/miniconda3/etc/profile.d/conda.sh} && conda activate mummer_syri"

FASTA = {}  # accession -> fasta path (filled from main fasta dir + Genomes)


def load_fastas():
    d1 = ROOT / "results_revision_20261001" / "01_retrieval_taxonomy" / "fasta"
    d2 = ROOT / "Genomes"
    for d in (d1, d2):
        if d.exists():
            for p in d.glob("*.fna"):
                FASTA[p.stem] = p
    # alias historical Genomes/ files (BioSample-named) to their assembly accessions
    import csv as _csv
    meta_tsv = ROOT / "neisseria_subflava_associated_group_metadata.tsv"
    with open(meta_tsv, encoding="utf-8-sig") as fh:
        for r in _csv.DictReader(fh, delimiter="\t"):
            acc = (r.get("assembly_accession") or "").strip()
            fn = (r.get("fasta_filename") or "").strip()
            if acc and fn:
                p = ROOT / "Genomes" / fn
                if p.exists():
                    FASTA.setdefault(acc, p)


def ref_len(fasta):
    n = 0
    for line in open(fasta):
        if not line.startswith(">"):
            n += len(line.strip())
    return n


def n_contigs(fasta):
    return sum(1 for line in open(fasta) if line.startswith(">"))


def pick_ref(members, meta):
    def key(acc):
        m = meta.get(acc, {})
        level = m.get("assembly_level", "")
        rank = {"Complete Genome": 0, "Chromosome": 1, "Scaffold": 2, "Contig": 3}.get(level, 4)
        n50 = -(int(m.get("contig_n50", 0) or 0))
        nc = int(m.get("n_contigs", 9e9) or 9e9)
        return (rank, n50, nc, acc)
    return sorted(members, key=key)[0]


def run(cmd):
    return subprocess.run(["bash", "-c", f"{CONDA} && {cmd}"],
                          capture_output=True, text=True, timeout=3600)


ACGT = set("ACGT")


def count_snp(snps_file):
    """Count ACGT<->ACGT substitutions from dnadiff .snps; return (snp, sites_ok)."""
    snp = 0
    if not Path(snps_file).exists():
        return 0, 0
    for line in open(snps_file):
        f = line.split("\t")
        if len(f) < 3:
            continue
        a, b = f[0].strip().upper(), f[1].strip().upper()
        if a in ACGT and b in ACGT:
            snp += 1
    return snp, 1


def aligned_bases(report_file):
    """Parse dnadiff report: total aligned bases (1-1 + 1-M? use 1-1 only)."""
    aligned = None
    total_snps = None
    try:
        txt = Path(report_file).read_text()
        m = re.search(r"AlignedBases\s+(\d+)\S*\s+(\d+)", txt)
        if m:
            aligned = int(m.group(1))  # ref-side aligned bases
        m2 = re.search(r"TotalSNPs\s+(\d+)\s+(\d+)", txt)
        if m2:
            total_snps = int(m2.group(1))
    except Exception:
        pass
    return aligned, total_snps


def main():
    load_fastas()
    meta = {}
    meta_tsv = ROOT / "results_revision_20261001" / "01_retrieval_taxonomy" / "screening_scope_accessions.tsv"
    if meta_tsv.exists():
        for r in csv.DictReader(open(meta_tsv), delimiter="\t"):
            meta[r["accession"]] = {"assembly_level": r["assembly_level"],
                                    "contig_n50": r["contig_n50"],
                                    "n_contigs": r["n_contigs"]}
    # historical 189 metadata for those not in scope file
    for r in csv.DictReader(open(ROOT / "neisseria_subflava_associated_group_metadata.tsv", encoding="utf-8-sig"),
                            delimiter="\t"):
        acc = r["assembly_accession"]
        meta.setdefault(acc, {"assembly_level": r.get("assembly_level", ""),
                              "contig_n50": r.get("contig_n50_bp", ""),
                              "n_contigs": r.get("number_of_contigs", "")})

    groups = list(csv.DictReader(open(OUT / "candidate_groups.tsv"), delimiter="\t"))
    pair_rows = []
    for g in groups:
        members = g["members"].split(";")
        present = [m for m in members if m in FASTA]
        if len(present) < 2:
            continue
        ref = pick_ref(present, meta)
        gdir = ALN / g["group_id"]
        gdir.mkdir(exist_ok=True)
        (gdir / "ref.txt").write_text(ref)
        rl = ref_len(FASTA[ref])
        for acc in present:
            if acc == ref:
                continue
            prefix = gdir / acc
            if not (Path(str(prefix) + ".report")).exists():
                run(f"nucmer --maxmatch -p {prefix} {FASTA[ref]} {FASTA[acc]} 2>/dev/null")
                run(f"dnadiff -d {prefix}.delta -p {prefix} 2>/dev/null")
            aligned, _ = aligned_bases(str(prefix) + ".report")
            snp, _ = count_snp(str(prefix) + ".snps")
            common = aligned if aligned else 0
            ratio = (common / rl) if rl else 0
            resolvable = (common >= 1_000_000) and (ratio >= 0.80)
            pair_rows.append({
                "group_id": g["group_id"], "group_type": g["group_type"],
                "reference": ref, "member": acc,
                "ref_length": rl, "common_aligned_bases": common,
                "common_ratio": round(ratio, 4),
                "pairwise_snp_ACGT": snp,
                "resolvable_by_snp": "yes" if resolvable else "CLONE_UNRESOLVED",
            })
    fields = ["group_id", "group_type", "reference", "member", "ref_length",
              "common_aligned_bases", "common_ratio", "pairwise_snp_ACGT", "resolvable_by_snp"]
    with open(OUT / "pairs_snp.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=fields)
        w.writeheader()
        w.writerows(pair_rows)
    n_un = sum(1 for p in pair_rows if p["resolvable_by_snp"] == "CLONE_UNRESOLVED")
    print(f"pairs: {len(pair_rows)} | unresolved: {n_un} -> {OUT/'pairs_snp.tsv'}")


if __name__ == "__main__":
    main()
