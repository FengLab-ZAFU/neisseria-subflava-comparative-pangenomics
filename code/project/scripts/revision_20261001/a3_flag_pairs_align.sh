#!/usr/bin/env bash
# A3: parallel nucmer/dnadiff for flag_pairs_list.tsv, then summarize to flag_pairs_snp.tsv.
# Usage: bash a3_flag_pairs_align.sh <pairs_list.tsv> <out_tsv> [threads]
set -u
source ${CONDA_SH:-$HOME/miniconda3/etc/profile.d/conda.sh}
conda activate mummer_syri
ROOT=${NEISSERIA_PROJECT_ROOT:?Set NEISSERIA_PROJECT_ROOT}
OUT=$ROOT/results_revision_20261001/02_independence
LIST=${1:-$OUT/flag_pairs_list.tsv}
TSV=${2:-$OUT/flag_pairs_snp.tsv}
PAR=${3:-24}
ALN=$OUT/alignments_screen
mkdir -p "$ALN"

# fasta lookup: retrieval dir (GCA) + Genomes/ with assembly->biosample alias
export ROOT OUT ALN
python3 - "$LIST" > "$ALN/.pairs_with_paths.tsv" <<'PY'
import csv, sys
from pathlib import Path
ROOT = Path("${NEISSERIA_PROJECT_ROOT:?Set NEISSERIA_PROJECT_ROOT}")
FASTA = {}
for d in (ROOT/"results_revision_20261001/01_retrieval_taxonomy/fasta", ROOT/"Genomes"):
    for p in d.glob("*.fna"):
        FASTA[p.stem] = p
with open(ROOT/"neisseria_subflava_associated_group_metadata.tsv", encoding="utf-8-sig") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        acc=(r.get("assembly_accession") or "").strip(); fn=(r.get("fasta_filename") or "").strip()
        if acc and fn and (ROOT/"Genomes"/fn).exists():
            FASTA.setdefault(acc, ROOT/"Genomes"/fn)
for line in open(sys.argv[1]):
    a,b = line.split()
    fa, fb = FASTA.get(a), FASTA.get(b)
    if fa and fb:
        print(f"{a}\t{b}\t{fa}\t{fb}")
    else:
        print(f"{a}\t{b}\tMISSING\tMISSING", file=sys.stderr)
PY

align_one() {
  a=$1; b=$2; fa=$3; fb=$4
  d="$ALN/${a}__${b}"
  mkdir -p "$d"
  [ -s "$d/aln.report" ] && return 0
  nucmer --maxmatch -p "$d/aln" "$fa" "$fb" 2>/dev/null
  dnadiff -d "$d/aln.delta" -p "$d/aln" 2>/dev/null
}
export -f align_one

xargs -a "$ALN/.pairs_with_paths.tsv" -P "$PAR" -L 1 bash -c 'align_one "$@"' _

# summarize
python3 - "$ALN" "$TSV" <<'PY'
import csv, re, sys
from pathlib import Path
ALN, TSV = Path(sys.argv[1]), Path(sys.argv[2])
ACGT = set("ACGT")
def aligned_bases(rp):
    try: txt = Path(rp).read_text()
    except Exception: return None
    m = re.search(r"AlignedBases\s+(\d+)\S*\s+(\d+)", txt)
    return int(m.group(1)) if m else None
def count_snp(sp):
    n=0
    if not Path(sp).exists(): return 0
    for line in open(sp):
        f=line.split("\t")
        if len(f)<3: continue
        a,b=f[0].strip().upper(),f[1].strip().upper()
        if a in ACGT and b in ACGT: n+=1
    return n
def ref_len(fa):
    n=0
    for line in open(fa):
        if not line.startswith(">"): n+=len(line.strip())
    return n
# recover ani/af per pair from list + allvsall
import collections
best = {}
av = TSV.parent/"allvsall_189"/"allvsall.tsv"
for line in open(av):
    q,r,ani,m,t = line.split("\t")[:5]
    if q==r: continue
    try: ani_f, af = float(ani), int(m)/int(t)
    except: continue
    if ani_f < 99.5 or af < 0.8: continue
    a,b = sorted((Path(q).stem, Path(r).stem))
    cur = best.get((a,b))
    if cur is None: best[(a,b)]=[ani_f,af]
    else: cur[0]=min(cur[0],ani_f); cur[1]=min(cur[1],af)
rows=[]
for line in open(ALN/".pairs_with_paths.tsv"):
    a,b,fa,fb = line.split()
    if fa=="MISSING":
        rows.append({"pair":f"{a};{b}","ani_min":"","af_min":"","ref_length":"","common_aligned_bases":"",
                     "common_ratio":"","pairwise_snp_ACGT":"","resolvable_by_snp":"","status":"FASTA_MISSING"})
        continue
    d = ALN/f"{a}__{b}"
    rl = ref_len(fa)
    common = aligned_bases(d/"aln.report") or 0
    snp = count_snp(d/"aln.snps")
    ratio = common/rl if rl else 0
    res = common>=1_000_000 and ratio>=0.80
    ani_m, af_m = best.get((a,b), ["",""])
    rows.append({"pair":f"{a};{b}","ani_min":round(ani_m,4) if ani_m!="" else "",
                 "af_min":round(af_m,4) if af_m!="" else "",
                 "ref_length":rl,"common_aligned_bases":common,"common_ratio":round(ratio,4),
                 "pairwise_snp_ACGT":snp,"resolvable_by_snp":"yes" if res else "CLONE_UNRESOLVED","status":"OK"})
fields=["pair","ani_min","af_min","ref_length","common_aligned_bases","common_ratio","pairwise_snp_ACGT","resolvable_by_snp","status"]
with open(TSV,"w",newline="") as fh:
    w=csv.DictWriter(fh,delimiter="\t",fieldnames=fields); w.writeheader(); w.writerows(rows)
print(f"wrote {len(rows)} rows -> {TSV}")
PY
