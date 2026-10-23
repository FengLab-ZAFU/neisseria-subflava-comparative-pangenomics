#!/usr/bin/env bash
# A2 download: per-accession datasets CLI downloads, highly parallel.
# Robust against server-side batch-zip assembly failures.
# Usage: bash a2_download_individual.sh [threads]
set -u
cd ${NEISSERIA_PROJECT_ROOT:?Set NEISSERIA_PROJECT_ROOT}
BASE=results_revision_20261001/01_retrieval_taxonomy
FASTA=$BASE/fasta
ZI=$BASE/zips_ind
LOGS=$BASE/logs/datasets_ind
mkdir -p "$FASTA" "$ZI" "$LOGS"
PAR=${1:-20}

source ${CONDA_SH:-$HOME/miniconda3/etc/profile.d/conda.sh}
conda activate ncbi-datasets

mapfile -t ACCS < <(tail -n +2 "$BASE/screening_scope_accessions.tsv" | cut -f1)
> "$BASE/.missing_ind.txt"
for acc in "${ACCS[@]}"; do
  [ -s "$FASTA/$acc.fna" ] || echo "$acc" >> "$BASE/.missing_ind.txt"
done
echo "missing: $(wc -l < "$BASE/.missing_ind.txt")"

export FASTA ZI LOGS BASE
dl_one() {
  acc=$1
  zp="$ZI/$acc.zip"
  for attempt in 1 2 3 4; do
    if datasets download genome accession "$acc" --include genome --filename "$zp" --no-progressbar > "$LOGS/$acc.log" 2>&1; then
      d="$ZI/ex_$acc"
      rm -rf "$d"; mkdir -p "$d"
      unzip -q -o "$zp" -d "$d" 2>/dev/null
      f=$(find "$d" -path "*/ncbi_dataset/data/*/*_genomic.fna" | head -1)
      if [ -n "$f" ] && [ ! -s "$FASTA/$acc.fna" ]; then
        cp "$f" "$FASTA/$acc.fna"
        echo "OK $acc"
        return 0
      fi
    fi
    sleep $((attempt*8))
  done
  echo "FAIL $acc"
  return 1
}
export -f dl_one

xargs -a "$BASE/.missing_ind.txt" -P "$PAR" -L 1 bash -c 'dl_one "$@"' _ > "$LOGS/progress.log" 2>&1

left=0
> "$BASE/download_failed_final.tsv"
for acc in "${ACCS[@]}"; do
  [ -s "$FASTA/$acc.fna" ] || { echo "$acc" >> "$BASE/download_failed_final.tsv"; left=$((left+1)); }
done
echo "DONE. fasta=$(ls "$FASTA" | wc -l) missing=$left"
grep -c "^OK" "$LOGS/progress.log" || true
grep "^FAIL" "$LOGS/progress.log" || true
