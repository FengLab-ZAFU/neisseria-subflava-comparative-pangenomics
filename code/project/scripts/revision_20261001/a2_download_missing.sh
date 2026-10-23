#!/bin/bash
# A2 download pass 2: smaller batches (50) for accessions still missing after pass 1.
# Usage: bash a2_download_missing.sh  (safe to re-run; skips existing fasta)
set -u
cd ${NEISSERIA_PROJECT_ROOT:?Set NEISSERIA_PROJECT_ROOT}
BASE=results_revision_20261001/01_retrieval_taxonomy
FASTA=$BASE/fasta
ZIPS=$BASE/zips2
LOGS=$BASE/logs/datasets_dl2
mkdir -p "$FASTA" "$ZIPS" "$LOGS"
LANES=${LANES:-6}
BATCH=50

source ${CONDA_SH:-$HOME/miniconda3/etc/profile.d/conda.sh}
conda activate ncbi-datasets

mapfile -t ACCS < <(tail -n +2 "$BASE/screening_scope_accessions.tsv" | cut -f1)
MISSING=()
for acc in "${ACCS[@]}"; do
  [ -s "$FASTA/$acc.fna" ] || MISSING+=("$acc")
done
TOTAL=${#MISSING[@]}
echo "Missing: $TOTAL"
[ "$TOTAL" -eq 0 ] && exit 0

dl_batch() { # batchfile
  local bf=$1
  local zname="b$(basename "$bf").zip"
  local zpath="$ZIPS/$zname"
  [ -s "$zpath.ok" ] && return 0
  local accs
  accs=$(paste -sd, "$bf")
  local attempt
  for attempt in 1 2 3 4 5; do
    if datasets download genome accession "$accs" --include genome --filename "$zpath" --no-progressbar > "$LOGS/$zname.log" 2>&1; then
      : > "$zpath.ok"; return 0
    fi
    sleep $((attempt*15))
  done
  echo "BATCHFAIL $zname"
  return 1
}
export -f dl_batch
export ZIPS LOGS

i=0; b=0
pids=()
while [ $i -lt $TOTAL ]; do
  bf="$BASE/.dl2_batch_$$.$b"
  printf '%s\n' "${MISSING[@]:$i:$BATCH}" > "$bf"
  dl_batch "$bf" > "$LOGS/lane_$b.status" 2>&1 &
  pids+=($!)
  if [ $(( (b+1) % LANES )) -eq 0 ] || [ $((i+BATCH)) -ge $TOTAL ]; then
    wait "${pids[@]}"; pids=()
  fi
  i=$((i+BATCH)); b=$((b+1))
done
rm -f "$BASE"/.dl2_batch_$$.*

# extract
cd "$ZIPS"
for z in b*.zip; do
  [ -f "$z" ] || continue
  tag=${z%.zip}
  [ -d "ex_$tag" ] && continue
  mkdir -p "ex_$tag"
  unzip -q -o "$z" -d "ex_$tag" 2>/dev/null
done
cd - >/dev/null

found=0
while IFS= read -r -d '' f; do
  acc=$(basename "$(dirname "$f")")
  [ -s "$FASTA/$acc.fna" ] && continue
  cp "$f" "$FASTA/$acc.fna"; found=$((found+1))
done < <(find "$ZIPS" -path "*/ncbi_dataset/data/*/*_genomic.fna" -print0)
echo "Newly extracted: $found"

left=0
> "$BASE/download_failed_final.tsv"
for acc in "${ACCS[@]}"; do
  [ -s "$FASTA/$acc.fna" ] || { echo "$acc" >> "$BASE/download_failed_final.tsv"; left=$((left+1)); }
done
echo "FINAL: missing=$left"
