#!/bin/bash
# A2 download of screening-scope genomes via NCBI datasets CLI (FTP was 503-throttled).
# Input: 01_retrieval_taxonomy/screening_scope_accessions.tsv (col1 = GCA accession)
# Output: 01_retrieval_taxonomy/fasta/<acc>.fna ; logs per batch; retry pass for failures.
set -u
cd ${NEISSERIA_PROJECT_ROOT:?Set NEISSERIA_PROJECT_ROOT}
BASE=results_revision_20261001/01_retrieval_taxonomy
SCOPE=$BASE/screening_scope_accessions.tsv
FASTA=$BASE/fasta
ZIPS=$BASE/zips
LOGS=$BASE/logs/datasets_dl
mkdir -p "$FASTA" "$ZIPS" "$LOGS"
LANES=${LANES:-5}
BATCH=200

source ${CONDA_SH:-$HOME/miniconda3/etc/profile.d/conda.sh}
conda activate ncbi-datasets

mapfile -t ACCS < <(tail -n +2 "$SCOPE" | cut -f1)
TOTAL=${#ACCS[@]}
echo "Total accessions: $TOTAL"

download_batch() { # $1=batchfile $2=laneid
  local batchfile=$1 lane=$2
  local n zname zpath
  n=$(wc -l < "$batchfile")
  zname="batch_${batchfile##*/}.zip"
  zpath="$ZIPS/$zname"
  if [ -s "$zpath.ok" ]; then echo "[lane$lane] skip done $zname"; return 0; fi
  local accs
  accs=$(paste -sd, "$batchfile")
  local attempt
  for attempt in 1 2 3; do
    if datasets download genome accession "$accs" --include genome --filename "$zpath" --no-progressbar > "$LOGS/${zname}.log" 2>&1; then
      : > "$zpath.ok"
      echo "[lane$lane] OK $zname ($n acc)"
      return 0
    fi
    echo "[lane$lane] attempt $attempt failed for $zname; sleeping"
    sleep $((attempt*20))
  done
  echo "[lane$lane] FAILED $zname"
  return 1
}

# --- pass 1: batch download in parallel lanes ---
i=0; b=0
pids=()
while [ $i -lt $TOTAL ]; do
  bf="$BASE/.dl_batch_$$.$b"
  printf '%s\n' "${ACCS[@]:$i:$BATCH}" > "$bf"
  lane=$(( b % LANES ))
  download_batch "$bf" "$lane" > "$LOGS/lane_${b}.status" 2>&1 &
  pids+=($!)
  if [ $(( (b+1) % LANES )) -eq 0 ] || [ $((i+BATCH)) -ge $TOTAL ]; then
    wait "${pids[@]}"
    pids=()
  fi
  i=$((i+BATCH)); b=$((b+1))
done
rm -f "$BASE"/.dl_batch_$$.*

# --- unzip everything ---
cd "$ZIPS"
for z in batch_*.zip; do
  [ -f "$z" ] || continue
  tag=${z%.zip}
  [ -d "ex_$tag" ] && continue
  mkdir -p "ex_$tag"
  unzip -q -o "$z" -d "ex_$tag" 2>/dev/null
done
cd - >/dev/null

# collect fna -> fasta/<acc>.fna
found=0
while IFS= read -r -d '' f; do
  acc=$(basename "$(dirname "$f")")
  [ -s "$FASTA/$acc.fna" ] && continue
  cp "$f" "$FASTA/$acc.fna"
  found=$((found+1))
done < <(find "$ZIPS" -path "*/ncbi_dataset/data/*/*_genomic.fna" -print0)
echo "Newly extracted: $found"

# --- pass 2: individual retry for any missing accession ---
missing=0
> "$BASE/download_failed2.tsv"
for acc in "${ACCS[@]}"; do
  if [ ! -s "$FASTA/$acc.fna" ]; then
    echo "$acc" >> "$BASE/download_failed2.tsv"
    missing=$((missing+1))
  fi
done
echo "After batch pass: missing=$missing (listed in download_failed2.tsv)"

if [ "$missing" -gt 0 ] && [ "${INDIV_RETRY:-1}" = "1" ]; then
  echo "Individual retry pass..."
  retry_one() {
    local acc=$1
    local zp="$ZIPS/ind_$acc.zip"
    for attempt in 1 2 3; do
      if datasets download genome accession "$acc" --include genome --filename "$zp" --no-progressbar > "$LOGS/ind_$acc.log" 2>&1; then
        local d="$ZIPS/ex_ind_$acc"
        mkdir -p "$d" && unzip -q -o "$zp" -d "$d" 2>/dev/null
        local f
        f=$(find "$d" -path "*/ncbi_dataset/data/*/*_genomic.fna" | head -1)
        if [ -n "$f" ]; then cp "$f" "$FASTA/$acc.fna"; echo "RETRY-OK $acc"; return 0; fi
      fi
      sleep $((attempt*15))
    done
    echo "RETRY-FAIL $acc"
    return 1
  }
  export -f retry_one
  export ZIPS FASTA
  # modest parallelism for individual retries
  xargs -a "$BASE/download_failed2.tsv" -P 4 -I{} bash -c 'retry_one {}' >> "$LOGS/ind_retry.status" 2>&1
fi

final_missing=0
> "$BASE/download_failed_final.tsv"
for acc in "${ACCS[@]}"; do
  [ -s "$FASTA/$acc.fna" ] || { echo "$acc" >> "$BASE/download_failed_final.tsv"; final_missing=$((final_missing+1)); }
done
echo "FINAL: total=$TOTAL downloaded=$((TOTAL-final_missing)) missing=$final_missing"
