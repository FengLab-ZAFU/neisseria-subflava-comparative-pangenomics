#!/usr/bin/env bash
# A2: parallel download of in-scope GenBank FASTAs from NCBI FTP.
# Usage: bash a2_download_scope.sh
# Resumable: skips accessions whose .fna already exists and is non-empty.
set -u
source ${CONDA_SH:-$HOME/miniconda3/etc/profile.d/conda.sh}
ROOT=${NEISSERIA_PROJECT_ROOT:?Set NEISSERIA_PROJECT_ROOT}
OUT=$ROOT/results_revision_20261001/01_retrieval_taxonomy
FASTA_DIR=$OUT/fasta
mkdir -p "$FASTA_DIR" "$OUT/logs"
PAR=48

dl_one() {
  acc="$1"; url="$2"
  fna="$FASTA_DIR/${acc}.fna"
  gz="$FASTA_DIR/${acc}.fna.gz"
  if [ -s "$fna" ]; then echo "SKIP $acc"; return 0; fi
  # retry up to 3 times
  for i in 1 2 3; do
    if curl -sfL --connect-timeout 20 --max-time 900 -o "$gz" "$url"; then
      if gunzip -f "$gz" 2>/dev/null && [ -s "$fna" ]; then
        echo "OK $acc"; return 0
      fi
    fi
    rm -f "$gz"; sleep $((i*5))
  done
  echo "FAIL $acc"; return 1
}
export -f dl_one
export FASTA_DIR

date -u +"start_utc %Y-%m-%dT%H:%M:%SZ" | tee "$OUT/logs/download_start.txt"
xargs -a "$OUT/download_urls.tsv" -P $PAR -L 1 bash -c 'dl_one "$@"' _ \
  > "$OUT/logs/download_progress.log" 2>&1
date -u +"end_utc %Y-%m-%dT%H:%M:%SZ" | tee "$OUT/logs/download_end.txt"
echo "OK:   $(grep -c '^OK' $OUT/logs/download_progress.log)"
echo "SKIP: $(grep -c '^SKIP' $OUT/logs/download_progress.log)"
echo "FAIL: $(grep -c '^FAIL' $OUT/logs/download_progress.log)"
grep '^FAIL' $OUT/logs/download_progress.log > "$OUT/download_failed.tsv" || true
ls "$FASTA_DIR"/*.fna | wc -l
