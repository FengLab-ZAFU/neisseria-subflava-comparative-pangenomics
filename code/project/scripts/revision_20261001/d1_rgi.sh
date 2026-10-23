#!/usr/bin/env bash
# D1: RGI (CARD, local DB) via Docker for a list of genomes - tet(M)/AMR detection.
# Usage: bash d1_rgi.sh <list_file> <out_subdir>
# list_file: one "analysis_id fasta_path" pair per line (TSV)
# Mirrors historical scripts/run_rgi.py exactly: rgi 6.0.3 image,
# db/rgi_local mounted at /carddb as the local DB workdir, -a DIAMOND.
set -u
LIST=$1
SUB=$2
ROOT=${NEISSERIA_PROJECT_ROOT:?Set NEISSERIA_PROJECT_ROOT}
OUT=$ROOT/results_revision_20261001/08_amr_distribution/$SUB
LOGS=$ROOT/results_revision_20261001/08_amr_distribution/logs
DB_DIR=$ROOT/db/rgi_local
mkdir -p "$OUT" "$LOGS"
IMG=quay.io/biocontainers/rgi:6.0.3--pyha8f3691_0
PAR=6

run_one() {
  aid="$1"; fasta="$2"
  [ -s "$OUT/$aid/$aid.txt" ] && { echo "SKIP $aid"; return 0; }
  mkdir -p "$OUT/$aid"
  rel=${fasta#${NEISSERIA_PROJECT_ROOT:?Set NEISSERIA_PROJECT_ROOT}/}
  if docker run --rm -v "$DB_DIR":/carddb -v "$ROOT":/work -w /carddb "$IMG" \
     rgi main --local -i "/work/$rel" -o "/work/results_revision_20261001/08_amr_distribution/$SUB/$aid/$aid" \
     -t contig -a DIAMOND -n 4 > "$LOGS/${aid}.log" 2>&1; then
    echo "OK $aid"
  else
    echo "FAIL $aid"
  fi
}
export -f run_one; export OUT ROOT IMG SUB LOGS DB_DIR

date -u +"start_utc %Y-%m-%dT%H:%M:%SZ"
xargs -a "$LIST" -P $PAR -L 1 bash -c 'run_one "$@"' _ | tee "$OUT/_progress.log"
date -u +"end_utc %Y-%m-%dT%H:%M:%SZ"
echo "OK: $(grep -c '^OK' $OUT/_progress.log) | SKIP: $(grep -c '^SKIP' $OUT/_progress.log) | FAIL: $(grep -c '^FAIL' $OUT/_progress.log)"
