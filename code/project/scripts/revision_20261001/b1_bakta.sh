#!/usr/bin/env bash
# B1: uniform Bakta annotation for frozen representative genomes.
# Input: 03_annotation/bakta_input.tsv (analysis_id<TAB>fasta_path)
# Output: 03_annotation/bakta/<analysis_id>/ (GFF3/GBFF/FAA/FFN/TSV)
# Method note: historical Bakta version unverifiable (no version string in old logs)
# -> all representatives re-annotated with bakta 1.12.0 + DB v6.0 (2025-02-24), documented.
set -u
source ${CONDA_SH:-$HOME/miniconda3/etc/profile.d/conda.sh}
conda activate bakta
ROOT=${NEISSERIA_PROJECT_ROOT:?Set NEISSERIA_PROJECT_ROOT}
OUT=$ROOT/results_revision_20261001/03_annotation
DB=${BAKTA_DB:?Set BAKTA_DB to the Bakta v6.0 full database}
mkdir -p "$OUT/bakta" "$OUT/logs"
PAR=32
THREADS=4

run_one() {
  aid="$1"; fasta="$2"
  od="$OUT/bakta/$aid"
  [ -s "$od/$aid.gff3" ] && [ -s "$od/$aid.faa" ] && { echo "SKIP $aid"; return 0; }
  mkdir -p "$od"
  # INSDC-compliant deterministic locus tag: N + digits of assembly number (<=12 alnum, starts with letter)
  digits=$(echo "$aid" | grep -oE '[0-9]+' | head -1)
  ltag="N${digits}"
  if bakta --db "$DB" --threads $THREADS --prefix "$aid" --locus-tag "$ltag" \
      --compliant --force \
      --output "$od" "$fasta" > "$OUT/logs/$aid.log" 2>&1; then
    echo "OK $aid"
  else
    echo "FAIL $aid"
  fi
}
export -f run_one; export OUT DB THREADS

date -u +"start_utc %Y-%m-%dT%H:%M:%SZ"
tail -n +2 "$OUT/bakta_input.tsv" | xargs -P $PAR -L 1 bash -c 'run_one "$@"' _ | tee "$OUT/logs/bakta_progress.log"
date -u +"end_utc %Y-%m-%dT%H:%M:%SZ"
echo "OK:   $(grep -c '^OK' $OUT/logs/bakta_progress.log)"
echo "SKIP: $(grep -c '^SKIP' $OUT/logs/bakta_progress.log)"
echo "FAIL: $(grep -c '^FAIL' $OUT/logs/bakta_progress.log)"
grep '^FAIL' "$OUT/logs/bakta_progress.log" > "$OUT/bakta_failed.tsv" || true
