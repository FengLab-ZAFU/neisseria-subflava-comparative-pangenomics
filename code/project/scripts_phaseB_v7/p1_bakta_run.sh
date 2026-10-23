#!/usr/bin/env bash
# P1 step 3b: run Bakta 1.12.0 + DB v6.0 for the accessions in bakta_runlist.tsv
# (19 pending + any cache downgrades). Parameters identical to the historical
# scripts/revision_20261001/b1_bakta.sh (--compliant --force, prefix=accession,
# locus-tag=N<digits>); only threads differ (8 per genome, 19 parallel = 152
# cores <= 192).
set -u
source ${CONDA_SH:-$HOME/miniconda3/etc/profile.d/conda.sh}
conda activate bakta
ROOT=${NEISSERIA_PROJECT_ROOT:?Set NEISSERIA_PROJECT_ROOT}
OUT=$ROOT/results_revision_20261010_phaseB_v7/01_annotation
DB=${BAKTA_DB:?Set BAKTA_DB to the Bakta v6.0 full database}
RUNLIST=${1:-$OUT/bakta_runlist.tsv}
STATUS=${2:-$OUT/bakta_new_status.tsv}
mkdir -p "$OUT/bakta_new" "$OUT/bakta_new/logs"
PAR=19
THREADS=8

bakta --version
echo -e "accession\texit_code" > "$STATUS"

run_one() {
  acc="$1"; fasta="$2"
  od="$OUT/bakta_new/$acc"
  mkdir -p "$od"
  digits=$(echo "$acc" | grep -oE '[0-9]+' | head -1)
  ltag="N${digits}"
  if bakta --db "$DB" --threads $THREADS --prefix "$acc" --locus-tag "$ltag" \
      --compliant --force \
      --output "$od" "$fasta" > "$OUT/bakta_new/logs/$acc.log" 2>&1; then
    echo -e "$acc\t0" >> "$STATUS"
    echo "OK $acc"
  else
    rc=$?
    echo -e "$acc\t$rc" >> "$STATUS"
    echo "FAIL $acc rc=$rc"
  fi
}
export -f run_one; export OUT DB THREADS STATUS

date -u +"start_utc %Y-%m-%dT%H:%M:%SZ"
awk -F'\t' 'NR>1 && $2!=""{print $1"\t"$2}' "$RUNLIST" | \
  xargs -P $PAR -L 1 bash -c 'run_one "$@"' _ | tee "$OUT/bakta_new/logs/progress.log"
date -u +"end_utc %Y-%m-%dT%H:%M:%SZ"
echo "OK:   $(grep -c '^OK' "$OUT/bakta_new/logs/progress.log" || true)"
echo "FAIL: $(grep -c '^FAIL' "$OUT/bakta_new/logs/progress.log" || true)"
