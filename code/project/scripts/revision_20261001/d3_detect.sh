#!/usr/bin/env bash
# D3: MGE detection on the revision annotation set (fasta-only tools).
# Tools: geNomad (end-to-end), ISEScan, IntegronFinder, MOB-suite (mob_recon).
# Same invocation style as historical runners (scripts/run_*.py), recorded per tool.
set -u
cd ${NEISSERIA_PROJECT_ROOT:?Set NEISSERIA_PROJECT_ROOT}
BASE=results_revision_20261001/08_mge_gi
INPUT=results_revision_20261001/03_annotation/bakta_input.tsv
mkdir -p "$BASE/genomad" "$BASE/isescan" "$BASE/integron_finder" "$BASE/mob_suite" "$BASE/logs"
GDB=${GENOMAD_DB:?Set GENOMAD_DB}
ISE=${ISESCAN_SCRIPT:?Set ISESCAN_SCRIPT to isescan.py}

run_genomad() {
  aid="$1"; fa="$2"
  od="$BASE/genomad/$aid"
  [ -s "$od/${aid}_summary/${aid}_virus_summary.tsv" ] && { echo "SKIP $aid"; return 0; }
  rm -rf "$od"
  if genomad end-to-end --threads 8 --cleanup "$fa" "$od" "$GDB" > "$BASE/logs/genomad_$aid.log" 2>&1; then
    echo "OK $aid"; else echo "FAIL $aid"; fi
}
run_isescan() {
  aid="$1"; fa="$2"
  od="$BASE/isescan/$aid"
  [ -n "$(ls "$od"/*.tsv 2>/dev/null)" ] && { echo "SKIP $aid"; return 0; }
  rm -rf "$od"; mkdir -p "$od"
  if python3 "$ISE" --seqfile "$fa" --output "$od" --nthread 2 > "$BASE/logs/isescan_$aid.log" 2>&1; then
    echo "OK $aid"; else echo "FAIL $aid"; fi
}
run_integron() {
  aid="$1"; fa="$2"
  od="$BASE/integron_finder/$aid"
  [ -n "$(find "$od" -name '*.integrons' 2>/dev/null)" ] && { echo "SKIP $aid"; return 0; }
  rm -rf "$od"
  if integron_finder --cpu 2 --outdir "$od" "$fa" > "$BASE/logs/integron_$aid.log" 2>&1; then
    echo "OK $aid"; else echo "FAIL $aid"; fi
}
run_mob() {
  aid="$1"; fa="$2"
  od="$BASE/mob_suite/$aid"
  [ -s "$od/mobtyper_results.txt" ] && { echo "SKIP $aid"; return 0; }
  rm -rf "$od"
  if mob_recon -i "$fa" -o "$od" -n 4 -s "$aid" -f > "$BASE/logs/mob_$aid.log" 2>&1; then
    echo "OK $aid"; else echo "FAIL $aid"; fi
}
export -f run_genomad run_isescan run_integron run_mob
export BASE GDB ISE

tail -n +2 "$INPUT" > "$BASE/.d3_input.tsv"
N=$(wc -l < "$BASE/.d3_input.tsv")

echo "== geNomad ($N genomes, 8 workers x 8 threads) =="
(
  source ${CONDA_SH:-$HOME/miniconda3/etc/profile.d/conda.sh} && conda activate genomad
  xargs -a "$BASE/.d3_input.tsv" -P 8 -L 1 bash -c 'run_genomad "$@"' _ > "$BASE/logs/genomad_progress.log" 2>&1
) &
PID_G=$!

echo "== ISEScan =="
(
  source ${CONDA_SH:-$HOME/miniconda3/etc/profile.d/conda.sh} && conda activate isescan
  xargs -a "$BASE/.d3_input.tsv" -P 24 -L 1 bash -c 'run_isescan "$@"' _ > "$BASE/logs/isescan_progress.log" 2>&1
) &
PID_I=$!

echo "== IntegronFinder =="
(
  source ${CONDA_SH:-$HOME/miniconda3/etc/profile.d/conda.sh} && conda activate integron_finder
  xargs -a "$BASE/.d3_input.tsv" -P 16 -L 1 bash -c 'run_integron "$@"' _ > "$BASE/logs/integron_progress.log" 2>&1
) &
PID_IN=$!

echo "== MOB-suite =="
(
  source ${CONDA_SH:-$HOME/miniconda3/etc/profile.d/conda.sh} && conda activate mob_suite
  xargs -a "$BASE/.d3_input.tsv" -P 12 -L 1 bash -c 'run_mob "$@"' _ > "$BASE/logs/mob_progress.log" 2>&1
) &
PID_M=$!

wait $PID_G $PID_I $PID_IN $PID_M
for t in genomad isescan integron mob; do
  echo "$t: OK=$(grep -c '^OK' $BASE/logs/${t}_progress.log) SKIP=$(grep -c '^SKIP' $BASE/logs/${t}_progress.log) FAIL=$(grep -c '^FAIL' $BASE/logs/${t}_progress.log)"
done
