#!/usr/bin/env bash
# D2 revision: DefenseFinder on Bakta protein FASTA for all revision genomes.
# Usage: bash d2_defense.sh <list_file> [n_jobs]
set -u
LIST=${1:-/tmp/rgi_list.txt}
NJOBS=${2:-24}
ROOT=${NEISSERIA_PROJECT_ROOT:?Set NEISSERIA_PROJECT_ROOT}
OUT=$ROOT/results_revision_20261001/09_defense_systems/per_genome
DF=${DEFENSE_FINDER:-defense-finder}
mkdir -p "$OUT" "$ROOT/results_revision_20261001/09_defense_systems/logs"
export ROOT OUT DF

run_one() {
    acc=$1
    faa=$ROOT/results_revision_20261001/03_annotation/bakta/$acc/$acc.faa
    [ -s "$faa" ] || { echo "MISSING_FAA $acc" >> "$ROOT/results_revision_20261001/09_defense_systems/logs/df_failures.txt"; return 0; }
    [ -s "$OUT/$acc/${acc}_defense_finder_systems.tsv" ] && return 0
    $DF run --out-dir "$OUT/$acc" --workers 4 --coverage 0.4 "$faa" \
        > "$ROOT/results_revision_20261001/09_defense_systems/logs/$acc.log" 2>&1
    rc=$?
    if [ $rc -ne 0 ] || [ ! -s "$OUT/$acc/${acc}_defense_finder_systems.tsv" ]; then
        echo "FAIL $acc rc=$rc" >> "$ROOT/results_revision_20261001/09_defense_systems/logs/df_failures.txt"
    else
        echo "OK $acc" >> "$ROOT/results_revision_20261001/09_defense_systems/logs/df_progress.txt"
    fi
    return 0
}
export -f run_one
: > "$ROOT/results_revision_20261001/09_defense_systems/logs/df_progress.txt"
grep -v '^accession' "$LIST" | cut -d' ' -f1 | xargs -P "$NJOBS" -I{} bash -c 'run_one {}'
echo "defense-finder done: $(grep -c '^OK' $ROOT/results_revision_20261001/09_defense_systems/logs/df_progress.txt) OK, $(grep -c . $ROOT/results_revision_20261001/09_defense_systems/logs/df_failures.txt 2>/dev/null || echo 0) FAIL"
