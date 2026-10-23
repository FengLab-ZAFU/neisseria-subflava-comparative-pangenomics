#!/usr/bin/env bash
# M7c: VFDB candidate homolog screen for the frozen main152 cohort.
# ABRicate 1.4.0, DB vfdb (4592 nucl, 2026-Apr-3), --minid 80 --mincov 80.
# Input: 152 frozen FASTA paths from 00_preflight/input_check/input_binding.tsv
# Output: raw per-genome TSV in 08_defense_mge_vfdb/vfdb/raw/<acc>.tsv
#         per-genome exit code ledger: 08_defense_mge_vfdb/vfdb/raw/run_exit_codes.tsv
set -u

PROJECT_ROOT=${NEISSERIA_PROJECT_ROOT:?Set NEISSERIA_PROJECT_ROOT}
RESULTS="$PROJECT_ROOT/results_revision_20261015_main152_species"
BINDING="$RESULTS/00_preflight/input_check/input_binding.tsv"
ORDER="$RESULTS/00_freeze/input_packet_snapshot/00_input/cohort/main152_accessions.txt"
OUTDIR="$RESULTS/08_defense_mge_vfdb/vfdb/raw"
mkdir -p "$OUTDIR"

source ${CONDA_SH:-$HOME/miniconda3/etc/profile.d/conda.sh}
conda activate abricate

EXIT_LEDGER="$OUTDIR/run_exit_codes.tsv"
printf 'accession\texit_code\tn_hit_lines\traw_tsv\n' > "$EXIT_LEDGER"

# Build accession -> path map in cohort order
mapfile -t ACCS < "$ORDER"

run_one() {
    acc="$1"
    path=$(awk -F'\t' -v a="$acc" 'NR>1 && $1==a {print $3; exit}' "$BINDING")
    raw="$OUTDIR/${acc}.tsv"
    if [ ! -f "$path" ]; then
        printf '%s\t%s\t%s\t%s\n' "$acc" "INPUT_MISSING" "NA" "$raw" >> "$EXIT_LEDGER.tmp.$acc"
        return
    fi
    abricate --db vfdb --minid 80 --mincov 80 "$path" > "$raw" 2> "$raw.stderr"
    ec=$?
    nh=$(awk -F'\t' 'NR>1 && NF>1' "$raw" | wc -l)
    printf '%s\t%s\t%s\t%s\n' "$acc" "$ec" "$nh" "$raw" >> "$EXIT_LEDGER.tmp.$acc"
    rm -f "$raw.stderr"
}
export -f run_one
export BINDING OUTDIR EXIT_LEDGER

printf '%s\n' "${ACCS[@]}" | xargs -P 12 -I{} bash -c 'run_one "$@"' _ {}

cat "$EXIT_LEDGER".tmp.* >> "$EXIT_LEDGER" 2>/dev/null
rm -f "$EXIT_LEDGER".tmp.*

echo "DONE: $(tail -n +2 "$EXIT_LEDGER" | wc -l) genomes processed"
echo "exit code distribution:"
tail -n +2 "$EXIT_LEDGER" | cut -f2 | sort | uniq -c
