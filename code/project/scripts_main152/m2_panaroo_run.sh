#!/usr/bin/env bash
# M2b: Panaroo strict method-comparison run for the main152 cohort.
# Mirrors historical invocation (panaroo 1.8.0):
#   panaroo -i <gff list> -o <out> --clean-mode strict --remove-invalid-genes -t <threads>
set -eu
source ${CONDA_SH:-$HOME/miniconda3/etc/profile.d/conda.sh}
set +u; conda activate panaroo; set -u
THREADS=$1
ROOT=${NEISSERIA_PROJECT_ROOT:?Set NEISSERIA_PROJECT_ROOT}
GFF_DIR=$ROOT/results_revision_20261015_main152_species/01_annotation/gff_main152
OUT=$ROOT/results_revision_20261015_main152_species/02_pangenome/panaroo_main152_strict
mkdir -p "$OUT"
LIST=$OUT/input_gff_list.txt
ls "$GFF_DIR"/*.gff3 > "$LIST"

{
date -u +"start_utc %Y-%m-%dT%H:%M:%SZ"
echo "gff_dir=$GFF_DIR n_gff=$(wc -l < "$LIST") threads=$THREADS"
echo "run_id=main152_20261015 seed=20261015"
echo "panaroo_version=$(panaroo --version 2>&1)"
echo "cmd: panaroo -i $LIST -o $OUT --clean-mode strict --remove-invalid-genes -t $THREADS"
} > "$OUT/run_info.txt"

set +e
panaroo -i "$LIST" -o "$OUT" --clean-mode strict --remove-invalid-genes -t "$THREADS" \
  > "$OUT/panaroo_main152_strict.log" 2>&1
RC=$?
set -e
{
date -u +"end_utc %Y-%m-%dT%H:%M:%SZ"
echo "exit_code=$RC"
} >> "$OUT/run_info.txt"
echo "panaroo exit=$RC"
exit $RC
