#!/usr/bin/env bash
# M2: PGAP2 pangenome run for the main152 cohort (2026-10-15 run).
# Usage: bash m2_pgap2_run.sh <run_name> <threads> [orth_id]
#   run_name: main152 (orth_id 0.98) | main152_orth95 (0.95)
# Settings mirror the historically verified run (pgap2 2.4.1):
#   --type-filter CDS --id-attr-key ID --aligner diamond --clust_method mmseqs2
set -eu
source ${CONDA_SH:-$HOME/miniconda3/etc/profile.d/conda.sh}
set +u; conda activate pgap2; set -u
RUN=$1
THREADS=$2
ORTH=${3:-0.98}
ROOT=${NEISSERIA_PROJECT_ROOT:?Set NEISSERIA_PROJECT_ROOT}
GFF_DIR=$ROOT/results_revision_20261015_main152_species/01_annotation/gff_main152
OUT=$ROOT/results_revision_20261015_main152_species/02_pangenome/pgap2_${RUN}
mkdir -p "$OUT"

{
date -u +"start_utc %Y-%m-%dT%H:%M:%SZ"
echo "run=$RUN orth_id=$ORTH gff_dir=$GFF_DIR n_gff=$(ls $GFF_DIR/*.gff3 | wc -l) threads=$THREADS"
echo "run_id=main152_20261015 seed=20261015"
echo "pgap2_version=$(pgap2 --version 2>&1)"
echo "diamond_version=$(diamond version 2>&1)"
echo "mmseqs_version=$(mmseqs version 2>&1)"
echo "cmd: pgap2 main -i $GFF_DIR -o $OUT -t $THREADS --orth_id $ORTH --type-filter CDS --id-attr-key ID --aligner diamond --clust_method mmseqs2"
} > "$OUT/run_info.txt"

set +e
pgap2 main -i "$GFF_DIR" -o "$OUT" -t "$THREADS" \
  --orth_id "$ORTH" \
  --type-filter CDS --id-attr-key ID --aligner diamond --clust_method mmseqs2 \
  > "$OUT/pgap2_${RUN}.log" 2>&1
RC=$?
set -e
{
date -u +"end_utc %Y-%m-%dT%H:%M:%SZ"
echo "exit_code=$RC"
} >> "$OUT/run_info.txt"
echo "$RUN exit=$RC"
exit $RC
