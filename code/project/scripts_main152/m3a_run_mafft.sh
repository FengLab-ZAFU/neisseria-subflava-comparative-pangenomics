#!/usr/bin/env bash
# M3a: per-gene MAFFT nucleotide alignments (parallel).
# Usage: bash m3a_run_mafft.sh <parallel_jobs>
set -u
source ${CONDA_SH:-$HOME/miniconda3/etc/profile.d/conda.sh}
set +u; conda activate pgap2; set -u
P=${1:-32}
ROOT=${NEISSERIA_PROJECT_ROOT:?Set NEISSERIA_PROJECT_ROOT}
CORE=$ROOT/results_revision_20261015_main152_species/03_phylogeny/strict_singlecopy_core
IN=$CORE/per_gene_fastas
OUTDIR=$CORE/per_gene_alignments
LOGD=$CORE/mafft_logs
mkdir -p "$OUTDIR" "$LOGD"
{
date -u +"start_utc %Y-%m-%dT%H:%M:%SZ"
echo "mafft_version=$(mafft --version 2>&1 | head -1)"
echo "parallel_jobs=$P"
echo "cmd: mafft --auto --thread 1 --inputorder <in> > <out>  (per gene)"
} > "$CORE/mafft_run_info.txt"

ls "$IN"/*.fna | xargs -P "$P" -I {} bash -c '
  b=$(basename {} .fna)
  mafft --auto --thread 1 --inputorder "{}" > "'"$OUTDIR"'/$b.aln.fna" 2> "'"$LOGD"'/$b.log"
  echo "$b $?"
' > "$CORE/mafft_exit_codes.txt"

{
date -u +"end_utc %Y-%m-%dT%H:%M:%SZ"
echo "n_alignments=$(ls $OUTDIR/*.aln.fna 2>/dev/null | wc -l)"
echo "n_failed=$(awk '$2!=0' $CORE/mafft_exit_codes.txt | wc -l)"
} >> "$CORE/mafft_run_info.txt"
awk '$2!=0' "$CORE/mafft_exit_codes.txt" | head
echo "MAFFT done: $(ls $OUTDIR/*.aln.fna | wc -l) alignments, $(awk '$2!=0' $CORE/mafft_exit_codes.txt | wc -l) failed"
