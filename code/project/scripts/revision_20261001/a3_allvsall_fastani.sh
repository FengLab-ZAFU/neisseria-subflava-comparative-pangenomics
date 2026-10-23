#!/usr/bin/env bash
# A3: all-vs-all FastANI within cohort (for possible-duplicate flagging + distance distribution).
# Usage: bash a3_allvsall_fastani.sh <list_file> <outdir>
# list_file: one fasta path per line. Output: <outdir>/allvsall.tsv (query ref ani af)
set -eu
source ${CONDA_SH:-$HOME/miniconda3/etc/profile.d/conda.sh}
conda activate fastani
LIST=$1
OUTDIR=$2
mkdir -p "$OUTDIR"
THREADS=32
fastANI --queryList "$LIST" --refList "$LIST" --output "$OUTDIR/allvsall.tsv" \
  -t "$THREADS" 2> "$OUTDIR/allvsall.log"
echo "pairs: $(wc -l < $OUTDIR/allvsall.tsv)"
