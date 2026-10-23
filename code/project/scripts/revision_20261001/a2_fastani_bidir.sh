#!/usr/bin/env bash
# A2: bidirectional FastANI of every downloaded genome vs the 3 type references.
# Rule (task doc A2): membership requires, per reference, bidirectional mean ANI >=95
# and min bidirectional aligned fraction >=0.70; AF per direction = matched/total query frags.
# Outputs raw per-direction results to fastani_raw/, parsed table built by a2_parse_fastani.py.
set -u
source ${CONDA_SH:-$HOME/miniconda3/etc/profile.d/conda.sh}
conda activate fastani
ROOT=${NEISSERIA_PROJECT_ROOT:?Set NEISSERIA_PROJECT_ROOT}
OUT=$ROOT/results_revision_20261001/01_retrieval_taxonomy
RAW=$OUT/fastani_raw
mkdir -p "$RAW"
PAR=64

REF_SUB=$ROOT/N_subflava_reference.fna
REF_PERF=$ROOT/N_perflava_reference.fna
REF_FLAV=$ROOT/N_flavescens_reference.fna

run_pair() {
  genome="$1"
  acc=$(basename "$genome" .fna)
  for refname in subflava perflava flavescens; do
    case $refname in
      subflava) ref=$REF_SUB;; perflava) ref=$REF_PERF;; flavescens) ref=$REF_FLAV;;
    esac
    o1="$RAW/${acc}_vs_${refname}.txt"   # genome as query
    o2="$RAW/${refname}_vs_${acc}.txt"   # ref as query
    [ -s "$o1" ] || fastANI -q "$genome" -r "$ref" -o "$o1" -t 1 2>"$RAW/${acc}_vs_${refname}.err" || true
    [ -s "$o2" ] || fastANI -q "$ref" -r "$genome" -o "$o2" -t 1 2>"$RAW/${refname}_vs_${acc}.err" || true
  done
  echo "done $acc"
}
export -f run_pair
export RAW REF_SUB REF_PERF REF_FLAV

date -u +"start_utc %Y-%m-%dT%H:%M:%SZ"
ls "$OUT/fasta"/*.fna | xargs -P $PAR -I{} bash -c 'run_pair "$@"' _ {}
date -u +"end_utc %Y-%m-%dT%H:%M:%SZ"
echo "raw result files: $(ls $RAW/*.txt | wc -l)"
