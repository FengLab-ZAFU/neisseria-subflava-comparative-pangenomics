#!/usr/bin/env bash
# fix3-H1 (Phase A / A4 reference calibration):
# The second N. perflava CCUG 17915 type-culture assembly (GCA_042649725.1,
# "perflava_type_dup") had no per-genome bidirectional FastANI results in the
# v3 delivery. Compute them here for all 1,303 screened genomes, both
# directions, with the same fastANI version/parameters as the v2 rescreen
# (fastANI 1.34, default k=16, frag=3000). Also completes ref-vs-ref pairs
# involving perflava_type_dup, and runs skani as the pre-agreed second method
# for the ANI_UNRESOLVED record GCA_051132615.1.
set -u
source ${CONDA_SH:-$HOME/miniconda3/etc/profile.d/conda.sh}
conda activate fastani
ROOT=${NEISSERIA_PROJECT_ROOT:?Set NEISSERIA_PROJECT_ROOT}
OLD=$ROOT/results_revision_20261001/01_retrieval_taxonomy
V2REFS=$ROOT/results_revision_20261006_fix/02_reference_check/refs
OUT=$ROOT/results_revision_20261009_fix3/01_reference_v4
RAW=$OUT/fastani_raw_dup
mkdir -p "$RAW" "$OUT/ref_vs_ref"
PAR=64
DUPREF=$V2REFS/GCA_042649725.1.fna
DUPNAME=perflava_type_dup

run_pair() {
  genome="$1"
  acc=$(basename "$genome" .fna)
  o1="$RAW/${acc}_vs_${DUPNAME}.txt"
  o2="$RAW/${DUPNAME}_vs_${acc}.txt"
  [ -s "$o1" ] || fastANI -q "$genome" -r "$DUPREF" -o "$o1" -t 1 2>"$RAW/${acc}_vs_${DUPNAME}.err" || echo "FAIL $acc g2r" >> "$RAW/failures.log"
  [ -s "$o2" ] || fastANI -q "$DUPREF" -r "$genome" -o "$o2" -t 1 2>"$RAW/${DUPNAME}_vs_${acc}.err" || echo "FAIL $acc r2q" >> "$RAW/failures.log"
  # no-hit pairs legitimately produce empty files; mark them so re-runs skip
  [ -f "$o1" ] || touch "$o1"
  [ -f "$o2" ] || touch "$o2"
}
export -f run_pair
export RAW DUPREF DUPNAME

date -u +"start_utc %Y-%m-%dT%H:%M:%SZ"
ls "$OLD/fasta"/*.fna | xargs -P $PAR -I{} bash -c 'run_pair "$@"' _ {}

# ref-vs-ref pairs involving perflava_type_dup (both directions)
declare -A R=( [subflava]=$V2REFS/N_subflava_ATCC49275.fna [flavescens]=$V2REFS/N_flavescens_ATCC13120.fna [perflava_type]=$V2REFS/GCA_023472875.1.fna [perflava_historical]=$V2REFS/N_perflava_27098_8_142_historical.fna [perflava_type_dup]=$DUPREF )
for a in perflava_type_dup; do
  for b in subflava flavescens perflava_type perflava_historical; do
    o="$OUT/ref_vs_ref/${a}_vs_${b}.txt"
    [ -s "$o" ] || fastANI -q "${R[$a]}" -r "${R[$b]}" -o "$o" -t 4 2>"$OUT/ref_vs_ref/${a}_vs_${b}.err" || echo "FAIL refpair $a $b" >> "$RAW/failures.log"
    o="$OUT/ref_vs_ref/${b}_vs_${a}.txt"
    [ -s "$o" ] || fastANI -q "${R[$b]}" -r "${R[$a]}" -o "$o" -t 4 2>"$OUT/ref_vs_ref/${b}_vs_${a}.err" || echo "FAIL refpair $b $a" >> "$RAW/failures.log"
  done
done
date -u +"fastani_done_utc %Y-%m-%dT%H:%M:%SZ"

# second method (skani) for the unresolved record vs all 5 references
SKANI=${SKANI:-skani}
if [ -x "$SKANI" ]; then
  G=$OLD/fasta/GCA_051132615.1.fna
  for b in subflava flavescens perflava_type perflava_type_dup perflava_historical; do
    "$SKANI" dist "$G" "${R[$b]}" > "$OUT/skani_GCA_051132615.1_vs_${b}.tsv" 2>> "$OUT/skani_GCA_051132615.1.log" || echo "skani $b failed" >> "$OUT/skani_GCA_051132615.1.log"
  done
fi
date -u +"end_utc %Y-%m-%dT%H:%M:%SZ"
echo "dup-ref result files: $(ls $RAW/*.txt 2>/dev/null | wc -l) (expect 2606)"
