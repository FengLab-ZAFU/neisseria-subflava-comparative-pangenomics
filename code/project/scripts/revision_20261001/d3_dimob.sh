#!/usr/bin/env bash
# D3 revision: Bakta GFF3 -> per-contig GBK -> IslandPath-DIMOB (per contig).
# Usage: bash d3_dimob.sh [n_jobs]
set -u
NJOBS=${1:-32}
ROOT=${NEISSERIA_PROJECT_ROOT:?Set NEISSERIA_PROJECT_ROOT}
OUT=$ROOT/results_revision_20261001/08_mge_gi/islandpath_dimob
LIST=${DIMOB_LIST:?Set DIMOB_LIST to the accession and FASTA list}
CONVERTER=${GFF3_TO_GBK:?Set GFF3_TO_GBK to the original per-contig conversion helper}
[ -f "$CONVERTER" ] || { echo "Missing GFF3_TO_GBK helper" >&2; exit 1; }
GFFDIR=$ROOT/results_revision_20261001/03_annotation/gff_all301

source ${CONDA_SH:-$HOME/miniconda3/etc/profile.d/conda.sh}
mkdir -p "$OUT/gbk_per_contig" "$OUT/islands" "$OUT/logs/convert"
export ROOT OUT GFFDIR CONVERTER

convert_one() {
    acc=$1
    source ${CONDA_SH:-$HOME/miniconda3/etc/profile.d/conda.sh}
    conda activate pgap2
    python3 "$CONVERTER" "$acc" \
        "$GFFDIR" \
        "$OUT/gbk_per_contig" \
        > "$OUT/logs/convert/$acc.json" 2> "$OUT/logs/convert/$acc.err"
    echo "$acc $?"
}
export -f convert_one
grep -v '^accession' "$LIST" | cut -d' ' -f1 | xargs -P 16 -I{} bash -c 'convert_one {}' > "$OUT/logs/convert_status.txt"
echo "GBK files: $(ls $OUT/gbk_per_contig/*.gbk 2>/dev/null | wc -l)"

dimob_one() {
    gbk=$1
    b=$(basename "$gbk" .gbk)
    [ -s "$OUT/islands/$b.gff" ] && return 0
    source ${CONDA_SH:-$HOME/miniconda3/etc/profile.d/conda.sh}
    conda activate islandpath
    islandpath "$gbk" "$OUT/islands/$b.gff" > /dev/null 2> "$OUT/logs/$b.stderr"
    rc=$?
    if [ $rc -ne 0 ] || [ ! -f "$OUT/islands/$b.gff" ]; then
        echo "FAIL $b rc=$rc" >> "$OUT/logs/dimob_failures.txt"
    fi
    return 0
}
export -f dimob_one
: > "$OUT/logs/dimob_failures.txt"
ls "$OUT/gbk_per_contig"/*.gbk | xargs -P "$NJOBS" -I{} bash -c 'dimob_one {}'
echo "DIMOB done: $(ls $OUT/islands/*.gff 2>/dev/null | wc -l) gff, $(grep -c . $OUT/logs/dimob_failures.txt 2>/dev/null || echo 0) failed"
