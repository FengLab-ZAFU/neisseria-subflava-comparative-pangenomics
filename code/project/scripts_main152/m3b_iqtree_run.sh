#!/usr/bin/env bash
# M3b: main ML tree, main152 strict single-copy core concatenated alignment.
# Prespecified model (plan): GTR+F+I+R4, single partition, UFBoot 1000 + SH-aLRT 1000,
# seed 20261015, keep-ident. IQ-TREE 2.4.0 (conda env revision2026).
# Checkpoint resume enabled: rerun the same command to resume (no -redo).
set -u
source ${CONDA_SH:-$HOME/miniconda3/etc/profile.d/conda.sh}
set +u; conda activate gubbins; set -u
THREADS=${1:-64}
ROOT=${NEISSERIA_PROJECT_ROOT:?Set NEISSERIA_PROJECT_ROOT}
CORE=$ROOT/results_revision_20261015_main152_species/03_phylogeny/strict_singlecopy_core
OUT=$CORE/iqtree
mkdir -p "$OUT"
cd "$OUT"

{
date -u +"start_utc %Y-%m-%dT%H:%M:%SZ"
echo "threads=$THREADS seed=20261015"
echo "iqtree_version=$(iqtree2 --version 2>&1 | head -1)"
echo "cmd: iqtree2 -s ../core_alignment.fasta -m GTR+F+I+R4 -B 1000 -alrt 1000 -seed 20261015 -T $THREADS -keep-ident --prefix main152"
} > run_info.txt

set +e
iqtree2 -s ../core_alignment.fasta -m GTR+F+I+R4 \
  -B 1000 -alrt 1000 -seed 20261015 \
  -T $THREADS -keep-ident --prefix main152 \
  > iqtree_main152.stdout.log 2>&1
RC=$?
set -e
{
date -u +"end_utc %Y-%m-%dT%H:%M:%SZ"
echo "exit_code=$RC"
} >> run_info.txt
echo "iqtree exit=$RC"
exit $RC
