#!/bin/bash
# M3b resume: continue UFBoot until convergence (acceptance A14), from checkpoint.
# Original run: iqtree2 -s ../core_alignment.fasta -m GTR+F+I+R4 -B 1000 -alrt 1000
#               -seed 20261015 -T 96 -keep-ident --prefix main152
# 2026-10-15 19:58 run ended at 1000 reps with correlation 0.944 (WARNING).
# This resume adds --nmax 3000 (default cap is 1000, hence the warning);
# checkpoint main152.ckp.gz is picked up automatically via --prefix.
# 1000-rep outputs preserved as *.b1000.
set -euo pipefail
THREADS=${1:-96}
DIR="$(cd "$(dirname "$0")/.." && pwd)/results_revision_20261015_main152_species/03_phylogeny/strict_singlecopy_core/iqtree"
cd "$DIR"
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate gubbins
iqtree2 -s ../core_alignment.fasta -m GTR+F+I+R4 -B 1000 -alrt 1000 \
  -seed 20261015 -T "$THREADS" -keep-ident --prefix main152 \
  --nmax 3000 --undo
echo "exit_code=$?" > run_info_resume.txt
date -u +%Y-%m-%dT%H:%M:%SZ >> run_info_resume.txt
