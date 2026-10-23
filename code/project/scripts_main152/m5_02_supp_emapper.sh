#!/usr/bin/env bash
# M5 supplementary eggNOG annotation for the 2,854 main152 protein sequences
# absent from the historical emapper-3.0.0-beta6 annotation bridge.
# Same software (emapper-3.0.0-beta6 SIF) and same DB snapshot (eggnog DB 7.0.0,
# supplied through EGGNOG_DB); <=16 CPU.
set -eu
ROOT=${NEISSERIA_PROJECT_ROOT:?Set NEISSERIA_PROJECT_ROOT}
OUT=$ROOT/results_revision_20261015_main152_species/06_function/supplement_emapper
SIF=${EGGNOG_SIF:?Set EGGNOG_SIF to the emapper-3.0.0-beta6 image}
EGGNOG_DB=${EGGNOG_DB:?Set EGGNOG_DB to eggNOG database 7.0.0}
APPTAINER=${APPTAINER:-apptainer}
mkdir -p "$OUT"
date -u +"start_utc %Y-%m-%dT%H:%M:%SZ" | tee "$OUT/run_info.txt"
$APPTAINER exec \
  --bind "$ROOT:/work" \
  --bind "$EGGNOG_DB:/data" \
  --bind /tmp:/tmp \
  "$SIF" /usr/local/bin/emapper.py \
  -i /work/results_revision_20261015_main152_species/06_function/unmatched_proteins.fasta \
  -o m5_supp --itype proteins --data_dir /data \
  --output_dir /work/results_revision_20261015_main152_species/06_function/supplement_emapper \
  --cpu 16 --temp_dir /tmp --dmnd_block_size 8 --dmnd_index_chunks 1 \
  > "$OUT/emapper.log" 2>&1
date -u +"end_utc %Y-%m-%dT%H:%M:%SZ" | tee -a "$OUT/run_info.txt"
echo "exit_code=$?" | tee -a "$OUT/run_info.txt"
ls -la "$OUT"/*.annotations
