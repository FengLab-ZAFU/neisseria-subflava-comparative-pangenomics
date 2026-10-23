#!/usr/bin/env bash
# A2: CheckM2 QC for ANI/coverage candidate genomes only (per task doc A2).
# Input: 01_retrieval_taxonomy/ani_candidates.txt (one accession per line)
# Output: 01_retrieval_taxonomy/checkm2/<acc>/ + merged quality_report
set -u
source ${CONDA_SH:-$HOME/miniconda3/etc/profile.d/conda.sh}
conda activate checkm2
ROOT=${NEISSERIA_PROJECT_ROOT:?Set NEISSERIA_PROJECT_ROOT}
OUT=$ROOT/results_revision_20261001/01_retrieval_taxonomy
CM=$OUT/checkm2
mkdir -p "$CM"
DB=${CHECKM2_DB:?Set CHECKM2_DB to uniref100.KO.1.dmnd}
PAR=24

run_one() {
  acc="$1"
  [ -s "$CM/$acc/quality_report.tsv" ] && { echo "SKIP $acc"; return 0; }
  fna=$OUT/fasta/$acc.fna
  mkdir -p "$CM/$acc"
  if checkm2 predict --threads 2 --input "$fna" --output-directory "$CM/$acc" \
      --force --database_path "$DB" > "$CM/$acc.log" 2>&1; then
    echo "OK $acc"
  else
    echo "FAIL $acc"
  fi
}
export -f run_one; export CM OUT DB

date -u +"start_utc %Y-%m-%dT%H:%M:%SZ"
xargs -a "$OUT/ani_candidates.txt" -P $PAR -L 1 bash -c 'run_one "$@"' _ | tee "$OUT/logs/checkm2_progress.log"
date -u +"end_utc %Y-%m-%dT%H:%M:%SZ"
# merge
python3 - <<'PY'
import csv, glob
rows = []
for d in sorted(glob.glob('${NEISSERIA_PROJECT_ROOT:?Set NEISSERIA_PROJECT_ROOT}/results_revision_20261001/01_retrieval_taxonomy/checkm2/*/quality_report.tsv')):
    acc = d.split('/')[-2]
    with open(d) as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            r['accession'] = acc
            rows.append(r)
if rows:
    fields = ['accession'] + [f for f in rows[0] if f != 'accession']
    with open('${NEISSERIA_PROJECT_ROOT:?Set NEISSERIA_PROJECT_ROOT}/results_revision_20261001/01_retrieval_taxonomy/checkm2_merged.tsv','w',newline='') as fh:
        w = csv.DictWriter(fh, delimiter='\t', fieldnames=fields); w.writeheader(); w.writerows(rows)
    print('merged', len(rows), 'rows')
PY
