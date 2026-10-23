#!/usr/bin/env bash
# A2 orchestrator: run after genus freeze completes.
# Chain: tier -> build urls -> download FASTAs -> bidirectional FastANI -> parse -> CheckM2 on candidates
set -u
cd ${NEISSERIA_PROJECT_ROOT:?Set NEISSERIA_PROJECT_ROOT}
echo "=== tiering ==="
python3 scripts/revision_20261001/a2_tier_screen.py || exit 1
echo "=== build download urls ==="
python3 scripts/revision_20261001/a2_build_download_list.py || exit 1
echo "=== download FASTAs (resumable) ==="
bash scripts/revision_20261001/a2_download_scope.sh || true
echo "=== bidirectional FastANI ==="
bash scripts/revision_20261001/a2_fastani_bidir.sh || true
echo "=== parse FastANI ==="
python3 scripts/revision_20261001/a2_parse_fastani.py || exit 1
echo "=== build ANI candidate list for QC ==="
python3 - <<'PY'
import csv
ani = list(csv.DictReader(open('results_revision_20261001/01_retrieval_taxonomy/ani_screen_results.tsv'), delimiter='\t'))
cands = [r['accession'] for r in ani if r['associated_group_member'] == 'True']
with open('results_revision_20261001/01_retrieval_taxonomy/ani_candidates.txt','w') as fh:
    fh.write('\n'.join(cands) + '\n')
print('ANI candidates for QC:', len(cands))
PY
echo "=== CheckM2 on candidates ==="
bash scripts/revision_20261001/a2_qc_checkm2.sh || true
echo "=== A2 chain DONE ==="
