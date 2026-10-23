#!/usr/bin/env bash
# A2: curl+jq paginated freeze of Neisseria genus (482) current assemblies.
# Resumable via state file. Appends JSONL pages to genus_482_current_frozen.jsonl.
set -u
RAW=${NEISSERIA_PROJECT_ROOT:?Set NEISSERIA_PROJECT_ROOT}/results_revision_20261001/01_retrieval_taxonomy/raw
DEST=$RAW/genus_482_current_frozen.jsonl
STATE=$RAW/genus_482_current_frozen.state.json
TMP=$RAW/.page_tmp.json
URL="https://api.ncbi.nlm.nih.gov/datasets/v2/genome/dataset_report"
KEY="${NCBI_API_KEY:?set NCBI_API_KEY in the environment before running}"

TOK=$(python3 -c "import json;print(json.load(open('$STATE'))['page_token'])" 2>/dev/null || echo "")
while true; do
  BODY="{\"taxons\": [\"482\"], \"page_size\": 1000, \"filters\": {\"assembly_version\": \"current\"}"
  [ -n "$TOK" ] && BODY="$BODY, \"page_token\": \"$TOK\""
  BODY="$BODY}"
  ok=0
  for i in 1 2 3 4 5 6 7 8 9 10; do
    code=$(curl -s -m 90 --retry 3 --retry-delay 5 -X POST "$URL" \
      -H "Accept: application/json" -H "Content-Type: application/json" -H "api-key: $KEY" \
      -d "$BODY" -o "$TMP" -w "%{http_code}")
    if [ "$code" = "200" ] && jq -e '.reports' "$TMP" >/dev/null 2>&1; then ok=1; break; fi
    echo "retry $i (http=$code)"; sleep $((i*8))
  done
  [ "$ok" = "1" ] || { echo "PAGE FAILED"; exit 1; }
  n=$(jq '.reports | length' "$TMP")
  jq -c '.reports[]' "$TMP" >> "$DEST"
  TOK=$(jq -r '.next_page_token // ""' "$TMP")
  tot=$(jq -r '.total_count // ""' "$TMP")
  echo "$(date -u +%H:%M:%S) +$n total_lines=$(wc -l < $DEST)/$tot"
  python3 -c "import json;json.dump({'page_token':'$TOK','fetched':$(wc -l < $DEST)},open('$STATE','w'))"
  [ -z "$TOK" ] && break
  [ "$n" = "0" ] && break
done
echo "FREEZE DONE: $(wc -l < $DEST) records"
