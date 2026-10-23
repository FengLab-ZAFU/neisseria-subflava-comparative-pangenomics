#!/usr/bin/env bash
# A2: parallel year-window chunk fetchers for the Neisseria genus freeze.
# Complements the main token-paginated fetcher (a2_fetch_genus_curl.sh); union + dedup
# by accession at merge time. Each window is independent (fresh token per window).
# Filter field 'last_release_date' verified via CLI --debug request body.
set -u
RAW=${NEISSERIA_PROJECT_ROOT:?Set NEISSERIA_PROJECT_ROOT}/results_revision_20261001/01_retrieval_taxonomy/raw
mkdir -p "$RAW/chunks"
URL="https://api.ncbi.nlm.nih.gov/datasets/v2/genome/dataset_report"
KEY="${NCBI_API_KEY:?set NCBI_API_KEY in the environment before running}"

fetch_window() {
  name="$1"; before="$2"   # records with release_date < before
  DEST="$RAW/chunks/window_${name}.jsonl"
  STATE="$RAW/chunks/window_${name}.state"
  TOK=""
  [ -f "$STATE" ] && TOK=$(cat "$STATE")
  while true; do
    BODY="{\"taxons\": [\"482\"], \"page_size\": 1000, \"filters\": {\"assembly_version\": \"current\", \"last_release_date\": \"${before}T00:00:00Z\"}"
    [ -n "$TOK" ] && BODY="$BODY, \"page_token\": \"$TOK\""
    BODY="$BODY}"
    ok=0
    for i in $(seq 1 15); do
      code=$(curl -s -m 75 --retry 2 --retry-delay 4 -X POST "$URL" \
        -H "Accept: application/json" -H "Content-Type: application/json" -H "api-key: $KEY" \
        -d "$BODY" -o "$RAW/chunks/.tmp_${name}.json" -w "%{http_code}")
      if [ "$code" = "200" ] && jq -e '.reports' "$RAW/chunks/.tmp_${name}.json" >/dev/null 2>&1; then ok=1; break; fi
      sleep $((i*6))
    done
    [ "$ok" = "1" ] || { echo "window $name PAGE FAILED"; return 1; }
    jq -c '.reports[]' "$RAW/chunks/.tmp_${name}.json" >> "$DEST"
    TOK=$(jq -r '.next_page_token // ""' "$RAW/chunks/.tmp_${name}.json")
    TOT=$(jq -r '.total_count // ""' "$RAW/chunks/.tmp_${name}.json")
    echo "$TOK" > "$STATE"
    n=$(wc -l < "$DEST")
    echo "window $name: $n/$TOT"
    [ -z "$TOK" ] && break
    [ "$n" -ge "${TOT:-0}" ] && [ -n "$TOT" ] && break
  done
  echo "window $name DONE: $(wc -l < $DEST)"
}

export -f fetch_window
export RAW URL KEY

# disjoint release-date windows (release_date < before)
windows=(
  "pre2012|2012-01-01"
  "y2012_2013|2014-01-01"
  "y2014|2015-01-01"
  "y2015|2016-01-01"
  "y2016|2017-01-01"
  "y2017|2018-01-01"
  "y2018|2019-01-01"
  "y2019|2020-01-01"
  "y2020|2021-01-01"
  "y2021|2022-01-01"
  "y2022|2023-01-01"
  "y2023|2024-01-01"
  "y2024|2025-01-01"
  "y2025|2026-01-01"
  "y2026|2027-01-01"
)

# run up to 4 windows concurrently to be kind to the flaky network
pids=()
for w in "${windows[@]}"; do
  name="${w%%|*}"; before="${w##*|}"
  fetch_window "$name" "$before" > "$RAW/chunks/window_${name}.log" 2>&1 &
  pids+=($!)
  while [ "$(jobs -rp | wc -l)" -ge 4 ]; do sleep 5; done
done
wait
echo "ALL WINDOWS DONE"
for w in "${windows[@]}"; do name="${w%%|*}"; echo "$name: $(wc -l < $RAW/chunks/window_${name}.jsonl 2>/dev/null || echo 0)"; done
