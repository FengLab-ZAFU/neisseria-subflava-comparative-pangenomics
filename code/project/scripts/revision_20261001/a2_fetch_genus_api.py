#!/usr/bin/env python3
"""A2: resumable paginated freeze of Neisseria genus (taxid 482) current assembly metadata.

The datasets CLI attempt (recorded in genus_482_current_unfiltered.query_info.txt)
failed with 'connection reset by peer' at ~22000/unknown records. This fetcher uses
the documented NCBI Datasets v2 API with the same filters (assembly_version=current
only; no source/MAG/atypical filters) and per-page retries + resume.

Output: raw/genus_482_current_frozen.jsonl (+ state file for resume)
"""
import json
import time
import urllib.request
import urllib.error
from pathlib import Path

OUT = Path(__file__).resolve().parent.parent.parent / "results_revision_20261001" / "01_retrieval_taxonomy" / "raw"
OUT.mkdir(parents=True, exist_ok=True)
import os
DEST = Path(os.environ.get("FETCH_DEST", OUT / "genus_482_current_frozen.jsonl"))
STATE = Path(os.environ.get("FETCH_STATE", OUT / "genus_482_current_frozen.state.json"))
LAST_RELEASE = os.environ.get("FETCH_LAST_RELEASE", "")
URL = "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/dataset_report"
API_KEY = os.environ["NCBI_API_KEY"]  # set NCBI_API_KEY in the environment; no credentials in code
PAGE = 1000

state = {"page_token": None, "fetched": 0, "started_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
if STATE.exists():
    state = json.loads(STATE.read_text())
    print(f"resuming: {state['fetched']} records already fetched")


def post(body, tries=40):
    data = json.dumps(body).encode()
    delay = 2
    for t in range(tries):
        try:
            req = urllib.request.Request(URL, data=data, method="POST", headers={
                "Accept": "application/json", "Content-Type": "application/json",
                "api-key": API_KEY, "User-Agent": "revision-fetch/1.0"})
            with urllib.request.urlopen(req, timeout=45) as resp:
                return json.loads(resp.read())
        except Exception as e:
            wait = min(delay * (1 + t // 5), 30)
            print(f"  request failed ({e}); retry {t+1}/{tries} in {wait}s", flush=True)
            time.sleep(wait)
    raise SystemExit("page failed after all retries")


total = None
mode = "a" if state["fetched"] else "w"
with open(DEST, mode) as out:
    while True:
        filters = {"assembly_version": "current"}
        if LAST_RELEASE:
            filters["last_release_date"] = LAST_RELEASE
        body = {"taxons": ["482"], "page_size": PAGE, "filters": filters}
        if state["page_token"]:
            body["page_token"] = state["page_token"]
        r = post(body)
        reports = r.get("reports", [])
        total = r.get("total_count", total)
        for rep in reports:
            out.write(json.dumps(rep) + "\n")
        state["fetched"] += len(reports)
        tok = r.get("next_page_token")
        state["page_token"] = tok
        STATE.write_text(json.dumps(state))
        print(f"fetched={state['fetched']}/{total} token={'yes' if tok else 'no'}", flush=True)
        if not tok or not reports:
            break

state["finished_utc"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
STATE.write_text(json.dumps(state, indent=2))
print(f"DONE: {state['fetched']} records, total_count={total}")
