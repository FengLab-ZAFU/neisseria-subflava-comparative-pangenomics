#!/usr/bin/env python3
"""M5: annotation coverage at four levels (cluster / CDS copy / genome / protein)
plus run_info and sha256 ledger for 06_function/.
"""
from pathlib import Path
import csv, hashlib, json, os, subprocess, sys
from collections import defaultdict
from datetime import datetime, timezone

ROOT = str(Path(os.environ["NEISSERIA_PROJECT_ROOT"]))
RES = f"{ROOT}/results_revision_20261015_main152_species"
OUT = f"{RES}/06_function"

def main():
    # protein level (unique sequences)
    prot = defaultdict(int)
    for r in csv.DictReader(open(f"{OUT}/protein_hash_bridge.tsv"), delimiter="\t"):
        prot[r["annotation_state"]] += 1
    # copy level + per-genome
    copy = defaultdict(int)
    per_genome = defaultdict(lambda: defaultdict(int))
    for r in csv.DictReader(open(f"{OUT}/copy_cog_map.tsv"), delimiter="\t"):
        copy[r["annotation_state"]] += 1
        per_genome[r["accession"]][r["annotation_state"]] += 1
    # cluster level (representative annotation)
    clus = defaultdict(int)
    for r in csv.DictReader(open(f"{OUT}/cluster_cog_map.tsv"), delimiter="\t"):
        clus[r["annotation_state"]] += 1

    states = ["annotated", "matched_no_cog", "matched_no_emapper_row", "no_sequence_match"]
    with open(f"{OUT}/per_genome_annotation_coverage.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["accession", "n_copies", "n_annotated", "n_matched_no_cog",
                    "n_matched_no_emapper_row", "n_no_sequence_match",
                    "annotated_fraction", "unannotated_fraction"])
        for acc in sorted(per_genome):
            c = per_genome[acc]
            tot = sum(c.values())
            ann = c["annotated"]
            w.writerow([acc, tot, ann, c["matched_no_cog"], c["matched_no_emapper_row"],
                        c["no_sequence_match"], f"{ann/tot:.6f}",
                        f"{(tot-ann)/tot:.6f}"])

    g_ann = [per_genome[a]["annotated"] / sum(per_genome[a].values()) for a in per_genome]
    with open(f"{OUT}/annotation_coverage.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["level", "total", "annotated", "matched_no_cog",
                    "matched_no_emapper_row", "no_sequence_match",
                    "annotated_fraction", "unannotated_fraction", "note"])
        for level, d, note in [
            ("cluster", clus, "cluster-level COG = representative protein annotation (longest protein, ties by accession/gene_id)"),
            ("cds_copy", copy, "one row per PGAP2 occurrence (cluster x accession x gene_id)"),
            ("protein", prot, "unique normalized protein sequences present in the main152 dictionary"),
        ]:
            tot = sum(d.values())
            ann = d["annotated"]
            w.writerow([level, tot, ann, d["matched_no_cog"], d["matched_no_emapper_row"],
                        d["no_sequence_match"], f"{ann/tot:.6f}", f"{(tot-ann)/tot:.6f}", note])
        w.writerow(["genome", len(per_genome), "NA", "NA", "NA", "NA",
                    f"{sum(g_ann)/len(g_ann):.6f}", f"{1-sum(g_ann)/len(g_ann):.6f}",
                    f"per-genome annotated copy fraction: mean reported here; min={min(g_ann):.6f}; max={max(g_ann):.6f}; per-genome rows in per_genome_annotation_coverage.tsv"])

    # ---- run_info
    def ver(cmd):
        try:
            return subprocess.run(cmd, shell=True, capture_output=True, text=True,
                                  timeout=120).stdout.strip() or subprocess.run(
                                      cmd, shell=True, capture_output=True, text=True,
                                      timeout=120).stderr.strip()
        except Exception as e:
            return f"unavailable: {e}"
    run_info = {
        "module": "M5 species-association clusters + COG functional comparison",
        "run_id": "main152_20261015", "seed": 20261015,
        "hostname": os.uname().nodename,
        "finished_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "python": sys.version.split()[0],
        "R": ver("conda run -n r_visual Rscript --version 2>&1 | tail -1"),
        "fisher_method": "R stats::fisher.test two-sided, conditional MLE OR/CI",
        "hash_normalization": "uppercase, whitespace removed, exactly one trailing '*' removed; internal '*' retained",
        "bridge_source": "emapper-3.0.0-beta6 annotations of revision_20261001 PGAP2 proteins (609,982), COG IDs translated via db/cog-24.def.tab; md5 cross-check vs results_revision_20261010_phaseB_v7 hash_bridge_db.tsv fully concordant",
        "supplement": "supplement_emapper/ : emapper-3.0.0-beta6 SIF + eggNOG DB 7.0.0 (same snapshot as historical run) on the 2,854 unique sequences absent from the historical bridge (input preserved as supplement_emapper/input_unmatched_proteins.fasta), 16 CPU, 2026-10-15T04:22:48Z-04:54:38Z; 2,730 queries scanned, 2,666 annotation rows, 188 queries without emapper row -> matched_no_emapper_row; working file 06_function/unmatched_proteins.fasta is empty after the merge",
        "scripts": ["m5_01_build_bridge.py", "m5_02_supp_emapper.sh",
                    "m5_03_association_prep.py", "m5_04_fisher.R",
                    "m5_05_finalize_association.py", "m5_06_cog_enrichment.py",
                    "m5_07_coverage.py"],
    }
    json.dump(run_info, open(f"{OUT}/run_info.json", "w"), indent=1, ensure_ascii=False)

    # ---- sha256 ledger over delivered files (excluding the ledger itself)
    ledger = []
    for dirpath, _dirs, files in os.walk(OUT):
        for fn in sorted(files):
            if fn in ("sha256_ledger.tsv",):
                continue
            p = os.path.join(dirpath, fn)
            rel = os.path.relpath(p, OUT)
            h = hashlib.sha256(open(p, "rb").read()).hexdigest()
            ledger.append((rel, os.path.getsize(p), h))
    with open(f"{OUT}/sha256_ledger.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["relative_path", "bytes", "sha256"])
        w.writerows(ledger)
    sys.stderr.write(f"coverage + ledger done; {len(ledger)} files hashed\n")

if __name__ == "__main__":
    main()
