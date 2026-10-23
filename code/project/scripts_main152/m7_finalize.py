#!/usr/bin/env python3
"""M7a/M7b finalize (main152, 2026-10-15): run_info JSON + per-file SHA-256
for the module outputs. VFDB results are recorded separately.

Usage: python3 scripts_main152/m7_finalize.py
"""
import csv
import json
import os
import platform
import subprocess
import sys
from datetime import datetime, timezone

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import ROOT, RESULTS, PACKET, sha256_file, RUN_ID, PLAN_ID

OUT = os.path.join(RESULTS, "08_defense_mge_vfdb")

MODULE_OUTPUTS = [
    "tool_coverage.tsv",
    "m7a_evidence_files.tsv",
    "defense_per_genome.tsv",
    "intermediate/defense_type_counts.tsv",
    "per_genome_features.tsv",
    "m7a_coordinate_sources.tsv",
    "m7a_qc_notes.tsv",
    "defense_type_prevalence.tsv",
    "defense_stats.tsv",
    "defense_family_tests.tsv",
    "defense_assembly_model.tsv",
    "model_diagnostics.tsv",
    "sensitivity_defense.tsv",
]

MODULE_SCRIPTS = [
    "scripts_main152/m7a_aggregate.py",
    "scripts_main152/m7b_stats.R",
    "scripts_main152/m7_finalize.py",
]

tool_versions = {}
tv = os.path.join(ROOT, "results_revision_20261010_phaseB_v7/06_annotation_summaries/per_genome_tool_versions.tsv")
for r in csv.DictReader(open(tv), delimiter="\t"):
    tool_versions[r["tool"]] = {"software_version": r["software_version"],
                                "db_version": r["db_version"], "evidence": r["evidence"]}

rver = subprocess.run(["conda", "run", "-n", "r_visual", "Rscript", "--version"],
                      capture_output=True, text=True).stdout.strip()

run_info = {
    "run_id": RUN_ID,
    "plan_id": PLAN_ID,
    "module": "M7a/M7b defense + MGE/GI per-genome summary and group statistics (VFDB recorded separately)",
    "hostname": platform.node(),
    "utc_finished": datetime.now(timezone.utc).isoformat(timespec="seconds"),
    "project_root": ROOT,
    "scripts": {s: sha256_file(os.path.join(ROOT, s)) for s in MODULE_SCRIPTS},
    "software": {
        "python3": platform.python_version(),
        "R": rver,
        "detection_tools": {t: tool_versions.get(t) for t in
                            ["defense_finder", "genomad", "islandpath", "mob_suite",
                             "integron_finder", "isescan"]},
    },
    "inputs": {
        "main152_metadata.tsv": sha256_file(os.path.join(PACKET, "00_input/cohort/main152_metadata.tsv")),
        "main152_accessions.txt": sha256_file(os.path.join(PACKET, "00_input/cohort/main152_accessions.txt")),
        "per_tool_coverage_historical_main152.tsv": sha256_file(os.path.join(PACKET, "00_input/reference_sources/per_tool_coverage_historical_main152.tsv")),
        "assembly_metrics.tsv": sha256_file(os.path.join(RESULTS, "00_preflight/assembly_metrics.tsv")),
        "participant_subset_manifest.tsv": sha256_file(os.path.join(OUT, "participant_subset_manifest.tsv")),
        "per_genome_tool_versions.tsv": sha256_file(tv),
    },
    "randomness": "none (KW, Fisher exact, Wilson, OLS are deterministic; no permutation in this module)",
    "statistical_conventions": {
        "run_status": "success_hits/success_empty = parsed cleanly (hits or zero hits); failed = evidence present but parse error; missing = no evidence file. Failures are never written as 0.",
        "wilson_ci": "z=qnorm(0.975), standard Wilson score interval on the evaluable binomial denominator",
        "kw": "stats::kruskal.test; BH q within exactly the three prespecified burden omnibus tests per analysis_id; epsilon2=max(0,(H-k+1)/(n-k)) with untruncated value also reported",
        "fisher": "stats::fisher.test two-sided, conditional MLE OR/CI; rows=target/rest, cols=present/absent; a=target present etc.; carrier_total>=5 prespecified filter; BH family defense_prevalence (main) and defense_prevalence_participant_subset (sensitivity, corrected separately)",
        "assembly_model": "rank(n_complete_systems, ties.average) ~ species_group + assembly_class + log1p(n_contigs) + genome_size_Mb; coefficients in rank units; partial F of species_group vs assembly-covariate-only model; does NOT adjust for phylogenetic or study-source confounding",
        "group_sizes_main152": {"subflava": 31, "perflava": 114, "flavescens": 7},
        "group_sizes_participant_subset": {"subflava": 24, "perflava": 90, "flavescens": 7},
        "rest_denominators_main152_evaluable": {"subflava": 121, "perflava": 38, "flavescens": 145},
        "flavescens_n7": "no strong negative/specificity claims for the 7-genome group",
    },
    "coordinate_provenance_for_M8": "m7a_coordinate_sources.tsv",
    "outputs": {o: {"sha256": sha256_file(os.path.join(OUT, o)),
                    "n_bytes": os.path.getsize(os.path.join(OUT, o))} for o in MODULE_OUTPUTS},
}

with open(os.path.join(OUT, "run_info_m7ab.json"), "w", encoding="utf-8") as fh:
    json.dump(run_info, fh, indent=2, ensure_ascii=False)
    fh.write("\n")

with open(os.path.join(OUT, "checksums_m7ab.sha256"), "w", encoding="utf-8") as fh:
    for o in MODULE_OUTPUTS + ["run_info_m7ab.json"]:
        fh.write(f"{sha256_file(os.path.join(OUT, o))}  {o}\n")

print("wrote run_info_m7ab.json and checksums_m7ab.sha256 covering",
      len(MODULE_OUTPUTS) + 1, "files")
