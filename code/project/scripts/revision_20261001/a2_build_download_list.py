#!/usr/bin/env python3
"""A2: build FTP download URLs for in-scope accessions and emit a download list.

Usage: python3 a2_build_download_list.py
Output: 01_retrieval_taxonomy/download_urls.tsv (accession<TAB>url)
"""
import csv
import sys
from pathlib import Path

OUT = Path(__file__).resolve().parent.parent.parent / "results_revision_20261001" / "01_retrieval_taxonomy"
scope = OUT / "screening_scope_accessions.tsv"

rows = list(csv.DictReader(open(scope), delimiter="\t"))
n_url = 0
with open(OUT / "download_urls.tsv", "w") as fh:
    for r in rows:
        acc = r["accession"]            # GCA_012345678.1
        asm = r["assembly_name"]        # ASM880699v1
        digits = acc.split("_")[1].split(".")[0]
        g1, g2, g3 = digits[0:3], digits[3:6], digits[6:9]
        prefix = acc.split("_")[0]      # GCA
        base = f"https://ftp.ncbi.nlm.nih.gov/genomes/all/{prefix}/{g1}/{g2}/{g3}/{acc}_{asm}/{acc}_{asm}_genomic.fna.gz"
        fh.write(f"{acc}\t{base}\n")
        n_url += 1
print(f"wrote {n_url} urls -> {OUT/'download_urls.tsv'}")
