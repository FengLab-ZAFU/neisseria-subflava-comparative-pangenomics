#!/usr/bin/env python3
"""Prepare PERMDISP index matrices from saved matrices or an expanded draw table.

Each matrix is indexed by permutation and accession. The R stage independently
checks the supplied matrices against their recorded random seeds.
"""
from pathlib import Path
import csv, hashlib, os, sys
import shutil
from collections import defaultdict

ROOT = str(Path(os.environ["NEISSERIA_PROJECT_ROOT"]))
S = f"{ROOT}/results_revision_20261019_boundary8_analysis"
F = f"{ROOT}/results_revision_20261019_boundary8_permdisp_fix"
OUTD = f"{F}/05_structure/P_matrices"
DRAW = f"{S}/05_structure/permutation_draws.tsv"


def use_saved_matrices():
    """Restore the saved index matrices when the expanded draw table is absent."""
    rows=[]
    for scope in ['primary144','sensitivity152']:
        for design in ['free','restricted_within_study_namespace']:
            name=f'P_{scope}_{design}.tsv'
            source=Path(S)/'05_structure/P_matrices'/name
            with source.open(newline='',encoding='utf-8') as stream:
                reader=csv.reader(stream,delimiter='\t');header=next(reader);matrix=list(reader)
            n=len(header)-1
            assert n== (144 if scope=='primary144' else 152)
            assert len(set(header[1:]))==n and len(matrix)==999
            for index,row in enumerate(matrix,1):
                assert int(row[0])==index and sorted(map(int,row[1:]))==list(range(1,n+1))
            target=Path(OUTD)/name
            shutil.copy2(source,target)
            rows.append(dict(scope=scope,design=design,n_accessions=n,n_permutations=len(matrix),
                file='05_structure/P_matrices/'+name,sha256=hashlib.sha256(target.read_bytes()).hexdigest()))
    with open(f'{F}/05_structure/P_matrix_manifest.tsv','w',encoding='utf-8',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n')
        writer.writeheader();writer.writerows(rows)
    print('Restored and checked four saved permutation matrices')

def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()

os.makedirs(OUTD, exist_ok=True)
if not os.path.isfile(DRAW):
    use_saved_matrices()
    sys.exit(0)
# draws[(scope, design)][perm_idx] = {accession: position}
draws = defaultdict(dict)
acc_order = defaultdict(list)
with open(DRAW) as fh:
    rdr = csv.DictReader(fh, delimiter="\t")
    for r in rdr:
        if r["test_use"] != "PERMDISP":
            continue
        key = (r["scope"], r["design"])
        pi = int(r["permutation_index"])
        d = draws[key].setdefault(pi, {})
        d[r["accession"]] = int(r["permuted_position"])
        if pi == 1:
            acc_order[key].append(r["accession"])

manifest = []
for (scope, design), perms in sorted(draws.items()):
    accs = acc_order[(scope, design)]
    n = len(accs)
    assert len(set(accs)) == n, f"{scope}/{design}: duplicate accessions in draw 1"
    for pi, d in perms.items():
        assert set(d) == set(accs), f"{scope}/{design} draw {pi}: accession set mismatch"
    path = f"{OUTD}/P_{scope}_{design}.tsv"
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["permutation_index"] + accs)
        for pi in sorted(perms):
            d = perms[pi]
            w.writerow([pi] + [d[a] for a in accs])
    manifest.append({"scope": scope, "design": design, "n_accessions": n,
                     "n_permutations": len(perms),
                     "file": f"05_structure/P_matrices/P_{scope}_{design}.tsv",
                     "sha256": sha256(path)})
    print(f"{scope}/{design}: n={n} B={len(perms)}")

with open(f"{F}/05_structure/P_matrix_manifest.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=["scope", "design", "n_accessions", "n_permutations", "file", "sha256"],
                       delimiter="\t", lineterminator="\n")
    w.writeheader()
    w.writerows(manifest)
print("P matrix manifest written.")
