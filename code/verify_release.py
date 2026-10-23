#!/usr/bin/env python3
"""Check release hashes, source syntax, cohort membership and large-file limits offline."""
from pathlib import Path
import argparse
import ast
import csv
import hashlib
import json
import sys

ROOT=Path(__file__).resolve().parents[1]


def digest(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()


def read(path):
    with path.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f,delimiter='\t'))


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--no-hashes',action='store_true');args=parser.parse_args()
    errors=[];n=0;listed=set()
    for line in (ROOT/'docs/SHA256SUMS.txt').read_text(encoding='utf-8').splitlines():
        expected,relative=line.split('  ',1);p=(ROOT/relative).resolve()
        if not p.is_relative_to(ROOT) or relative in listed:
            errors.append('Invalid or duplicate checksum path: '+relative)
            continue
        listed.add(relative)
        if not args.no_hashes:
            n+=1
            if not p.is_file() or digest(p)!=expected:errors.append('Checksum: '+relative)
    excluded={'.git','work','outputs','__pycache__','.Rproj.user'}
    actual={p.relative_to(ROOT).as_posix() for p in ROOT.rglob('*') if p.is_file()
            and not excluded.intersection(p.relative_to(ROOT).parts) and p.suffix!='.pyc'}
    actual.discard('docs/SHA256SUMS.txt')
    errors.extend('Unlisted release file: '+p for p in sorted(actual-listed))
    errors.extend('Missing release file: '+p for p in sorted(listed-actual))
    py=0
    for p in (ROOT/'code').rglob('*.py'):
        py+=1
        try:ast.parse(p.read_text(encoding='utf-8-sig'),filename=str(p))
        except SyntaxError as exc:errors.append(str(exc))
    meta=read(ROOT/'metadata/overall152.tsv');comp=read(ROOT/'metadata/comparison144.tsv')
    ids={r['strain'] for r in meta};cids={r['strain'] for r in comp}
    uids=set((ROOT/'metadata/unresolved8.txt').read_text().split())
    if len(meta)!=len(ids) or len(ids)!=152 or len(comp)!=len(cids) or len(cids)!=144 or len(uids)!=8 or cids&uids or cids|uids!=ids:errors.append('Cohort membership')
    screening=read(ROOT/'metadata/screening1303.tsv')
    if len(screening)!=1303 or len({r['assembly_accession'] for r in screening})!=1303:errors.append('Screening ledger')
    provenance=read(ROOT/'docs/SCRIPT_PROVENANCE.tsv')
    native={p.relative_to(ROOT).as_posix() for p in (ROOT/'code/project').rglob('*') if p.is_file() and p.suffix in ['.py','.R','.sh']}
    recorded=[r['release_path'] for r in provenance]
    if set(recorded)!=native or len(recorded)!=len(set(recorded)):errors.append('Script provenance inventory')
    if not args.no_hashes:
        for row in provenance:
            if row['release_path'] in native and digest(ROOT/row['release_path'])!=row['release_sha256']:
                errors.append('Script provenance checksum: '+row['release_path'])
    frozen=read(ROOT/'data/frozen_manifest.tsv')
    destinations=set()
    for row in frozen:
        relative=row['project_path'];destination=Path(relative)
        if destination.is_absolute() or '..' in destination.parts or relative in destinations:
            errors.append('Invalid or duplicate frozen destination: '+relative)
        destinations.add(relative)
        parts=json.loads(row['parts_json']) if row['compression']=='gzip_parts' else [row['stored_path']]
        if row['compression'] not in ['none','gzip','gzip_parts']:errors.append('Unknown compression: '+relative)
        for part in parts:
            candidate=(ROOT/part).resolve()
            if not candidate.is_relative_to(ROOT/'data/frozen') or not candidate.is_file():
                errors.append('Missing or invalid frozen input: '+part)
    large=[str(p.relative_to(ROOT)) for p in ROOT.rglob('*') if p.is_file() and not excluded.intersection(p.relative_to(ROOT).parts) and p.stat().st_size>=100*1024**2]
    errors.extend('GitHub file-size limit: '+p for p in large)
    report={'status':'FAIL' if errors else 'PASS','hashed_files':n,'frozen_files':len(frozen),'python_sources_parsed':py,'cohort':[152,144,8],'errors':errors}
    print(json.dumps(report,indent=2));return bool(errors)


if __name__=='__main__':sys.exit(main())
