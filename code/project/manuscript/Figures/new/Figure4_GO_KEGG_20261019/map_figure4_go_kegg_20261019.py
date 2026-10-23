"""Exact representative-protein sequence mapping of existing GO/KEGG annotations.

No old cluster identifiers are reused. No new sequence annotation is performed.
"""
import os
from pathlib import Path
import csv, hashlib, json, re
from collections import defaultdict, Counter

ROOT = Path(os.environ["NEISSERIA_PROJECT_ROOT"])
BASE = ROOT / 'results_revision_20261015_main152_species'
SRC = ROOT / 'results_revision_20261019_boundary8_analysis'
OUT = ROOT / 'manuscript/Figures/new/Figure4_GO_KEGG_20261019'
OUT.mkdir(parents=True, exist_ok=True)

def table(path):
    with path.open(encoding='utf-8') as f:
        return list(csv.DictReader(f, delimiter='\t'))

def write(name, rows, fields=None):
    if fields is None:
        fields = list(rows[0]) if rows else []
    with (OUT / name).open('w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter='\t', lineterminator='\n')
        w.writeheader(); w.writerows(rows)

def norm(s):
    s = ''.join(s.split()).upper()
    return s[:-1] if s.endswith('*') else s

def parse_annot(path):
    with path.open(encoding='utf-8') as f:
        for line in f:
            if line.startswith('#'): continue
            a = line.rstrip('\n').split('\t')
            if len(a) < 13: continue
            yield a[0], {'GO': set(re.findall(r'GO:\d{7}', a[9])),
                         'KO': set(re.findall(r'K\d{5}', a[11])),
                         'pathway': set(re.findall(r'\d{5}', a[12])),
                         'preferred_name': a[8], 'EC': a[10]}

current = table(BASE / '06_function/cluster_cog_map.tsv')
target_hashes = {r['repre_protein_sha256'] for r in current}
old_ann_path = OUT / 'reference_inputs/pgap2_revision_primary.emapper.annotations'
old_fasta = OUT / 'reference_inputs/total.involved_prot.fa'
assert old_ann_path.exists() and old_fasta.exists(), 'Full historical source pair is required'
# First recover only exact-sequence query-to-hash links needed by current representatives.
# Do not infer matches from historical cluster or gene indices.
query_hash = {}
name = None; chunks = []; n_old_proteins = 0
def map_query(name, seq):
    global n_old_proteins
    n_old_proteins += 1
    digest = hashlib.sha256(norm(''.join(seq)).encode()).hexdigest()
    if digest in target_hashes: query_hash[name] = digest
with old_fasta.open(encoding='utf-8') as f:
    for line in f:
        if line.startswith('>'):
            if name is not None: map_query(name, chunks)
            name=line[1:].split()[0];chunks=[]
        else: chunks.append(line.strip())
    if name is not None: map_query(name,chunks)
old_ann = {q:a for q,a in parse_annot(old_ann_path) if q in query_hash}
by_hash = defaultdict(list)
old_ids_seen = set(old_ann)
for q,a in old_ann.items():
    by_hash[query_hash[q]].append(('historical_revision_exact_sequence',q,a))

supp_path = BASE / '06_function/supplement_emapper/m5_supp.emapper.annotations'
supp_ann = dict(parse_annot(supp_path))
known_hashes=set(query_hash.values()) | set(supp_ann)
supp_input=BASE/'06_function/unmatched_proteins.fasta'
with supp_input.open(encoding='utf-8') as f:
    for line in f:
        if line.startswith('>'): known_hashes.add(line[1:].split()[0])
missing_source_hashes=target_hashes-known_hashes
if missing_source_hashes:
    write('representatives_without_local_full_input.tsv',
          [r for r in current if r['repre_protein_sha256'] in missing_source_hashes])
    print('Representative hashes without downloaded input provenance:',len(missing_source_hashes))
for digest, a in supp_ann.items():
    assert re.fullmatch('[0-9a-f]{64}', digest), digest
    if digest in target_hashes:
        by_hash[digest].append(('main152_supplement_sha256', digest, a))

candidates = table(SRC / '06_function/associated_candidates_ranked.tsv')
ledger = table(SRC / '06_function/cluster_filter_ledger.tsv')
eligible = {r['cluster'] for r in ledger if r['status'] == 'ELIGIBLE'}
candidate_sets = {g: {r['cluster'] for r in candidates if r['target_group'] == g}
                  for g in ['subflava', 'perflava', 'flavescens']}

mapped = []; membership = []; conflicts = []
for r in current:
    digest = r['repre_protein_sha256']; found = by_hash.get(digest, [])
    row = {'cluster': r['cluster'], 'representative_sha256': digest,
           'product': r['product'], 'cog_letters': r['cog_letters'],
           'mapping_status': ('EXACT_MATCH' if found else
                              'EXACT_SEQUENCE_INPUT_WITHOUT_EMAPPER_ROW' if digest in known_hashes else
                              'NO_LOCAL_FULL_ANNOTATION_MATCH'),
           'sources': ','.join(sorted({x[0] for x in found})),
           'source_queries': ','.join(sorted({x[1] for x in found})),
           'eligible_for_association': r['cluster'] in eligible,
           'candidate_groups': ','.join(g for g,s in candidate_sets.items() if r['cluster'] in s)}
    for kind in ['GO','KO','pathway']:
        values = set().union(*(x[2][kind] for x in found)) if found else set()
        row[kind] = ','.join(sorted(values)) or '-'
        for term in sorted(values):
            membership.append({'cluster': r['cluster'], 'annotation_type': kind, 'term_id': term})
        distinct = {tuple(sorted(x[2][kind])) for x in found}
        if len(distinct)>1:
            conflicts.append({'cluster': r['cluster'], 'kind': kind,
                              'distinct_exact_sequence_annotation_sets': len(distinct)})
    row['preferred_name'] = ';'.join(sorted({x[2]['preferred_name'] for x in found if x[2]['preferred_name']!='-'})) or '-'
    row['EC'] = ';'.join(sorted({x[2]['EC'] for x in found if x[2]['EC']!='-'})) or '-'
    mapped.append(row)

write('cluster_GO_KEGG_map.tsv', mapped)
write('cluster_term_membership.tsv', membership)
write('exact_sequence_annotation_conflicts.tsv', conflicts,
      ['cluster','kind','distinct_exact_sequence_annotation_sets'])
coverage = []
all_sets = {'all13215': {r['cluster'] for r in mapped}, 'eligible3165': eligible, **candidate_sets}
for label, ids in all_sets.items():
    z = [r for r in mapped if r['cluster'] in ids]
    coverage.append({'set': label, 'n_clusters':len(z),
                     'n_sequence_provenance_verified':sum(r['representative_sha256'] in known_hashes for r in z),
                     'n_exact_sequence_mapped':sum(r['mapping_status']=='EXACT_MATCH' for r in z),
                     'n_GO':sum(r['GO']!='-' for r in z),
                     'n_KO':sum(r['KO']!='-' for r in z),
                     'n_pathway':sum(r['pathway']!='-' for r in z)})
write('annotation_coverage.tsv', coverage)
summary = {'sources':[str(old_ann_path.relative_to(ROOT)),str(old_fasta.relative_to(ROOT)),str(supp_path.relative_to(ROOT))],
           'historical_proteins_scanned':n_old_proteins,
           'historical_exact_sequence_queries':len(query_hash),
           'old_annot_rows':len(old_ann), 'old_annotation_queries_with_fasta':len(old_ids_seen),
           'old_annotation_queries_missing_fasta':len(set(old_ann)-old_ids_seen),
           'conflicts':len(conflicts), 'coverage':coverage,
           'mapping_rule':'exact normalized representative protein SHA256; union GO/KO/pathway only across identical protein sequences',
           'missing_source_hashes':len(missing_source_hashes),
           'unmatched_status':'absence of an emapper row or term is not inferred functional absence; missing source provenance is separately flagged'}
(OUT/'mapping_summary.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8')
print(json.dumps(summary,indent=2))
