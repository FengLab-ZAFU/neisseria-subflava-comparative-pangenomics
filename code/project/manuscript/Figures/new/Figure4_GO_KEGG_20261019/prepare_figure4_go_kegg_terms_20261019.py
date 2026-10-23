import os
"""Prepare GO ancestry and KEGG pathway sets from exact-sequence annotations."""
from pathlib import Path
import csv,json,re
from functools import lru_cache
from collections import defaultdict

ROOT=Path(os.environ["NEISSERIA_PROJECT_ROOT"])
OUT=ROOT/'manuscript/Figures/new/Figure4_GO_KEGG_20261019'
def read(p):
    with p.open(encoding='utf-8') as f:return list(csv.DictReader(f,delimiter='\t'))
def write(name,rows,fields):
    with (OUT/name).open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)

ontology={};alts={};replacements={};term=None;section=None;version=None
def finish(t):
    if not t or 'id' not in t:return
    if t.get('obsolete'):
        if len(t.get('replaced_by',[]))==1:replacements[t['id']]=t['replaced_by'][0]
        return
    ontology[t['id']]=t
    for a in t.get('alt',[]):alts[a]=t['id']
obo=OUT/'reference_inputs/go-basic.obo'
with obo.open(encoding='utf-8') as f:
    for line in f:
        line=line.rstrip()
        if line.startswith('data-version:'):version=line.split(': ',1)[1]
        if line=='[Term]':
            finish(term);term={'parents':[],'alt':[],'replaced_by':[]};section='term';continue
        if line.startswith('['):finish(term);term=None;section='other';continue
        if section!='term' or term is None:continue
        if line.startswith('id: '):term['id']=line[4:]
        elif line.startswith('name: '):term['name']=line[6:]
        elif line.startswith('namespace: '):term['namespace']=line[11:]
        elif line.startswith('alt_id: '):term['alt'].append(line[8:])
        elif line.startswith('is_a: '):term['parents'].append(line[6:].split()[0])
        elif line.startswith('relationship: part_of '):term['parents'].append(line.split()[2])
        elif line=='is_obsolete: true':term['obsolete']=True
        elif line.startswith('replaced_by: '):term['replaced_by'].append(line[13:])
finish(term)
assert len(ontology)>30000
def canonical(x):
    return alts.get(x,replacements.get(x,x))
@lru_cache(None)
def ancestry(x):
    x=canonical(x)
    if x not in ontology:return frozenset()
    p=ontology[x]['parents']
    return frozenset({x}).union(*(ancestry(canonical(y)) for y in p))

path_names={}
with (OUT/'reference_inputs/kegg_pathway_names.tsv').open(encoding='utf-8') as f:
    for l in f:
        a=l.rstrip().split('\t',1)
        if len(a)==2:path_names[re.search(r'\d{5}',a[0]).group()]=a[1]
rows=read(OUT/'cluster_GO_KEGG_map.tsv')
membership=[];metadata={};missing=defaultdict(set);candidate_sets=defaultdict(set)
roots={'GO:0008150','GO:0003674','GO:0005575'}
for r in rows:
    direct=set() if r['GO']=='-' else set(r['GO'].split(','))
    go=set()
    for g in direct:
        c=canonical(g)
        if c not in ontology:missing[g].add(r['cluster']);continue
        go.update(ancestry(c))
    go-=roots
    paths=set() if r['pathway']=='-' else set(r['pathway'].split(','))
    ko=set() if r['KO']=='-' else set(r['KO'].split(','))
    for kind,terms in [('GO',go),('KEGG_pathway',paths),('KO',ko)]:
        for t in sorted(terms):
            membership.append({'cluster':r['cluster'],'annotation_type':kind,'term_id':t})
            if kind=='GO':
                metadata[(kind,t)]={'annotation_type':kind,'term_id':t,'description':ontology[t]['name'],
                                    'namespace':ontology[t]['namespace']}
            elif kind=='KEGG_pathway':
                metadata[(kind,t)]={'annotation_type':kind,'term_id':t,'description':path_names.get(t,'UNRESOLVED_PATHWAY_NAME'),
                                    'namespace':'KEGG_pathway'}
            else:
                metadata[(kind,t)]={'annotation_type':kind,'term_id':t,'description':t,'namespace':'KEGG_orthology'}
    for g in filter(None,r['candidate_groups'].split(',')):candidate_sets[g].add(r['cluster'])
write('cluster_term_membership_expanded.tsv',membership,['cluster','annotation_type','term_id'])
write('term_metadata.tsv',list(metadata.values()),['annotation_type','term_id','description','namespace'])
write('GO_unresolved_terms.tsv',[{'term_id':g,'n_clusters':len(x)} for g,x in missing.items()],['term_id','n_clusters'])
paired=[]
for i,g in enumerate(['subflava','perflava','flavescens']):
    for h in ['subflava','perflava','flavescens'][i+1:]:
        for k in ['KO','pathway']:
            z={x['cluster']:set() if x[k]=='-' else set(x[k].split(',')) for x in rows}
            a=set().union(*(z[c] for c in candidate_sets[g]));b=set().union(*(z[c] for c in candidate_sets[h]))
            paired.append({'group1':g,'group2':h,'annotation_type':k,'n_terms_group1':len(a),'n_terms_group2':len(b),
                           'n_shared_terms':len(a&b),'shared_terms':','.join(sorted(a&b))})
write('candidate_function_overlap.tsv',paired,['group1','group2','annotation_type','n_terms_group1','n_terms_group2','n_shared_terms','shared_terms'])
summary={'GO_ontology_version':version,'GO_relations':['is_a','part_of'],
         'GO_roots_excluded':sorted(roots),'active_GO_terms':len(ontology),
         'unresolved_source_GO_terms':len(missing),'membership_rows':len(membership),
         'KEGG_pathway_names':len(path_names),'annotation_mapping':'exact representative sequence SHA256',
         'scope':'gene-cluster-level annotation; no pathway completeness or activity inference'}
(OUT/'ontology_preparation_summary.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8')
print(json.dumps(summary,indent=2))
