import os
"""Independent finite-support Fisher/BH checks and reader-facing annotations."""
from pathlib import Path
from collections import defaultdict
import csv,hashlib,json,math,shutil

ROOT=Path(os.environ["NEISSERIA_PROJECT_ROOT"])
OUT=ROOT/'manuscript/Figures/new/Figure4_GO_KEGG_20261019'
SRC=ROOT/'results_revision_20261019_boundary8_analysis'
def read(p):
 with p.open(encoding='utf-8') as f:return list(csv.DictReader(f,delimiter='\t'))
def write(p,r,fields=None):
 fields=fields or list(r[0])
 with p.open('w',encoding='utf-8',newline='') as f:
  w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(r)
def bh(p):
 order=sorted(range(len(p)),key=p.__getitem__);q=[0.]*len(p);prev=1.
 for j in range(len(p)-1,-1,-1):
  i=order[j];prev=min(prev,p[i]*len(p)/(j+1));q[i]=prev
 return q
def fisher(a,b,c,d):
 N=a+b;K=a+c;M=a+b+c+d
 support=range(max(0,N-(M-K)),min(N,K)+1)
 den=math.comb(M,N)
 prob={x:math.comb(K,x)*math.comb(M-K,N-x)/den for x in support}
 observed=prob[a]
 return min(1.,sum(v for v in prob.values() if v<=observed*(1+1e-7)))
mapping=read(OUT/'cluster_GO_KEGG_map.tsv')
eligible={r['cluster'] for r in mapping if r['eligible_for_association']=='True'}
cand=read(SRC/'06_function/associated_candidates_ranked.tsv')
assert len(eligible)==3165 and len(cand)==181
candidate_sets={g:{r['cluster'] for r in cand if r['target_group']==g} for g in ['subflava','perflava','flavescens']}
members=defaultdict(set)
for r in read(OUT/'cluster_term_membership_expanded.tsv'):
 if r['cluster'] in eligible:members[(r['annotation_type'],r['term_id'])].add(r['cluster'])
audit=[]
for kind in ['GO','KEGG_pathway']:
 tests=read(OUT/f'{kind}_all_tests.tsv')
 universe=set().union(*(v for (k,t),v in members.items() if k==kind))
 pvals=[]
 for r in tests:
  fg=candidate_sets[r['target_group']]&universe;bg=universe-fg
  z=members[(kind,r['term_id'])]
  a=len(fg&z);c=len(bg&z);b=len(fg)-a;d=len(bg)-c
  assert (a,len(fg),c,len(bg))==tuple(int(r[k]) for k in ['foreground_count','foreground_total','background_count','background_total'])
  pvals.append(fisher(a,b,c,d))
 qvals=bh(pvals)
 maxp=max(abs(p-float(r['p'])) for p,r in zip(pvals,tests))
 maxq=max(abs(q-float(r['q'])) for q,r in zip(qvals,tests))
 assert maxp<1e-10 and maxq<1e-10
 audit.append({'family':kind,'n_tests':len(tests),'annotated_universe':len(universe),
               'n_q_lt_005':sum(q<.05 for q in qvals),'max_P_difference':maxp,'max_q_difference':maxq,'status':'PASS'})
write(OUT/'independent_GO_KEGG_test_audit.tsv',audit)
meta={(r['annotation_type'],r['term_id']):r['description'] for r in read(OUT/'term_metadata.tsv')}
bycl={r['cluster']:r for r in mapping}
assoc={(r['feature'],r['target_group']):r for r in read(SRC/'06_function/association_all_tests.tsv')}
rows=[]
for c in cand:
 r=bycl[c['cluster']];a=assoc[(c['cluster'],c['target_group'])]
 direct_go=[] if r['GO']=='-' else r['GO'].split(',')
 paths=[] if r['pathway']=='-' else r['pathway'].split(',')
 rows.append({'target_group':c['target_group'],'cluster':c['cluster'],
  'target_carriers':a['a'],'target_total':a['n_target_evaluable'],
  'rest_carriers':a['c'],'rest_total':a['n_rest_evaluable'],'gene_association_q':a['q'],
  'product':r['product'],'preferred_name':r['preferred_name'],'COG':r['cog_letters'],
  'GO_IDs':r['GO'],'GO_descriptions':'; '.join(meta.get(('GO',t),t) for t in direct_go) or '-',
  'KO_IDs':r['KO'],'KEGG_pathway_IDs':r['pathway'],
  'KEGG_pathway_names':'; '.join(meta.get(('KEGG_pathway',t),t) for t in paths) or '-',
  'EC':r['EC'],'mapping_status':r['mapping_status']})
write(OUT/'candidate181_GO_KEGG_annotations.tsv',rows)
sig=[r for r in read(OUT/'KEGG_pathway_all_tests.tsv') if float(r['q'])<.05]
write(OUT/'KEGG_significant_results.tsv',sig)
refs=[OUT/'reference_inputs/pgap2_revision_primary.emapper.annotations',OUT/'reference_inputs/total.involved_prot.fa',
      OUT/'reference_inputs/go-basic.obo',OUT/'reference_inputs/kegg_pathway_names.tsv',
      ROOT/'results_revision_20261015_main152_species/06_function/supplement_emapper/m5_supp.emapper.annotations']
expected={'pgap2_revision_primary.emapper.annotations':'1f28d682fb719e5b8f1510a8fde52e64cbc96897a1bdc6b02788eba58a05e236',
          'total.involved_prot.fa':'025d2095c4f0072ca898be16ff155e110143a4dcfeae8f7ff10b02a6e0d8ab16'}
manifest=[]
for p in refs:
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 digest=h.hexdigest()
 if p.name in expected:assert digest==expected[p.name]
 manifest.append({'source':str(p.relative_to(ROOT)),'bytes':p.stat().st_size,'sha256':digest,
                  'remote_hash_verified':p.name in expected})
write(OUT/'source_manifest.tsv',manifest)
print(json.dumps(audit,indent=2))
