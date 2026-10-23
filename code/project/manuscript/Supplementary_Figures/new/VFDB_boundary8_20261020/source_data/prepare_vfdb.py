import os
from pathlib import Path
import csv, json, math, hashlib, collections

ROOT = Path(os.environ["NEISSERIA_PROJECT_ROOT"])
WORK = Path(__file__).resolve().parent
BASE = ROOT / 'results_revision_20261015_main152_species/08_defense_mge_vfdb'
CURRENT = ROOT / 'results_revision_20261019_boundary8_analysis'
FIG = ROOT / 'manuscript/Supplementary_Figures/new/VFDB_boundary8_20261020'
SRC = FIG / 'source_data'
SRC.mkdir(exist_ok=True)

def read(p):
    with p.open(encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f, delimiter='\t'))

def write(name, rows):
    with (SRC/name).open('w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t',lineterminator='\n')
        w.writeheader(); w.writerows(rows)

def numeric(x):
    if x in ('NA', '', None): return None
    if x in ('Inf', '-Inf'): return x
    try:
        v = float(x)
        return int(v) if v.is_integer() else v
    except (ValueError, AttributeError): return x

meta = {r['strain']:r for r in read(CURRENT/'00_cohort/metadata_overall152.tsv')}
assert len(meta)==152
features = ['fbpA','fbpB','fbpC','hmbR','lbpA','ctrD','lipA','lipB','rfaC','lgtF','rfaK','fHbp','mntC','katA','msrA/B(pilB','recN','mtrC','mtrD','pilT']
display = {f:f for f in features}
display['msrA/B(pilB'] = 'msrA/B (pilB)'
groups = ['subflava','perflava','flavescens','unresolved']
ids = {g:sorted(k for k,v in meta.items() if (v['included_species_comparison144']=='0' if g=='unresolved' else v['display_group']==g)) for g in groups}
assert [len(ids[g]) for g in groups]==[23,114,7,8]
groupof = {k:g for g in groups for k in ids[g]}
matrix = read(BASE/'VFDB_presence_long.tsv')
assert len(matrix)==152*19 and set(r['feature'] for r in matrix)==set(features)
assert len({(r['accession'],r['feature']) for r in matrix})==152*19
det = {(r['accession'],r['feature']):int(r['detected']) for r in matrix}
hits=[]; rawsets=collections.defaultdict(set)
manifest=[]
paths=[CURRENT/'00_cohort/metadata_overall152.tsv', CURRENT/'08_defense_mge_vfdb/vfdb_prevalence.tsv', CURRENT/'08_defense_mge_vfdb/vfdb_feature_tests.tsv',BASE/'VFDB_presence_long.tsv',BASE/'vfdb_candidate_catalog.tsv',BASE/'vfdb/vfdb_run_status.tsv',BASE/'vfdb/vfdb_input_binding.tsv',BASE/'vfdb/vfdb_run_info.txt']
for acc in sorted(meta):
    p=BASE/'vfdb/raw'/f'{acc}.tsv';paths.append(p)
    for r in read(p):
        f=r['GENE']
        assert f in features
        assert float(r['%IDENTITY'])>=80 and float(r['%COVERAGE'])>=80
        rawsets[acc].add(f)
        hits.append(dict(accession=acc,display_group=groupof[acc],feature=f,vfdb_accession=r['ACCESSION'],sequence=r['SEQUENCE'],start=int(r['START']),end=int(r['END']),strand=r['STRAND'],identity_percent=float(r['%IDENTITY']),coverage_percent=float(r['%COVERAGE']),alignment_coverage=r['COVERAGE'],gaps=r['GAPS'],product=r['PRODUCT']))
    assert rawsets[acc]=={f for f in features if det[acc,f]}
for p in paths:
    manifest.append(dict(source=str(p.relative_to(ROOT)).replace('\\','/'),sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
write('input_checksums.tsv',manifest)

tests=read(CURRENT/'08_defense_mge_vfdb/vfdb_feature_tests.tsv')
testmap={(r['feature'],r['target_group']):r for r in tests}
tested=[r for r in tests if r['status']=='TESTED']
assert len(tests)==57 and len(tested)==42
def fisher(a,b,c,d):
    n=a+b+c+d;row=a+b;col=a+c
    lo=max(0,row-(n-col)); hi=min(row,col)
    def prob(x):return math.comb(col,x)*math.comb(n-col,row-x)/math.comb(n,row)
    observed=prob(a)
    return min(1,sum(prob(x) for x in range(lo,hi+1) if prob(x)<=observed*(1+1e-7)))
ps=[]
for r in tests:
    f,g=r['feature'],r['target_group'];target=ids[g];rest=[k for gg in groups[:3] if gg!=g for k in ids[gg]]
    a=sum(det[k,f] for k in target);c=sum(det[k,f] for k in rest)
    counts=[a,len(target)-a,c,len(rest)-c]
    for k,v in zip(['a','b','c','d'],counts):
        if r[k]=='NA': r[k]=str(v)
        else: assert int(r[k])==v
    if r['status']=='TESTED':
        p=fisher(a,len(target)-a,c,len(rest)-c)
        assert abs(p-float(r['p']))<1e-9
        ps.append(p)
order=sorted(range(len(ps)),key=lambda i:ps[i]); q=[0]*len(ps); running=1
for rank in range(len(ps),0,-1):
    i=order[rank-1];running=min(running,ps[i]*len(ps)/rank);q[i]=running
assert all(abs(float(r['q'])-v)<1e-9 for r,v in zip(tested,q))

prevalence=[]
for f in features:
    for g in ['overall']+groups:
        members=list(meta) if g=='overall' else ids[g]
        n=sum(det[k,f] for k in members)
        t=testmap.get((f,g))
        prevalence.append(dict(feature=f,display_feature=display[f],group=g,n_positive=n,n_evaluable=len(members),prevalence=n/len(members),q=numeric(t['q']) if t else None,test_status=t['status'] if t else 'DESCRIPTIVE_ONLY'))
for r in read(CURRENT/'08_defense_mge_vfdb/vfdb_prevalence.tsv'):
    v=next(v for v in prevalence if v['feature']==r['feature'] and v['group']==r['group'])
    assert (v['n_positive'],v['n_evaluable'])==(int(r['n_positive']),int(r['n_evaluable']))
write('VFDB_prevalence_all_groups.tsv',prevalence)
write('VFDB_association_tests_primary144.tsv',tests)
write('VFDB_hit_evidence.tsv',hits)
wide=[]
run={r['accession']:r for r in read(BASE/'vfdb/vfdb_run_status.tsv')}
binding={r['accession']:r for r in read(BASE/'vfdb/vfdb_input_binding.tsv')}
for g in groups:
    for acc in ids[g]:
        r=meta[acc]
        wide.append(dict(accession=acc,display_group=g,included_comparison144=int(r['included_species_comparison144']),assembly_level=r['assembly_level'],n_contigs=int(r['n_contigs']),n_detected_features=sum(det[acc,f] for f in features),**{f:det[acc,f] for f in features}))
write('VFDB_gene_matrix_all152.tsv',wide)
catalog=[]
original={r['feature']:r for r in read(BASE/'vfdb_candidate_catalog.tsv')}
for f in features:
    r=original[f]
    catalog.append(dict(feature=f,display_feature=display[f],vfdb_accession=r['vfdb_accession'],vfdb_product=r['vfdb_product'],n_positive_overall152=sum(det[k,f] for k in meta),min_identity_percent=float(r['min_identity_observed']),min_coverage_percent=float(r['min_coverage_observed'])))
write('VFDB_candidate_catalog.tsv',catalog)
records=[]
for g in groups:
    for acc in ids[g]:
        r=meta[acc];s=run[acc];b=binding[acc]
        assert s['run_status']=='SUCCESS_HITS' and b['expected_sha256']==b['actual_sha256']
        assert int(s['n_vfdb_hits'])==len(rawsets[acc])
        records.append(dict(accession=acc,biosample=r['biosample'],display_group=g,reference_assignment_original=r['reference_ANI_assignment_original'],taxonomy_status=r['taxonomy_status'],included_comparison144=int(r['included_species_comparison144']),assembly_level=r['assembly_level'],n_contigs=int(r['n_contigs']),genome_size_bp=int(r['genome_size']),run_status=s['run_status'],exit_code=int(s['exit_code']),n_detected_features=len(rawsets[acc]),n_sequence_hits=sum(h['accession']==acc for h in hits),assembly_sha256=b['actual_sha256']))
write('VFDB_genome_screening_records.tsv',records)
assert sum(r['n_sequence_hits'] for r in records)==len(hits)

print('VFDB tables verified: 152 genomes, 19 homologs, 42 eligible Fisher tests.')
