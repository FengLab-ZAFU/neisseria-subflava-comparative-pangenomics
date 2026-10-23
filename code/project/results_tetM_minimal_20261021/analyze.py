from pathlib import Path
import sys,csv,json,hashlib,re
OUT=Path(__file__).resolve().parent;ROOT=OUT.parent
from Bio import SeqIO,Align
def read(p):
 with open(p,encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
def write(rows,name):
 with (OUT/name).open('w',encoding='utf-8',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
refs={x:SeqIO.read(OUT/'public_records'/(x+'.gb'),'genbank') for x in ['GU479464.1','GU479466.1','U09422.1']}
tetrefs={};prot=[];manifest=[]
for acc,rec in refs.items():
 for i,ft in enumerate(rec.features):
  if ft.type!='CDS':continue
  q=ft.qualifiers;desc=' '.join(q.get('gene',[])+q.get('product',[])+q.get('note',[]))
  is_tet=bool(re.search(r'tet\(?m\)?|tetracycline',desc,re.I))
  if is_tet and acc.startswith('GU'):tetrefs[acc]=str(ft.extract(rec.seq)).upper()
  elif not is_tet and 'translation' in q:
   key=acc.replace('.','_')+'__'+str(i);seq=q['translation'][0]
   prot.append('>'+key+'\n'+seq+'\n');manifest.append(dict(marker_id=key,reference=acc,protein_id=';'.join(q.get('protein_id',[])),description=desc,start0=int(ft.location.start),end0=int(ft.location.end),length_aa=len(seq)))
assert len(tetrefs)==2
write(manifest,'backbone_reference_markers.tsv');(OUT/'backbone_markers.faa').write_text(''.join(prot),encoding='utf-8')
a=tetrefs['GU479464.1'];d=tetrefs['GU479466.1'];assert len(a)==len(d)==1935
sites=[i for i,(x,y) in enumerate(zip(a,d)) if x!=y]
hapfile=ROOT/'manuscript/reviews/revision_work/tetM_gene_tree_20261020/unique_haplotypes.fasta'
haps={r.id:str(r.seq).upper() for r in SeqIO.parse(hapfile,'fasta')};types=[]
# Describe reference-diagnostic alleles without inventing a validated typing cutoff.
for h,s in haps.items():
 assert len(s)==1935
 ia=100*sum(x==y for x,y in zip(s,a))/len(s);id=100*sum(x==y for x,y in zip(s,d))/len(s)
 na=sum(s[i]==a[i] for i in sites);nd=sum(s[i]==d[i] for i in sites);other=len(sites)-na-nd
 call='mixed_reference_profile' if na and nd else 'American_reference_profile' if na==len(sites) else 'Dutch_reference_profile' if nd==len(sites) else 'unresolved_reference_profile'
 types.append(dict(haplotype=h,length_nt=len(s),identity_American_pct=ia,identity_Dutch_pct=id,diagnostic_sites=len(sites),American_alleles=na,Dutch_alleles=nd,other_alleles=other,reference_like_assignment=call))
write(types,'haplotype_reference_comparison.tsv');print(json.dumps(types,indent=2))
loci=read(ROOT/'results_revision_20261015_main152_species/07_amr_tetm/tetm_locus_evidence.tsv')
allseq={r.id:str(r.seq).upper() for r in SeqIO.parse(ROOT/'manuscript/reviews/revision_work/tetM_gene_tree_20261020/all57_mapped_loci.fasta','fasta')}
subjects=[];features=[];windows=[]
for r in loci:
 p=OUT/'public_records'/(r['locus_id']+'.gb')
 if not p.exists():continue
 rec=SeqIO.read(p,'genbank');left=max(0,int(r['start'])-30000);start=int(r['start'])-left;end=int(r['end'])-left
 observed=rec.seq[start:end]
 if r['strand']=='-':observed=observed.reverse_complement()
 assert str(observed).upper()==allseq[r['locus_id']],(r['locus_id'],'locus_sequence_mismatch')
 subjects.append('>'+r['locus_id']+'\n'+str(rec.seq)+'\n')
 windows.append(dict(locus_id=r['locus_id'],accession=r['accession'],contig_id=r['contig_id'],window_start0=left,window_end0=left+len(rec),tetM_start_window0=start,tetM_end_window0=end,left_flank_bp=start,right_flank_bp=len(rec)-end,n_CDS=sum(f.type=='CDS' for f in rec.features)))
 for f in rec.features:
  if f.type!='CDS':continue
  q=f.qualifiers;features.append(dict(locus_id=r['locus_id'],start0=left+int(f.location.start),end0=left+int(f.location.end),strand=f.location.strand,gene=';'.join(q.get('gene',[])),locus_tag=';'.join(q.get('locus_tag',[])),protein_id=';'.join(q.get('protein_id',[])),product=';'.join(q.get('product',[]))))
(OUT/'flank_windows.fna').write_text(''.join(subjects),encoding='utf-8');write(windows,'window_verification.tsv')
if features:write(features,'flank_CDS_annotations.tsv')
print('Windows verified',len(windows),'CDS',len(features))

# Use the same reference records for the positive-control searches.
SeqIO.write(list(refs.values()), OUT/'reference_positive_controls.fna', 'fasta')
