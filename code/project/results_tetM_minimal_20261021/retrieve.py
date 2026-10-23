from pathlib import Path
import csv,json,time,urllib.request,urllib.parse,hashlib,sys
from Bio import SeqIO
from io import StringIO
OUT=Path(__file__).resolve().parent
ROOT=OUT.parent
RAW=OUT/'public_records';RAW.mkdir(exist_ok=True)
with (ROOT/'results_revision_20261015_main152_species/07_amr_tetm/tetm_locus_evidence.tsv').open() as f: loci=list(csv.DictReader(f,delimiter='\t'))
jobs=[]
for acc in ['GU479464.1','GU479466.1','U09422.1']:
 jobs.append(dict(key=acc,id=acc,rettype='gbwithparts',retmode='text'))
for r in loci:
 start=max(1,int(r['start'])+1-30000);end=min(int(r['contig_len']),int(r['end'])+30000)
 jobs.append(dict(key=r['locus_id'],id=r['contig_id'],rettype='gbwithparts',retmode='text',seq_start=start,seq_stop=end))
log=[]
for i,j in enumerate(jobs):
 p=RAW/(j['key']+'.gb');args={k:v for k,v in j.items() if k!='key'}
 url='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?'+urllib.parse.urlencode(dict(db='nuccore',**args))
 try:
  if not p.exists():
   req=urllib.request.Request(url,headers={'User-Agent':'academic-revision-sequence-verification/1.0'})
   with urllib.request.urlopen(req,timeout=45) as resp:data=resp.read()
   record=SeqIO.read(StringIO(data.decode()),'genbank')
   p.write_bytes(data);time.sleep(.4)
  record=SeqIO.read(p,'genbank')
  expected=j.get('seq_stop',len(record))-j.get('seq_start',1)+1
  assert len(record)==expected,(len(record),expected)
  log.append(dict(**j,url=url,length=len(record),features=len(record.features),sha256=hashlib.sha256(p.read_bytes()).hexdigest(),status='OK'))
 except Exception as e:log.append(dict(**j,url=url,status='ERROR',error=str(e)))
 (OUT/'retrieval_log.json').write_text(json.dumps(log,indent=2),encoding='utf-8')
 print(i+1,j['key'],log[-1]['status'],flush=True)
