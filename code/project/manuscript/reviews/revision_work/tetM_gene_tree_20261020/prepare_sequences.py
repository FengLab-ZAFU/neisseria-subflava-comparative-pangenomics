from pathlib import Path
import csv, hashlib, json
from urllib.parse import urlencode

ROOT = Path(__file__).resolve().parents[4]
OUT = Path(__file__).resolve().parent
RAW = OUT / 'public_intervals'
RAW.mkdir(exist_ok=True)

def read(p):
    with p.open(encoding='utf-8-sig',newline='') as f:
        return list(csv.DictReader(f,delimiter='\t'))

def fasta(p):
    seqs={};key=None
    with p.open(encoding='utf-8') as f:
        for line in f:
            if line.startswith('>'):
                key=line[1:].split()[0];assert key not in seqs
                seqs[key]=[]
            elif key is not None:
                seqs[key].append(line.strip().upper())
    return {key:''.join(value) for key,value in seqs.items()}

def write(rows,p):
    with p.open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n')
        w.writeheader();w.writerows(rows)

loci=read(ROOT/'results_revision_20261015_main152_species/07_amr_tetm/tetm_locus_evidence.tsv')
manifest={r['accession']:r for r in read(ROOT/'results_revision_20261015_main152_species/00_input/cohort/main152_fasta_manifest.tsv')}
meta={r['strain']:r for r in read(ROOT/'results_revision_20261019_boundary8_analysis/00_cohort/metadata_overall152.tsv')}
localfiles=list((ROOT/'Genomes').rglob('*.fna'))
assembly={};provenance={}
for acc in sorted(set(r['accession'] for r in loci)):
    m=manifest[acc]
    hits=[p for p in localfiles if p.stem in (acc,m['biosample'])]
    good=[p for p in hits if hashlib.sha256(p.read_bytes()).hexdigest()==m['expected_sha256']]
    if good:
        assembly[acc]=fasta(good[0])
        provenance[acc]=str(good[0].relative_to(ROOT))

records=[];downloads=[]
complement=str.maketrans('ACGTRYMKBDHVN','TGCAYRKMVHDBN')
for r in loci:
    acc=r['accession'];lid=r['locus_id'];start=int(r['start']);end=int(r['end'])
    seq=None;source='NCBI_versioned_contig_interval';header=''
    if acc in assembly:
        contig=assembly[acc][r['contig_id']]
        assert len(contig)==int(r['contig_len']),(lid,'contig_length_mismatch')
        seq=contig[start:end]
        if r['strand']=='-':seq=seq.translate(complement)[::-1]
        source='frozen_whole_assembly_SHA256_match'
    strand='2' if r['strand']=='-' else '1'
    url='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?'+urlencode(dict(db='nuccore',id=r['contig_id'],rettype='fasta',retmode='text',seq_start=start+1,seq_stop=end,strand=strand))
    # One local interval is retrieved independently as a sequence-identity check.
    crosscheck=(lid=='GCA_003044385.1_tetM_L1')
    if seq is None or crosscheck:
        p=RAW/(lid+'.fasta')
        if not p.exists():
            downloads.append(dict(locus_id=lid,url=url,path=str(p),contig_id=r['contig_id'],expected_length=end-start,crosscheck=crosscheck))
        else:
            public=fasta(p);assert len(public)==1
            header=next(iter(public));observed=next(iter(public.values()))
            assert header.split(':')[0]==r['contig_id'],(lid,'public_identifier_mismatch',header)
            assert len(observed)==end-start,(lid,'public_length_mismatch')
            if crosscheck:
                assert seq==observed,(lid,'public_local_sequence_mismatch')
            else:seq=observed
    span=end-start;coverage=max(float(r['cov_american_pct']),float(r['cov_dutch_pct']))
    ambiguity=None if seq is None else sum(x not in 'ACGT' for x in seq)
    reasons=[]
    if coverage<95:reasons.append('reference_coverage_below_95_percent')
    if not .95*1935<=span<=1.05*1935:reasons.append('mapped_length_outside_95_to_105_percent')
    if seq is not None and ambiguity/len(seq)>.01:reasons.append('ambiguous_bases_above_1_percent')
    records.append(dict(locus_id=lid,accession=acc,biosample=manifest[acc]['biosample'],
      display_group=meta[acc]['display_group'],mapped_length=span,max_reference_coverage_pct=coverage,
      ambiguous_bases=ambiguity if ambiguity is not None else 'pending',
      sequence_source=source,whole_assembly_sha256_verified=acc in assembly,
      public_crosscheck=crosscheck,sequence_sha256=hashlib.sha256(seq.encode()).hexdigest() if seq is not None else '',
      inclusion='INCLUDED' if not reasons else 'EXCLUDED',exclusion_reason=';'.join(reasons),
      source_contig=r['contig_id'],source_start0=start,source_end0=end,strand=r['strand'],sequence=seq))

(OUT/'download_requests.json').write_text(json.dumps(downloads,indent=2)+'\n',encoding='utf-8')
write([{k:v for k,v in r.items() if k!='sequence'} for r in records],OUT/'sequence_QC_ledger.tsv')
print(json.dumps(dict(mapped_loci=len(records),genomes=len(manifest.keys() & {r['accession'] for r in loci}),
 frozen_local_assemblies=len(assembly),pending_downloads=len(downloads)),indent=2))
if downloads:raise SystemExit(0)
assert all(r['sequence'] is not None for r in records)
included=[r for r in records if r['inclusion']=='INCLUDED']
assert len(included)==52 and len(set(r['accession'] for r in included))==38
for name,rs in [('all57_mapped_loci.fasta',records),('full_length52_loci.fasta',included)]:
    with (OUT/name).open('w',encoding='utf-8',newline='\n') as f:
        for r in rs:f.write('>'+r['locus_id']+'\n'+r['sequence']+'\n')
unique={}
for r in included:unique.setdefault(r['sequence'],[]).append(r)
groups=sorted(unique.items(),key=lambda item:(-len(set(r['accession'] for r in item[1])),-len(item[1]),hashlib.sha256(item[0].encode()).hexdigest()))
mapping=[]
with (OUT/'unique_haplotypes.fasta').open('w',encoding='utf-8',newline='\n') as f:
    for i,(seq,rs) in enumerate(groups,1):
        hid=f'H{i:02d}';f.write('>'+hid+'\n'+seq+'\n')
        for r in rs:mapping.append(dict(haplotype=hid,locus_id=r['locus_id'],accession=r['accession'],biosample=r['biosample'],display_group=r['display_group'],n_loci_in_haplotype=len(rs),n_genomes_in_haplotype=len(set(x['accession'] for x in rs)),sequence_sha256=r['sequence_sha256']))
write(mapping,OUT/'haplotype_locus_genome_mapping.tsv')
print(json.dumps(dict(included_loci=len(included),included_genomes=len(set(r['accession'] for r in included)),unique_haplotypes=len(groups),lengths=sorted(set(len(s) for s in unique))),indent=2))
