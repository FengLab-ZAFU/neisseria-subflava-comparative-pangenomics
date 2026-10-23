from pathlib import Path
import csv
import hashlib

OUT = Path(__file__).resolve().parent
ROOT = OUT.parent


def read(path):
    with path.open(encoding='utf-8-sig', newline='') as stream:
        return list(csv.DictReader(stream, delimiter='\t'))


def write(rows, name):
    with (OUT/name).open('w', encoding='utf-8', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)

from collections import defaultdict,Counter
markers={r['marker_id']:r for r in read(OUT/'backbone_reference_markers.tsv')}
windows={r['locus_id']:r for r in read(OUT/'window_verification.tsv')}
hits=defaultdict(list);qualified=[]
def coverage(intervals):
 points=set()
 for a,b in intervals:points.update(range(a,b+1))
 return len(points)
def merged_hits(path):
 groups=defaultdict(list)
 with path.open() as f:
  for r in csv.reader(f,delimiter='\t'):
   if float(r[2])>=30:groups[(r[0],r[1],int(r[7])<int(r[8]))].append(r)
 for (qid,sid,strand),rs in groups.items():
  rs.sort(key=lambda r:min(int(r[7]),int(r[8])))
  chunks=[]
  for r in rs:
   if not chunks or min(int(r[7]),int(r[8]))-max(max(int(x[7]),int(x[8])) for x in chunks[-1])>3*int(r[6])+200:chunks.append([])
   chunks[-1].append(r)
  for chunk in chunks:
   cov=coverage([(int(r[4]),int(r[5])) for r in chunk])/int(chunk[0][6])
   if cov<.7:continue
   yield qid,sid,cov,chunk
control=set(q for q,s,c,rs in merged_hits(OUT/'positive_control_tblastn.tsv') if s==markers[q]['reference'])
assert len(control)==75 and set(markers)-control=={'GU479464_1__2'}
write([dict(marker_id=k,description=v['description'],length_aa=int(v['length_aa']),self_match_detected=k in control) for k,v in markers.items()],'positive_control_summary.tsv')
for qid,lid,cov,rs in merged_hits(OUT/'backbone_tblastn.tsv'):
 w=windows[lid];start=min(min(int(r[7]),int(r[8])) for r in rs)-1;end=max(max(int(r[7]),int(r[8])) for r in rs)
 if start<int(w['tetM_end_window0']) and end>int(w['tetM_start_window0']):continue
 m=markers[qid];q=dict(locus_id=lid,marker_id=qid,reference=m['reference'],description=m['description'],identity_pct=sum(float(r[2])*int(r[3]) for r in rs)/sum(int(r[3]) for r in rs),query_coverage_pct=cov*100,aligned_aa=sum(int(r[3]) for r in rs),start0=int(w['window_start0'])+start,end0=int(w['window_start0'])+end,evalue=min(float(r[9]) for r in rs))
 hits[lid].append(q);qualified.append(q)
if qualified:write(qualified,'qualified_backbone_hits.tsv')
flanks=[]
for lid,w in windows.items():
 hs=hits[lid];names=sorted(set(h['description'] for h in hs));tn=sum(h['reference']=='U09422.1' for h in hs);pl=sum(h['reference']=='GU479466.1' for h in hs)
 flanks.append(dict(locus_id=lid,accession=w['accession'],left_flank_bp=int(w['left_flank_bp']),right_flank_bp=int(w['right_flank_bp']),CDS_annotations=int(w['n_CDS']),Tn916_non_tetM_marker_hits=tn,plasmid_reference_marker_hits=pl,matched_markers='; '.join(names) if names else 'none',backbone_interpretation='short_ORF_homology_only' if tn else 'no_qualifying_reference_backbone_match'))
write(flanks,'flank_summary.tsv')
sources=[ROOT/'results_revision_20261015_main152_species/03_phylogeny/strict_singlecopy_core/main152_ML_unrooted.nwk',ROOT/'results_revision_20261015_main152_species/07_amr_tetm/amr_tetm_per_genome.tsv',ROOT/'results_revision_20261015_main152_species/07_amr_tetm/tetm_locus_evidence.tsv']
sources+=list((ROOT/'manuscript/reviews/revision_work/tetM_gene_tree_20261020').glob('*.fasta'))
write([dict(path=str(p.relative_to(ROOT)),sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in sources],'input_sha256.tsv')
print('flanks',len(flanks),'genomes',len(set(r['accession'] for r in flanks)),'full_30kb_each',sum(min(r['left_flank_bp'],r['right_flank_bp'])==30000 for r in flanks),'short_flank_under2kb',sum(min(r['left_flank_bp'],r['right_flank_bp'])<2000 for r in flanks))
print('annotation_evaluable',sum(r['CDS_annotations']>0 for r in flanks),'no_annotation',sum(r['CDS_annotations']==0 for r in flanks))
