#!/usr/bin/env python3
"""Recompute cohort flow and pangenome/carriage counts from frozen analytical inputs."""
from pathlib import Path
from collections import Counter
import argparse
import csv
import json


def read(p):
    with p.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f,delimiter='\t'))


def main():
    a=argparse.ArgumentParser();a.add_argument('--project',required=True,type=Path);a.add_argument('--output',required=True,type=Path);x=a.parse_args()
    p=x.project;repo=Path(__file__).resolve().parents[1]
    meta=read(repo/'data/tables/S1_Genome_metadata.tsv');ledger=read(repo/'metadata/screening1303.tsv')
    counts=Counter(r['analysis_group'] for r in meta)
    assert counts=={'subflava':23,'perflava':114,'flavescens':7,'unresolved_species':8}
    ids={r['assembly_accession'] for r in meta};assert len(ids)==152
    dispositions=Counter(r['final_disposition'] for r in ledger)
    assert len(ledger)==1303 and dispositions['retained_overall152']==152
    candidates=len(ledger)-dispositions['ANI_screening_not_passed']
    eligible=candidates-sum(dispositions[k] for k in ['MAG_confirmed','source_role_review','QC_fail','experimental_derivative'])
    representatives=eligible-dispositions['non_representative_duplicate']
    final=representatives-dispositions['non_unique_reference_assignment']
    assert (candidates,eligible,representatives,final)==(327,243,169,152)
    strict=single=core=soft=shell=cloud=total=0
    filename=p/'results_revision_20261015_main152_species/02_pangenome/pgap2_main152/pav_tables/pav_copy.tsv'
    with filename.open(encoding='utf-8') as f:
        rd=csv.reader(f,delimiter='\t');header=next(rd);assert set(header[1:])==ids
        for r in rd:
            z=list(map(int,r[1:]));n=sum(v>0 for v in z);total+=1
            strict+=n==152;single+=all(v==1 for v in z)
            if n/152>=.99:core+=1
            elif n/152>=.95:soft+=1
            elif n/152>=.15:shell+=1
            else:cloud+=1
    assert (total,strict,single,core,soft,shell,cloud)==(13215,1118,962,1347,214,851,10803)
    amr=read(repo/'data/tables/S4_AMR_per_genome.tsv');positive=[r for r in amr if r['tetM_status']=='positive_complete']
    assert len(positive)==38 and {r['accession'] for r in amr}==ids
    carrier=Counter(r['analysis_group'] for r in positive);assert carrier=={'perflava':35,'unresolved_species':3}
    loci=read(repo/'data/tables/S4_tetM_loci.tsv');types=read(repo/'data/tables/S4_tetM_sequence_types.tsv')
    assert len(loci)==57 and len(types)==7 and sum(r['inclusion']=='INCLUDED' for r in loci)==52
    associations=read(repo/'data/tables/S3_Associated_clusters.tsv')
    assert len(associations)==181 and len({r['cluster'] for r in associations})==180
    result={'groups':dict(counts),'screening_dispositions':dict(dispositions),'screened':len(ledger),
        'screening_flow':[len(ledger),candidates,eligible,representatives,final],
        'pangenome':{'total':total,'strict_core':strict,'single_copy_core':single,'core_ge99':core,'soft_core':soft,'shell':shell,'cloud':cloud},
        'tetM':{'carriers':38,'by_group':dict(carrier),'loci':57,'included_sequences':52,'sequence_types':7},
        'associations':{'group_cluster_associations':181,'distinct_clusters':180},
        'cohort_reconstruction_basis':'Frozen screening/unit decisions; this check does not rerun whole-genome SNP calling or rederive clustering from unavailable raw historical alignments.'}
    x.output.parent.mkdir(parents=True,exist_ok=True);x.output.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
