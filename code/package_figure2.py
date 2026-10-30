#!/usr/bin/env python3
"""Check and export the supplied 157-tip Figure 2 iTOL package.

This packages the recorded tree and annotations; it does not infer a new tree.
"""
from pathlib import Path
from collections import defaultdict
import argparse
import csv
import hashlib
import json
import re
import shutil


class Node:
    def __init__(self,children=None,label='',length=0.):
        self.children=children or [];self.label=label;self.length=length


def parse_newick(s):
    tokens=re.findall(r"'(?:[^']|'')*'|\[.*?\]|[(),:;]|[^\s(),:;]+",s)
    i=0
    def go():
        nonlocal i
        ch=[]
        if tokens[i]=='(':
            i+=1;ch.append(go())
            while tokens[i]==',':i+=1;ch.append(go())
            assert tokens[i]==')';i+=1
        label='';length=0.
        if i<len(tokens) and tokens[i] not in {':',',',')',';'}:
            label=tokens[i].strip("'");i+=1
        if i<len(tokens) and tokens[i]==':':i+=1;length=float(tokens[i]);i+=1
        return Node(ch,label,length)
    n=go();assert tokens[i]==';';return n


def tips(n):
    return {n.label} if not n.children else set().union(*(tips(c) for c in n.children))


def walk(n):
    yield n
    for c in n.children:yield from walk(c)


def splitmaps(n):
    all_t=tips(n);lens=defaultdict(float);labs=defaultdict(set)
    for c in walk(n):
        if c is n:continue
        a=tuple(sorted(tips(c)));b=tuple(sorted(all_t-tips(c)))
        k=min((a,b),key=lambda x:(len(x),x));lens[k]+=c.length
        if c.children and c.label:labs[k].add(c.label)
    return lens,labs


def read(p):
    with p.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f,delimiter='\t'))


def main():
    a=argparse.ArgumentParser(description=__doc__)
    a.add_argument('--project',type=Path,required=True);x=a.parse_args()
    src=x.project/'revision2/figure2';out=x.project/'analysis_outputs/figure2_itol'
    tr=parse_newick((src/'01_Figure2_157_outgroup_rooted.nwk').read_text())
    original=parse_newick((x.project/'revision2/phylogeny/main157.treefile').read_text())
    L0,B0=splitmaps(original);L1,B1=splitmaps(tr)
    assert set(L0)==set(L1) and dict(B0)==dict(B1)
    assert max(abs(L0[k]-L1[k]) for k in L0)<1e-12
    meta=read(src/'annotation_metadata.tsv')
    cohort=read(x.project/'revision2/tables/S1_Genome_metadata.tsv')
    cohort_ids={r['assembly_accession'] for r in cohort}
    assert len(meta)==157 and {r['tip'] for r in meta}==tips(tr)
    outgroups=tips(tr)-cohort_ids
    assert len(outgroups)==5 and sorted(len(tips(c)) for c in tr.children)==[5,152]
    assert any(tips(c)==outgroups for c in tr.children)
    assert sum(r['tetM']=='positive' for r in meta)==38
    nodes={frozenset(tips(n)):n for n in walk(tr)}
    perflava=frozenset(r['assembly_accession'] for r in cohort if r['analysis_group']=='perflava')
    assert nodes[perflava].label=='77.6/71'
    for filename in ['02_assignment_class.txt','03_collection_region.txt','04_tetM_carriage.txt','08_strain_name_labels.txt']:
        lines=(src/filename).read_text(encoding='utf-8').splitlines()
        records=[line.split('\t')[0] for line in lines[lines.index('DATA')+1:] if line]
        assert len(records)==157 and set(records)==tips(tr),filename
    out.mkdir(parents=True,exist_ok=True)
    hashes=[]
    for p in sorted(src.iterdir()):
        if p.is_file():
            shutil.copy2(p,out/p.name)
            hashes.append(hashlib.sha256(p.read_bytes()).hexdigest()+'  '+p.name)
    (out/'SHA256SUMS.txt').write_text('\n'.join(hashes)+'\n',encoding='utf-8')
    report={'status':'PASS','tree_tips':157,'study_genomes':152,'outgroups':5,
            'original_unrooted_splits_lengths_support_preserved':True,
            'perflava_support':'77.6/71','tetM_positive_study_genomes':38,
            'output':str(out),'mode':'Package supplied tree and tracks; no ML rerun.'}
    (x.project.parent/'figure2_package_checks.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(report,indent=2))


if __name__=='__main__':main()
