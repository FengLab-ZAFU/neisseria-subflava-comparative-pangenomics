#!/usr/bin/env python3
"""Prepare supplied inputs and run one analysis stage."""
from pathlib import Path
import argparse
import csv
import gzip
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile

REPO = Path(__file__).resolve().parents[1]
B = 'results_revision_20261019_boundary8_analysis/scripts/'


def digest(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''): h.update(b)
    return h.hexdigest()


def prepare(project):
    protected = [REPO / p for p in ['code', 'data', 'docs', 'metadata', '.git']]
    if project == REPO or any(project.is_relative_to(p) or p.is_relative_to(project) for p in protected):
        raise ValueError('Choose a separate workspace, such as work/project')
    project.mkdir(parents=True,exist_ok=True)
    marker = project / '.prepared.json'
    marker.unlink(missing_ok=True)
    with (REPO/'data/frozen_manifest.tsv').open(encoding='utf-8',newline='') as f:
        manifest=list(csv.DictReader(f,delimiter='\t'))
    for r in manifest:
        src=contained(REPO,r['stored_path']);dst=contained(project,r['project_path'])
        if dst.exists() and digest(dst)==r['uncompressed_sha256']: continue
        dst.parent.mkdir(parents=True,exist_ok=True)
        if r['compression'] not in ['none','gzip','gzip_parts']:
            raise ValueError('Unknown compression: '+r['compression'])
        with tempfile.TemporaryDirectory(dir=dst.parent,prefix='.restore-') as tmp:
            tmp=Path(tmp)
            if r['compression']=='gzip_parts':
                joined=tmp/'joined.gz'
                with joined.open('wb') as fout:
                    for part in json.loads(r['parts_json']):
                        with contained(REPO,part).open('rb') as fin: shutil.copyfileobj(fin,fout)
                src=joined
            if digest(src)!=r['stored_sha256']: raise RuntimeError('Frozen input checksum failed: '+r['stored_path'])
            opener=gzip.open if r['compression'] in ['gzip','gzip_parts'] else open
            restored=tmp/'restored'
            with opener(src,'rb') as a,restored.open('wb') as b: shutil.copyfileobj(a,b)
            if digest(restored)!=r['uncompressed_sha256']: raise RuntimeError('Decompressed checksum failed: '+str(dst))
            restored.replace(dst)
    shutil.copytree(REPO/'code/project',project,dirs_exist_ok=True)
    marker.write_text(json.dumps({'manifest_sha256':digest(REPO/'data/frozen_manifest.tsv')}),encoding='utf-8')
    print(f'Prepared {len(manifest)} frozen files in {project}')


def contained(base,relative):
    path=(base/relative).resolve()
    if Path(relative).is_absolute() or not path.is_relative_to(base.resolve()):
        raise ValueError('Path outside workspace: '+relative)
    return path


def jobs(stage):
    b='results_revision_20261019_boundary8_analysis/scripts/'
    p='results_panstripe_minimal_20261020/scripts/'
    return {
      'summary':[('python','@summary')],
      'pangenome':[('python',b+'module_A_pangenome.py')],
      'structure':[('R',b+'module_B_structure.R'),('python',b+'f1_extract_P_matrices.py'),('R',b+'f2_permdisp_vegan.R'),('python',b+'f2b_merge_distance_stats.py')],
      'association':[('python',b+'module_C_association.py')],
      'functional':[('python','manuscript/Figures/new/Figure4_GO_KEGG_20261019/map_figure4_go_kegg_20261019.py'),('python','manuscript/Figures/new/Figure4_GO_KEGG_20261019/prepare_figure4_go_kegg_terms_20261019.py'),('R','manuscript/Figures/new/Figure4_GO_KEGG_20261019/analyze_figure4_go_kegg_20261019.R'),('python','manuscript/Figures/new/Figure4_GO_KEGG_20261019/verify_figure4_go_kegg_20261019.py')],
      'amr':[('python',b+'module_D1_amr.py')],
      'defense':[('python',b+'module_D2_defense.py')],
      'overlap':[('python',b+'module_E_intervals.py'),('python',b+'module_E2_augment_interval_tables.py')],
      'permdisp':[('python',b+'f1_extract_P_matrices.py'),('R',b+'f2_permdisp_vegan.R'),('python',b+'f2b_merge_distance_stats.py')],
      'panstripe':[('R',p+'preflight_and_pilot.R'),('R',p+'run_minimal.R'),('R',p+'diagnostic_plots.R')],
      'tetm-parsimony':[('R','results_tetM_minimal_20261021/parsimony.R')],
      'tetm-sequences':[('python','results_tetM_minimal_20261021/analyze.py'),('blast','flank_windows.fna'),('blast','reference_positive_controls.fna'),('python','results_tetM_minimal_20261021/summarize.py')],
      'tetm-cached-summary':[('python','results_tetM_minimal_20261021/summarize.py')],
      'mge-gi':[('R','manuscript/Figures/new/Figure6_MGE_GI_20261021/analyse_MGE_GI.R')],
      'figure1':[('R','manuscript/Figures/new/Figure1_source/plot_Figure1_Figure2_boundary8_20261019.R')],
      'figure3':[('R','manuscript/Figures/new/Figure3_current_20261021/plot_figure3_current_20261021.R')],
      'figure4':[('R','manuscript/Figures/new/Figure4_GO_KEGG_20261020_v20/plot_figure4_go_kegg_20261020_v20.R')],
      'figure5':[('R','manuscript/Figures/new/Figure5_ABC_20261020_v4/plot_Figure5_ABC.R')],
      'figure6':[('R','manuscript/Figures/new/Figure6_MGE_GI_20261021/plot_Figure6.R')],
      'supplementary1':[('R','manuscript/Supplementary_Figures/Supplementary_Figure_S1_source/plot_ANI_supplement.R')],
      'supplementary2':[('R','manuscript/Figures/new/Figure5_AB_20261020_v2/plot_Figure5_AB.R')],
      'supplementary3':[('python','manuscript/Supplementary_Figures/new/VFDB_boundary8_20261020/source_data/prepare_vfdb.py'),('R','manuscript/Supplementary_Figures/new/VFDB_boundary8_20261020/source_data/plot_vfdb.R')],
    }[stage]


def main():
    a=argparse.ArgumentParser(description=__doc__)
    a.add_argument('command',choices=['prepare','list','run'])
    a.add_argument('stage',nargs='?')
    a.add_argument('--project',type=Path,default=REPO/'work/project')
    a.add_argument('--rscript',default=os.environ.get('RSCRIPT','Rscript'))
    a.add_argument('--tblastn',default=os.environ.get('TBLASTN','tblastn'))
    a.add_argument('--dry-run',action='store_true')
    x=a.parse_args();project=x.project.resolve()
    if x.command=='prepare':
        if x.dry_run:a.error('--dry-run applies to run only')
        prepare(project);return
    stages=['summary','pangenome','structure','association','functional','amr','defense','overlap','permdisp','panstripe','tetm-parsimony','tetm-sequences','tetm-cached-summary','mge-gi','figure1','figure3','figure4','figure5','figure6','supplementary1','supplementary2','supplementary3']
    if x.command=='list': print('\n'.join(stages));return
    if x.stage not in stages: a.error('Choose a stage shown by list')
    marker=project/'.prepared.json'
    if not marker.is_file() or json.loads(marker.read_text())['manifest_sha256']!=digest(REPO/'data/frozen_manifest.tsv'):a.error('Run prepare for this release first')
    if x.stage=='tetm-sequences' and not x.dry_run and not shutil.which(x.tblastn):a.error('tetm-sequences requires BLAST+ tblastn; select it with --tblastn or TBLASTN')
    if x.stage=='overlap' and not x.dry_run:
        folder=project/'results_revision_20261015_main152_species/09_interval_overlap/tables'
        missing=[name for name in ['contigs.tsv','genes.tsv','mge_local.tsv','gi.tsv','whole_contig.tsv'] if not (folder/name).is_file()]
        if missing:a.error('Positional randomization requires the original normalized interval tables: '+', '.join(missing)+'. These inputs are not bundled; see docs/REPRODUCIBILITY.md.')
    env=os.environ.copy();env['NEISSERIA_PROJECT_ROOT']=str(project);env['RSCRIPT']=x.rscript
    env['PYTHONUTF8']='1';env['MPLBACKEND']='Agg';env['PYTHON']=sys.executable
    for language,script in jobs(x.stage):
        if language=='blast':
            folder=project/'results_tetM_minimal_20261021'
            output='backbone_tblastn.tsv' if script=='flank_windows.fna' else 'positive_control_tblastn.tsv'
            command=[x.tblastn,'-query',str(folder/'backbone_markers.faa'),'-subject',str(folder/script),'-evalue','1e-5','-seg','yes','-outfmt','6 qseqid sseqid pident length qstart qend qlen sstart send evalue bitscore','-out',str(folder/output)]
            print('Running:',subprocess.list2cmdline(command),flush=True)
            if not x.dry_run:subprocess.run(command,cwd=project,env=env,check=True)
            continue
        entry=REPO/'code/recompute_summary.py' if script=='@summary' else project/script
        command=[sys.executable if language=='python' else x.rscript,str(entry)]
        if script=='@summary':command+=['--project',str(project),'--output',str(project.parent/'recomputed_summary.json')]
        if x.stage=='figure1':command += [str(project/'manuscript/Figures/new'),'--fig1-only']
        if not entry.exists():raise FileNotFoundError(entry)
        print('Running:', ' '.join(command),flush=True)
        if not x.dry_run:
            cwd=project/'results_revision_20261019_boundary8_analysis' if script.startswith(B) else project
            subprocess.run(command,cwd=cwd,env=env,check=True)
    if x.stage in ['permdisp','structure'] and not x.dry_run:
        corrected=project/'results_revision_20261019_boundary8_permdisp_fix/05_structure'
        for name in ['distance_statistics.tsv','distance_stats_run_info.txt']:
            shutil.copy2(corrected/name,project/'results_revision_20261019_boundary8_analysis/05_structure'/name)


if __name__=='__main__':
    main()
