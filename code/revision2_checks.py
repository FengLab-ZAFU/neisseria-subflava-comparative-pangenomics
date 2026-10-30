#!/usr/bin/env python3
"""Reconcile retained revision results against frozen directional ANI and locus data.

This checks supplied results. It does not repeat fastANI, genome prediction,
recombination analysis, or epidemiological dereplication.
"""
from pathlib import Path
from collections import Counter, defaultdict
import argparse
import csv
import hashlib
import json
import math


def read(path):
    with path.open(encoding='utf-8-sig', newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--project', required=True, type=Path)
    args = parser.parse_args()
    project = args.project.resolve()
    repo = Path(__file__).resolve().parents[1]
    tables = project / 'revision2/tables'
    for source in sorted((repo / 'data/tables').glob('S[1-6]_*.tsv')):
        assert hashlib.sha256(source.read_bytes()).digest() == hashlib.sha256((tables / source.name).read_bytes()).digest(), source.name
    assert not list(tables.glob('S2_Panstripe_*.tsv'))
    assert not any('WITHDRAWN' in str(r) or 'Panstripe' in str(r) for r in read(tables / 'S2_Notes.tsv'))
    meta = read(tables / 'S1_Genome_metadata.tsv')
    groups = {r['assembly_accession']: r['analysis_group'] for r in meta}
    assert Counter(groups.values()) == {'subflava': 23, 'perflava': 114, 'flavescens': 7, 'unresolved_species': 8}
    raw = read(project / 'revision2/ani/raw_scores_long.tsv')
    expected = read(project / 'revision2/ani/membership_sensitivity.tsv')
    channels = defaultdict(dict)
    for row in raw:
        key = row['accession'], row['channel']
        assert row['direction'] not in channels[key], key
        channels[key][row['direction']] = row
    assert len(raw) == 912 and len(channels) == 456
    assert len(expected) == 152 and {r['accession'] for r in expected} == set(groups)
    perflava = []
    for row in expected:
        passed = []
        for channel, prefix in [('subflava','subflava'), ('perflava_single','perflava_single'), ('flavescens','flavescens')]:
            pair = channels[row['accession'], channel]
            assert len(pair) == 2
            mean = sum(float(r['ANI']) for r in pair.values()) / 2
            af = min(int(r['matched']) / int(r['total']) for r in pair.values())
            assert math.isclose(mean, float(row[prefix+'_mean_ANI']), abs_tol=1e-6)
            assert math.isclose(af, float(row[prefix+'_min_AF']), abs_tol=1e-6)
            passed.append(mean >= 95 and af >= .70)
            if groups[row['accession']] == 'perflava' and channel == 'perflava_single':
                assert passed[-1]
                perflava.append(mean)
        status = 'PASS' if any(passed) else 'FAIL'
        assert row['single_ref_membership'] == status == 'PASS'
    loci = read(tables / 'S4_tetM_loci.tsv')
    assert len(loci) == 57 and len({r['locus_id'] for r in loci}) == 57
    assert len({r['accession'] for r in loci}) == 38
    assert sum(r['inclusion'] == 'INCLUDED' for r in loci) == 52
    assert Counter(r['location_grade'] for r in loci) == {
        'chromosomal_verified': 3, 'assembly_supported_association': 1,
        'undetermined_small_contig': 39, 'undetermined_contig_boundary': 14}
    order = read(project / 'revision2/supplementary2/carrier_order.tsv')
    assert len(order) == 38 and {r['accession'] for r in order} == {r['accession'] for r in loci}
    assert [int(r['column']) for r in order] == list(range(1,39))
    totals = dict.fromkeys(['chromosome_verified','assembly_supported_association','same_contig_mge','interval_overlap'], 0)
    for row in order:
        selected = [r for r in loci if r['accession'] == row['accession']]
        flags = {'chromosome_verified': any(r['location_grade'] == 'chromosomal_verified' for r in selected),
                 'assembly_supported_association': any(r['location_grade'] == 'assembly_supported_association' for r in selected),
                 'same_contig_mge': any(r['MGE_same_contig'] == 'yes' for r in selected),
                 'interval_overlap': any(r['MGE_interval_overlap'] == 'yes' for r in selected)}
        assert len(selected) == int(row['n_loci'])
        for key, value in flags.items():
            assert value == (row[key] == 'TRUE'), (row['accession'], key)
            totals[key] += value
    assert list(totals.values()) == [3,1,22,0]
    report = {'status':'PASS','directional_ANI_records':len(raw),'single_reference_membership_pass':152,
              'perflava_single_reference_min_mean_ANI':min(perflava),
              'tetM_carriers':38,'tetM_loci':57,'included_sequences':52,
              'carrier_evidence_counts':totals,
              'scope':'Frozen-result reconciliation; no new prediction, filtering, or dereplication.'}
    output = project.parent / 'revision2_checks.json'
    output.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(report,indent=2))


if __name__ == '__main__':
    main()
