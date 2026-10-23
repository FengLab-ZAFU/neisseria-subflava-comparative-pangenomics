"""Common helpers for boundary8 analysis scripts."""
from pathlib import Path
import csv, os, json

P = str(Path(os.environ["NEISSERIA_PROJECT_ROOT"]) / 'boundary8_split_handoff_20261019')
R = str(Path(os.environ["NEISSERIA_PROJECT_ROOT"]) / 'results_revision_20261015_main152_species')
S = str(Path(os.environ["NEISSERIA_PROJECT_ROOT"]) / 'results_revision_20261019_boundary8_analysis')

def read_tsv(path):
    with open(path, newline='', encoding='utf-8') as fh:
        return list(csv.DictReader(fh, delimiter='\t'))

def write_tsv(path, rows, fieldnames):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'w', newline='', encoding='utf-8') as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, delimiter='\t', lineterminator='\n')
        w.writeheader()
        w.writerows(rows)

def load_json(path):
    with open(path) as fh:
        return json.load(fh)

def write_json(path, obj):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'w') as fh:
        json.dump(obj, fh, indent=2)
