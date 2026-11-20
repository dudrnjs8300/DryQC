#!/usr/bin/env python3
"""Synthetic tests for dryqc_assembly.py. Usage: python run_tests.py <work_dir> [threads]
Needs prodigal, diamond and hmmsearch on PATH. Exit code 0 when all checks pass."""
import csv, json, os, subprocess, sys
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
D = os.path.join(ROOT, 'data')
work = os.path.abspath(sys.argv[1]); thr = sys.argv[2] if len(sys.argv) > 2 else '2'
os.makedirs(work, exist_ok=True)
db = os.path.join(work, 'panel')
subprocess.run(['diamond', 'makedb', '--in', os.path.join(D, 'panel.faa'), '-d', db, '--quiet'], check=True)
sys.path.insert(0, HERE)
import make_synthetic
extra = os.path.join(D, 'test_proteins.json')
cases = make_synthetic.build(os.path.join(work, 'syn'), os.path.join(D, 'panel.faa'), extra if os.path.exists(extra) else None)
out = os.path.join(work, 'out')
for s, path in cases.items():
    subprocess.run([sys.executable, os.path.join(ROOT, 'scripts', 'dryqc_assembly.py'), '--sample', s, '--assembly', path,
                    '--out', out, '--panel-db', db + '.dmnd', '--hmm', os.path.join(D, 'dryqc_profiles.hmm'),
                    '--ranges', os.path.join(D, 'expected_length_ranges.tsv'), '--threads', thr],
                   check=True, capture_output=True)
def genes(s):
    return {r['gene']: r for r in csv.DictReader(open(os.path.join(out, f'{s}.genes.tsv')), delimiter='\t')}
def summ(s):
    return json.load(open(os.path.join(out, f'{s}.sample.json')))
checks = []
def chk(name, cond):
    checks.append((name, bool(cond))); print(('PASS ' if cond else 'FAIL ') + name)
b = genes('baseline'); chk('baseline: 27 single genes (both rule sets)', all(r['v2_state'] == 'PASS_SC' and r['v3_state'] == 'DETECTED_SINGLE' for r in b.values()))
g = genes('deletion_rpsD')['rpsD']; chk('deletion: rpsD missing', g['v2_state'] == 'FAIL_MISS' and g['v3_state'] == 'UNDETECTED_MISSING')
g = genes('identical_rplB')['rplB']; chk('identical copy: WARN_DUP / multi, duplicated', g['v2_state'] == 'WARN_DUP' and g['v3_state'] == 'DETECTED_MULTI' and g['dup_counted'] == 'True')
g = genes('split_rplC')['rplC']; chk('split gene: PASS_FRAG / split, not duplicated', g['v2_state'] == 'PASS_FRAG' and g['v3_state'] == 'DETECTED_SPLIT' and g['dup_counted'] == 'False')
g = genes('truncated_rplK')['rplK']; chk('truncated gene: not detected, reason recorded', g['v3_detected'] == 'False' and g['v3_state'] in ('UNDETECTED_TRUNCATED', 'UNDETECTED_SHORT'))
g = genes('long_rpsI')['rpsI']; chk('long lineage form: BADLEN -> profile detected', g['v2_state'] == 'FAIL_BADLEN' and g['v3_state'] == 'DETECTED_PROFILE')
if 'divergent_rpsF' in cases:
    g = genes('divergent_rpsF')['rpsF']; chk('divergent gene: not eligible -> profile detected', g['v2_state'] in ('AMBIGUOUS', 'FAIL_MISS') and g['v3_state'] == 'DETECTED_PROFILE')
    g = genes('second_rpsD_Clostridium')['rpsD']; chk('second, different rpsD: detected, duplicated', g['v3_detected'] == 'True' and g['dup_counted'] == 'True')
    g = genes('weak_identical_rpsD')['rpsD']; chk('weak identical partial copy: detected, not duplicated', g['v3_state'] == 'DETECTED_MULTI' and g['dup_counted'] == 'False')
    s = summ('mixture_all'); chk('distant second genome: duplication flag', s['v3']['duplication_flag'] and s['v3']['duplicated_genes'] >= 25)
chk('baseline: score 100, no flag', summ('baseline')['v3']['score'] == 100 and not summ('baseline')['v3']['duplication_flag'])
# stored-protein mode: proteins and CDS of a separate Prodigal run, no assembly
pm = os.path.join(work, 'protein_mode')
os.makedirs(pm, exist_ok=True)
same = True
for s in [x for x in ('baseline', 'identical_rplB', 'second_rpsD_Clostridium', 'mixture_all') if x in cases]:
    faa, fna = os.path.join(pm, s + '.faa'), os.path.join(pm, s + '.fna')
    subprocess.run(['prodigal', '-i', cases[s], '-a', faa, '-d', fna, '-p', 'meta', '-q', '-o', os.devnull],
                   check=True, capture_output=True)
    subprocess.run([sys.executable, os.path.join(ROOT, 'scripts', 'dryqc_assembly.py'), '--sample', s, '--proteins', faa,
                    '--cds', fna, '--out', os.path.join(pm, 'out'), '--panel-db', db + '.dmnd',
                    '--hmm', os.path.join(D, 'dryqc_profiles.hmm'), '--ranges', os.path.join(D, 'expected_length_ranges.tsv'),
                    '--threads', thr], check=True, capture_output=True)
    a_ = genes(s)
    b_ = {r['gene']: r for r in csv.DictReader(open(os.path.join(pm, 'out', f'{s}.genes.tsv')), delimiter='\t')}
    same &= all((a_[g]['v2_state'], a_[g]['v3_state'], a_[g]['dup_counted']) ==
                (b_[g]['v2_state'], b_[g]['v3_state'], b_[g]['dup_counted']) for g in a_)
chk('stored-protein mode: same gene states as the assembly run', same)
n = sum(1 for _, ok in checks if not ok)
print(f'{len(checks) - n}/{len(checks)} checks passed')
sys.exit(1 if n else 0)
