#!/usr/bin/env python3
"""Build small synthetic assemblies with known marker content for testing dryqc_assembly.py.

Genes are back-translated panel proteins (E. coli) or proteins given in a JSON file, placed between
random spacers with a ribosome-binding site. Each case changes one thing relative to the baseline.
"""
import json, os, random, sys

CODON = {'A': ['GCG', 'GCC'], 'R': ['CGC', 'CGT'], 'N': ['AAC', 'AAT'], 'D': ['GAT', 'GAC'], 'C': ['TGC', 'TGT'],
         'Q': ['CAG', 'CAA'], 'E': ['GAA', 'GAG'], 'G': ['GGC', 'GGT'], 'H': ['CAT', 'CAC'], 'I': ['ATT', 'ATC'],
         'L': ['CTG', 'TTA'], 'K': ['AAA', 'AAG'], 'M': ['ATG'], 'F': ['TTT', 'TTC'], 'P': ['CCG', 'CCA'],
         'S': ['AGC', 'TCT'], 'T': ['ACC', 'ACG'], 'W': ['TGG'], 'Y': ['TAT', 'TAC'], 'V': ['GTG', 'GTT'],
         'X': ['GCG']}
COMP = str.maketrans('ACGT', 'TGCA')


def rc(s):
    return s.translate(COMP)[::-1]


def back(aa, rnd):
    aa = aa.rstrip('*')
    if not aa.startswith('M'):
        aa = 'M' + aa[1:]
    return ''.join(rnd.choice(CODON.get(a, ['GCG'])) for a in aa) + 'TAA'


def spacer(n, rnd):
    # AT-rich random spacer with stops in all frames to limit long spurious ORFs
    out = []
    while len(out) * 3 < n:
        out.append(rnd.choice(['TAA', 'TGA', 'TAG', 'ATT', 'AAT', 'TTA', 'GCA', 'ACT', 'TCA', 'CTA']))
    return ''.join(out)[:n]


def block(cds, rnd, strand=1):
    # stops in all three frames, then RBS and a linker without start codons
    unit = 'TAGCTAGCTAG' + 'AGGAGG' + 'AAACAAC' + cds
    if strand == -1:
        unit = rc(unit)
    return spacer(rnd.randint(150, 300), rnd) + unit


def panel_proteins(panel_faa, assembly='GCF_000005845.2'):
    """first record per gene from one E. coli K-12 panel assembly (falls back to any E. coli)."""
    fam2gene = {}
    seqs, name = {}, None
    for line in open(panel_faa):
        if line.startswith('>'):
            name = line[1:].strip()
            seqs[name] = ''
        else:
            seqs[name] += line.strip()
    out = {}
    for n, s in seqs.items():
        fam, sp, acc, pid = n.split('|')
        g = ('rpl' if fam[6] == 'l' else 'rps') + fam[7]
        if 'Escherichia_coli' in sp and (g not in out or acc == assembly):
            out[g] = s
    return out


def build(out_dir, panel_faa, extra_json=None, seed=7):
    os.makedirs(out_dir, exist_ok=True)
    ec = panel_proteins(panel_faa)
    extra = json.load(open(extra_json)) if extra_json else {}
    genes = sorted(ec)
    assert len(genes) == 27, len(genes)
    cases = {}

    def assemble(gene_prot, mods=None, name='x', seed_=seed):
        rnd = random.Random(seed_)
        mods = mods or {}
        cds = {g: back(p, rnd) for g, p in gene_prot.items()}
        contigs, cur = [], ''
        for i, g in enumerate(genes):
            if g in mods.get('delete', []):
                continue
            if g in mods.get('split', []):
                c = cds[g]
                mid = (len(c) // 2) // 3 * 3
                cur += block(c[:mid], rnd)
                contigs.append(cur)                      # contig ends inside the gene
                cur = c[mid:] + spacer(200, rnd)          # next contig starts inside the gene
                continue
            if g in mods.get('truncate_at_end', []):
                # gene placed at the start of a contig, missing its first 40% (no start codon, partial)
                contigs.append(cur + spacer(300, rnd))
                c = cds[g]
                cut = int(len(c) * 0.4) // 3 * 3
                cur = c[cut:] + spacer(200, rnd)          # contig starts inside the gene
                continue
            cur += block(cds[g], rnd, strand=-1 if i % 4 == 3 else 1)
            if len(cur) > 6000:
                contigs.append(cur + spacer(200, rnd))
                cur = ''
        for g, prot in mods.get('add', []):
            cur += block(back(prot, rnd), rnd)
        for g in mods.get('identical_copy', []):
            cur += block(cds[g], rnd)
        for g in mods.get('add_cds_suffix', []):
            c = cds[g]
            cut = int(len(c) // 3 * 0.3) * 3
            cur += block('ATG' + c[cut + 3:], rnd)
        if cur:
            contigs.append(cur + spacer(200, rnd))
        path = os.path.join(out_dir, f'{name}.fna')
        with open(path, 'w') as fh:
            for i, c in enumerate(contigs, 1):
                fh.write(f'>{name}_ctg{i}_length_{len(c)}_cov_{50 if i < 99 else 5}.0\n')
                for k in range(0, len(c), 80):
                    fh.write(c[k:k + 80] + '\n')
        return path

    cases['baseline'] = assemble(ec, name='baseline')
    cases['deletion_rpsD'] = assemble(ec, {'delete': ['rpsD']}, 'deletion_rpsD')
    cases['identical_rplB'] = assemble(ec, {'identical_copy': ['rplB']}, 'identical_rplB')
    cases['split_rplC'] = assemble(ec, {'split': ['rplC']}, 'split_rplC')
    cases['truncated_rplK'] = assemble(ec, {'truncate_at_end': ['rplK']}, 'truncated_rplK')
    long_rpsI = ec['rpsI'] + 'AKEAAKKAAEEAKKAAEAPKAEKKAAPAAEKKAEAPKAEAKKEAAPKAEKKAEAAKKAEEAPKAEKK'
    cases['long_rpsI'] = assemble(dict(ec, rpsI=long_rpsI), name='long_rpsI')
    if extra:
        cd = extra.get('GCF_001457455.1', {})
        if 'rpsF' in cd:
            cases['divergent_rpsF'] = assemble(dict(ec, rpsF=cd['rpsF']), name='divergent_rpsF')
        cp = extra.get('GCF_000013285.1', {})
        if 'rpsD' in cp:
            cases['second_rpsD_Clostridium'] = assemble(ec, {'add': [('rpsD', cp['rpsD'])]}, 'second_rpsD_Clostridium')
        # weak copy: E. coli rpsD with 12% of residues changed (same length)
        rnd = random.Random(3)
        s = list(ec['rpsD'])
        for k in rnd.sample(range(1, len(s)), int(len(s) * 0.12)):
            s[k] = rnd.choice('AVLIKREDST')
        cases['weak_divergent_rpsD'] = assemble(ec, {'add': [('rpsD', ''.join(s))]}, 'weak_divergent_rpsD')
        # weak identical copy: identical C-terminal 70% of rpsD with its own start codon
        suf = ec['rpsD'][int(len(ec['rpsD']) * 0.3):]
        cases['weak_identical_rpsD'] = assemble(ec, {'add_cds_suffix': ['rpsD']}, 'weak_identical_rpsD')
        # all 27 genes duplicated with a divergent second genome (simulated mixture)
        mix = [(g, cp[g]) for g in genes if g in cp]
        cases['mixture_all'] = assemble(ec, {'add': mix}, 'mixture_all')
    with open(os.path.join(out_dir, 'manifest.tsv'), 'w') as fh:
        fh.write('sample\tdataset\tpath\n')
        for k, v in cases.items():
            fh.write(f'{k}\tsynthetic\t{os.path.abspath(v)}\n')
    return cases


if __name__ == '__main__':
    out, panel = sys.argv[1], sys.argv[2]
    extra = sys.argv[3] if len(sys.argv) > 3 else None
    for k, v in build(out, panel, extra).items():
        print(k, v)
