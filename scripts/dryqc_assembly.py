#!/usr/bin/env python3
"""DryQC assembly branch: earlier rules (Supplementary Table S8, "v2") and final rules ("v3").

One assembly in, three tables out:
  <out>/<sample>.genes.tsv   one row per marker gene, with the v2 state and the v3 state side by side
  <out>/<sample>.loci.tsv    one row per candidate locus (eligible panel hits, profile hits)
  <out>/<sample>.sample.json sample summary for both rule sets and the software versions

Earlier rules (v2) follow Supplementary Table S8 exactly, so that the earlier gene states can be
reproduced and compared (regression check). Final rules (v3) are the decisions of 2026-10-09/10:
  rule 1  copies are counted per locus (overlapping proteins on one contig are one locus)
  rule 2  weak copy = both coverages >=0.6 and bitscore >=50% of the best, not strong
  rule 3  a weak copy counts as duplicated when its CDS differs from the best copy
          (nucleotide identity <98% over >=60% of the shorter CDS)
  rule 4  a second strong copy counts as duplicated whatever its sequence
  rule 5  eligible primary + weak second candidate (earlier AMBIGUOUS) -> detected
  rule 6  genes not detected with the panel are searched with the profile HMMs
          (domain i-Evalue <=1e-5, model coverage >=0.65); a complete ORF whose length lies in
          the expected range of the gene is detected
  rule 7  expected range = type-strain minimum x0.9 to maximum x1.1 (data/expected_length_ranges.tsv)
  rule 8  copies found only with the profiles are recorded, not counted as duplicated
  rule 9  undetected genes carry a reason: missing / residual / truncated / short / long
  rule 10 score = max[0, 100 x (C + 0.5F)/27 - U]; copies do not change it
  rule 11 state names: DETECTED_* or UNDETECTED_*
  rule 14 sample passes with >=24 detected genes; duplication flag with >=25 duplicated genes

Stored-protein mode (--proteins, 2026-10-13): when the assembly FASTA is no longer available, the proteins
that Prodigal (-p meta) predicted in the original run are read from their FASTA with Prodigal headers, and
the same rules are applied. Coding sequences for the weak-copy test are taken from a second Prodigal CDS
file (--cds, e.g. BUSCO's predicted.fna) when that file predicts the same gene (same contig, strand and stop
position, translation verified); weak copies without a coding sequence are not counted as duplicated, and the
six-frame scan of undetected genes is not possible.

Requirements: prodigal 2.6.3, diamond (one version for all runs), hmmsearch (HMMER 3.x),
python >=3.8 with biopython. pyhmmer is not required.
"""
import argparse, csv, gzip, hashlib, io, json, os, re, shutil, subprocess, sys, tempfile, time
from collections import defaultdict

MODULE_VERSION = "dryqc_asm_v3 1.0 (2026-10-12)"

GENES = ['rplA', 'rplB', 'rplC', 'rplD', 'rplE', 'rplF', 'rplK', 'rplL', 'rplM', 'rplN', 'rplP', 'rplR', 'rplS',
         'rplT', 'rplV', 'rpsB', 'rpsC', 'rpsD', 'rpsE', 'rpsF', 'rpsG', 'rpsH', 'rpsI', 'rpsJ', 'rpsK', 'rpsL', 'rpsM']
FAMILY_TO_GENE = {('bDBv1_' + ('l' if g.startswith('rpl') else 's') + g[3]): g for g in GENES}

# ---- earlier rules (Supplementary Table S8) ----
MIN_ID = 0.40          # DIAMOND pident / 100
MIN_QCOV = 0.40        # alignment length / query length
LEN_TOL = 0.40         # query length within +-40% of the matched reference length
STRONG_COV = 0.80      # query and target coverage
STRONG_RATIO = 0.75    # bitscore / best bitscore
INTER_COV = 0.60
INTER_RATIO = 0.50
FRAG_UNION = 0.80
RESID = 0.32
N_DETECT = 24
# ---- final rules ----
CDS_ID_DIFF = 0.98
CDS_MIN_ALN = 0.60
HMM_E = 1e-5
HMM_COV = 0.65
N_DUPFLAG = 25
ORF_SCAN_MIN_AA = 40


def log(msg):
    sys.stderr.write(time.strftime('%H:%M:%S ') + msg + '\n')
    sys.stderr.flush()


def sha256(s):
    return hashlib.sha256(s.encode()).hexdigest()


def open_text(path):
    if path.endswith('.gz'):
        return io.TextIOWrapper(gzip.open(path, 'rb'))
    return open(path)


# ---------------------------------------------------------------- input
def read_assembly(path):
    """Return list of (contig_id, sequence). FASTA or GenBank flat file, optionally gzipped."""
    low = path.lower()
    is_gb = any(low.endswith(x) for x in ('.gbff', '.gbk', '.gb', '.gbff.gz', '.gbk.gz', '.gb.gz'))
    seqs = []
    if is_gb:
        from Bio import SeqIO
        with open_text(path) as fh:
            for rec in SeqIO.parse(fh, 'genbank'):
                seqs.append((rec.id, str(rec.seq).upper()))
    else:
        name, buf = None, []
        with open_text(path) as fh:
            for line in fh:
                if line.startswith('>'):
                    if name is not None:
                        seqs.append((name, ''.join(buf).upper()))
                    name = line[1:].strip().split()[0]
                    buf = []
                else:
                    buf.append(line.strip())
            if name is not None:
                seqs.append((name, ''.join(buf).upper()))
    if not seqs:
        raise ValueError(f'no sequences in {path}')
    return seqs


def contig_coverage(name):
    m = re.search(r'_cov_([0-9.]+)', name)
    return float(m.group(1)) if m else None


# ---------------------------------------------------------------- external tools
def tool_version(cmd):
    try:
        out = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
        txt = (out.stdout + out.stderr).strip().splitlines()
        for line in txt:
            if 'HMMER' in line:
                return line.strip('# ').strip()
        return txt[0] if txt else ''
    except Exception as e:
        return f'unavailable ({e})'


def run(cmd, log_path=None):
    res = subprocess.run(cmd, capture_output=True, text=True)
    if log_path:
        with open(log_path, 'a') as fh:
            fh.write(' '.join(cmd) + '\n' + res.stdout + res.stderr + '\n')
    if res.returncode != 0:
        raise RuntimeError(f"command failed ({res.returncode}): {' '.join(cmd)}\n{res.stderr[-2000:]}")
    return res


def run_prodigal(fasta, workdir, prodigal='prodigal'):
    faa = os.path.join(workdir, 'proteins.faa')
    fna = os.path.join(workdir, 'genes.fna')
    gff = os.path.join(workdir, 'genes.gff')
    run([prodigal, '-i', fasta, '-a', faa, '-d', fna, '-f', 'gff', '-o', gff, '-p', 'meta', '-q'],
        os.path.join(workdir, 'tools.log'))
    return faa, fna


def parse_prodigal(faa, fna):
    prots, order = {}, []
    name = None
    for line in open(faa):
        if line.startswith('>'):
            parts = [p.strip() for p in line[1:].split('#')]
            name = parts[0]
            start, end, strand = int(parts[1]), int(parts[2]), int(parts[3])
            attrs = dict(kv.split('=', 1) for kv in parts[4].split(';') if '=' in kv)
            contig = name.rsplit('_', 1)[0]
            prots[name] = dict(id=name, contig=contig, start=start, end=end, strand=strand,
                               partial=attrs.get('partial', '00'), start_type=attrs.get('start_type', ''),
                               aa='', cds='')
            order.append(name)
        else:
            prots[name]['aa'] += line.strip()
    name = None
    for line in open(fna):
        if line.startswith('>'):
            name = line[1:].split()[0]
        elif name in prots:
            prots[name]['cds'] += line.strip().upper()
    for p in prots.values():
        p['aa'] = p['aa'].rstrip('*')
        p['len'] = len(p['aa'])
    return prots, order


CODON = {}
_b = 'TCAG'
_aa = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
for _i, _x in enumerate(_b):
    for _j, _y in enumerate(_b):
        for _k, _z in enumerate(_b):
            CODON[_x + _y + _z] = _aa[16 * _i + 4 * _j + _k]


def translate(cds):
    return ''.join(CODON.get(cds[i:i + 3], 'X') for i in range(0, len(cds) - 2, 3))


def parse_stored_proteins(faa):
    """Proteins of an earlier Prodigal run (FASTA with Prodigal headers); no coding sequences."""
    prots, order = {}, []
    name = None
    with open_text(faa) as fh:
        for line in fh:
            if line.startswith('>'):
                parts = [x.strip() for x in line[1:].split('#')]
                if len(parts) < 5:
                    raise ValueError(f'{faa}: header without Prodigal coordinates: {line.strip()[:80]}')
                name = parts[0]
                start, end, strand = int(parts[1]), int(parts[2]), int(parts[3])
                attrs = dict(kv.split('=', 1) for kv in parts[4].split(';') if '=' in kv)
                prots[name] = dict(id=name, contig=name.rsplit('_', 1)[0], start=start, end=end, strand=strand,
                                   partial=attrs.get('partial', '00'), start_type=attrs.get('start_type', ''),
                                   aa='', cds='')
                order.append(name)
            elif name is not None:
                prots[name]['aa'] += line.strip()
    for x in prots.values():
        x['aa'] = x['aa'].rstrip('*')
        x['len'] = len(x['aa'])
    if not prots:
        raise ValueError(f'no proteins in {faa}')
    return prots, order


def attach_cds(prots, cds_fna):
    """Coding sequences from another Prodigal CDS file, matched by contig, strand and stop position.
    The CDS is trimmed at its 5' end when that prediction starts upstream, and kept only if its translation
    equals the protein (first residue ignored). Returns the number of proteins given a CDS."""
    ref = {}
    name, buf = None, []

    def keep(nm, seq):
        parts = [x.strip() for x in nm.split('#')]
        if len(parts) < 4:
            return
        cid = parts[0].rsplit('_', 1)[0]
        st, en, sd = int(parts[1]), int(parts[2]), int(parts[3])
        stop = en if sd == 1 else st
        ref[(cid, sd, stop)] = (st, en, seq.upper())
    with open_text(cds_fna) as fh:
        for line in fh:
            if line.startswith('>'):
                if name is not None:
                    keep(name, ''.join(buf))
                name, buf = line[1:].strip(), []
            else:
                buf.append(line.strip())
        if name is not None:
            keep(name, ''.join(buf))
    n = 0
    for x in prots.values():
        if x['cds']:
            continue
        stop = x['end'] if x['strand'] == 1 else x['start']
        hit = ref.get((x['contig'], x['strand'], stop))
        if hit is None:
            continue
        st, en, seq = hit
        trim = (x['start'] - st) if x['strand'] == 1 else (en - x['end'])
        if trim < 0 or trim % 3:
            continue
        cds = seq[trim:]
        aa = translate(cds).rstrip('*')
        if len(aa) == len(x['aa']) and aa[1:] == x['aa'][1:]:
            x['cds'] = cds
            n += 1
    return n


def write_proteins(prots, order, path):
    with open(path, 'w') as fh:
        for n in order:
            aa = prots[n]['aa'].replace('*', 'X')
            if aa:
                fh.write(f'>{n}\n{aa}\n')


DMND_FIELDS = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend',
               'sstart', 'send', 'evalue', 'bitscore', 'qlen', 'slen']


def run_diamond(prot_faa, db, workdir, threads=1, diamond='diamond'):
    out = os.path.join(workdir, 'diamond.tsv')
    run([diamond, 'blastp', '-q', prot_faa, '-d', db, '-o', out, '--more-sensitive', '--max-target-seqs', '20',
         '--evalue', '1e-5', '--threads', str(threads), '--outfmt', '6'] + DMND_FIELDS,
        os.path.join(workdir, 'tools.log'))
    hits = []
    for row in csv.reader(open(out), delimiter='\t'):
        if not row:
            continue
        h = dict(zip(DMND_FIELDS, row))
        for k in ('pident', 'evalue', 'bitscore'):
            h[k] = float(h[k])
        for k in ('length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'qlen', 'slen'):
            h[k] = int(h[k])
        fam = h['sseqid'].split('|')[0]
        h['gene'] = FAMILY_TO_GENE.get(fam)
        if h['gene'] is None:
            continue
        parts = h['sseqid'].split('|')
        h['panel_species'] = parts[1] if len(parts) > 1 else ''
        h['panel_assembly'] = parts[2] if len(parts) > 2 else ''
        h['id'] = h['pident'] / 100.0
        h['qcov'] = h['length'] / h['qlen'] if h['qlen'] else 0.0
        h['tcov'] = h['length'] / h['slen'] if h['slen'] else 0.0
        h['len_ok'] = abs(h['qlen'] - h['slen']) <= LEN_TOL * h['slen']
        h['eligible'] = h['id'] >= MIN_ID and h['qcov'] >= MIN_QCOV and h['len_ok']
        hits.append(h)
    return hits


def run_hmmsearch(prot_faa, hmm, workdir, threads=1, hmmsearch='hmmsearch'):
    dom = os.path.join(workdir, 'hmm.domtbl')
    run([hmmsearch, '--noali', '-o', os.devnull, '--domtblout', dom, '-E', '1', '--domE', '1',
         '--cpu', str(threads), hmm, prot_faa], os.path.join(workdir, 'tools.log'))
    return parse_domtbl(dom)


def parse_domtbl(path):
    """best domain per (protein, profile): i-Evalue, model coverage, score."""
    best = {}
    for line in open(path):
        if line.startswith('#'):
            continue
        f = line.split()
        if len(f) < 22:
            continue
        tname, tlen, qname, qlen = f[0], int(f[2]), f[3], int(f[5])
        i_e, dscore = float(f[12]), float(f[13])
        hfrom, hto = int(f[15]), int(f[16])
        cov = (hto - hfrom + 1) / qlen
        key = (tname, qname)
        cur = best.get(key)
        if cur is None or dscore > cur['score']:
            best[key] = dict(protein=tname, gene=qname, i_evalue=i_e, score=dscore, model_cov=cov,
                             hmm_from=hfrom, hmm_to=hto, model_len=qlen)
    return best


# ---------------------------------------------------------------- CDS identity
_ALIGNER = None


def _aligner():
    """protein aligner: BLOSUM62, gap open -11, extend -1, free end gaps."""
    global _ALIGNER
    if _ALIGNER is None:
        from Bio import Align
        from Bio.Align import substitution_matrices
        al = Align.PairwiseAligner()
        al.mode = 'global'
        al.substitution_matrix = substitution_matrices.load('BLOSUM62')
        al.open_gap_score, al.extend_gap_score = -11, -1
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            al.target_end_gap_score = 0.0   # end_insertion_score in recent Biopython
            al.query_end_gap_score = 0.0    # end_deletion_score in recent Biopython
        _ALIGNER = al
    return _ALIGNER


def cds_identity(aa_a, cds_a, aa_b, cds_b):
    """Nucleotide identity of two coding sequences along their protein alignment (codon by codon).
    Returns (identity, aligned fraction of the shorter protein).
    identity = identical nucleotides / nucleotide columns of the aligned region (internal gaps count
    as non-identical columns); fraction = aligned span on the shorter protein / its length."""
    if not aa_a or not aa_b or not cds_a or not cds_b:
        return None, 0.0
    if cds_a == cds_b:
        return 1.0, 1.0
    a = aa_a.replace('*', 'X')
    b = aa_b.replace('*', 'X')
    aln = _aligner().align(a, b)[0]
    blocks = list(zip(*aln.aligned))
    if not blocks:
        return None, 0.0
    ident = cols = 0
    for (ts, te), (qs, qe) in blocks:
        for k in range(te - ts):
            ca = cds_a[3 * (ts + k):3 * (ts + k) + 3]
            cb = cds_b[3 * (qs + k):3 * (qs + k) + 3]
            ident += sum(1 for x, y in zip(ca, cb) if x == y)
            cols += 3
    for i in range(1, len(blocks)):
        cols += 3 * ((blocks[i][0][0] - blocks[i - 1][0][1]) + (blocks[i][1][0] - blocks[i - 1][1][1]))
    span_a = blocks[-1][0][1] - blocks[0][0][0]
    span_b = blocks[-1][1][1] - blocks[0][1][0]
    frac = span_a / len(a) if len(a) <= len(b) else span_b / len(b)
    return (ident / cols if cols else None), frac


# ---------------------------------------------------------------- earlier rules (v2)
def best_per_query(hits):
    """one best alignment per query within each family (S8 ordering)."""
    by = {}
    for h in hits:
        key = (h['gene'], h['qseqid'])
        cur = by.get(key)
        rank = (h['bitscore'], h['id'], h['qcov'], h['tcov'], h['length'])
        if cur is None or rank > cur[0]:
            by[key] = (rank, h)
    out = defaultdict(list)
    for (g, q), (_, h) in by.items():
        out[g].append(h)
    for g in out:
        out[g].sort(key=lambda h: (h['bitscore'], h['id'], h['qcov'], h['tcov'], h['length']), reverse=True)
    return out


def frag_union(cands):
    """union of subject intervals / target length of the highest-bitscore fragment."""
    if len(cands) < 2:
        return 0.0
    top = cands[0]
    iv = sorted((min(h['sstart'], h['send']), max(h['sstart'], h['send'])) for h in cands)
    tot, cs, ce = 0, None, None
    for s, e in iv:
        if cs is None:
            cs, ce = s, e
        elif s <= ce + 1:
            ce = max(ce, e)
        else:
            tot += ce - cs + 1
            cs, ce = s, e
    tot += ce - cs + 1
    return tot / top['slen']


def classify_v2(qhits):
    """returns state and details for one gene (list of per-query best hits, sorted)."""
    elig = [h for h in qhits if h['eligible']]
    if len(elig) == 1:
        return 'PASS_SC', dict(n_eligible=1)
    if len(elig) >= 2:
        prim = elig[0]
        strong = [h for h in elig if h['qcov'] >= STRONG_COV and h['tcov'] >= STRONG_COV
                  and h['bitscore'] / prim['bitscore'] >= STRONG_RATIO]
        if len(strong) >= 3:
            return 'FAIL_MULTI', dict(n_eligible=len(elig), n_strong=len(strong))
        if len(strong) == 2:
            return 'WARN_DUP', dict(n_eligible=len(elig), n_strong=2)
        inter = [h for h in elig[1:] if h['qcov'] >= INTER_COV and h['tcov'] >= INTER_COV
                 and h['bitscore'] / prim['bitscore'] >= INTER_RATIO]
        if inter:
            return 'AMBIGUOUS', dict(n_eligible=len(elig), n_strong=len(strong), ambiguous_type='second_candidate')
        return 'PASS_SC', dict(n_eligible=len(elig), n_strong=len(strong))
    idq = [h for h in qhits if h['id'] >= MIN_ID]
    if len(idq) >= 2 and frag_union(idq) >= FRAG_UNION:
        return 'PASS_FRAG', dict(n_fragments=len(idq))
    if any(h['id'] >= MIN_ID and h['qcov'] >= MIN_QCOV and not h['len_ok'] for h in qhits):
        return 'FAIL_BADLEN', {}
    if qhits and (max(h['id'] for h in qhits) >= RESID or max(h['qcov'] for h in qhits) >= RESID):
        return 'AMBIGUOUS', dict(ambiguous_type='residual')
    return 'FAIL_MISS', {}


def v2_sample(states):
    c = defaultdict(int)
    for s in states.values():
        c[s] += 1
    S, F, W, B, H, M = (c['PASS_SC'], c['PASS_FRAG'], c['WARN_DUP'], c['FAIL_BADLEN'], c['FAIL_MULTI'],
                        c['FAIL_MISS'])
    detected = S + F + W
    score = max(0.0, 100.0 * (S + 0.5 * F) / 27 - (W + B + H + M))
    return dict(counts=dict(c), detected=detected, status='PASS' if detected >= N_DETECT else 'FAIL',
                score=round(score, 4))


# ---------------------------------------------------------------- final rules (v3)
def group_loci(hits, prots):
    """rule 1: eligible proteins with overlapping coordinates on one contig form one locus."""
    loci = []
    for h in hits:  # sorted by rank
        p = prots[h['qseqid']]
        placed = False
        for L in loci:
            r = prots[L['rep']['qseqid']]
            if r['contig'] == p['contig'] and not (p['end'] < r['start'] or p['start'] > r['end']):
                L['members'].append(h)
                placed = True
                break
        if not placed:
            loci.append(dict(rep=h, members=[h]))
    return loci


def classify_v3(gene, qhits, prots, hmm_by_gene, used, ranges):
    elig = [h for h in qhits if h['eligible']]
    rec = dict(gene=gene, v3_state='', v3_detected=False, v3_split=False, v3_reason='', dup_counted=False,
               n_loci=0, n_strong_other=0, n_weak=0, n_weak_diff=0, weak_cds_identity='', n_other_eligible=0,
               profile_loci=0)
    loci_rows = []
    lo, hi = ranges.get(gene, (None, None))
    if elig:
        loci = group_loci(elig, prots)
        best = loci[0]['rep']
        bp = prots[best['qseqid']]
        rec.update(n_loci=len(loci), best_protein=best['qseqid'], best_len=bp['len'],
                   best_identity=round(best['id'], 4), best_qcov=round(best['qcov'], 4),
                   best_tcov=round(best['tcov'], 4), best_bitscore=best['bitscore'],
                   best_partial=bp['partial'], best_aa_sha256=sha256(bp['aa']), best_cds_sha256=sha256(bp['cds']),
                   best_panel_species=best['panel_species'], best_panel_assembly=best['panel_assembly'])
        cls_list = []
        for i, L in enumerate(loci):
            h = L['rep']
            p = prots[h['qseqid']]
            ratio = h['bitscore'] / best['bitscore']
            if i == 0:
                cls, ident, frac = 'best', 1.0, 1.0
            elif h['qcov'] >= STRONG_COV and h['tcov'] >= STRONG_COV and ratio >= STRONG_RATIO:
                cls, ident, frac = 'strong', None, None
            elif h['qcov'] >= INTER_COV and h['tcov'] >= INTER_COV and ratio >= INTER_RATIO:
                cls, ident, frac = 'weak', None, None
            else:
                cls, ident, frac = 'other_eligible', None, None
            if cls in ('strong', 'weak'):
                ident, frac = cds_identity(bp['aa'], bp['cds'], p['aa'], p['cds'])
            cls_list.append((cls, ident, frac))
            loci_rows.append(dict(gene=gene, source='panel', locus_rank=i + 1, copy_class=cls, protein=h['qseqid'],
                                  contig=p['contig'], start=p['start'], end=p['end'], strand=p['strand'],
                                  partial=p['partial'], aa_len=p['len'], identity=round(h['id'], 4),
                                  qcov=round(h['qcov'], 4), tcov=round(h['tcov'], 4), bitscore=h['bitscore'],
                                  bitscore_ratio=round(ratio, 4),
                                  cds_identity_to_best='' if ident is None else round(ident, 5),
                                  cds_aligned_fraction='' if frac is None else round(frac, 4),
                                  contig_coverage=contig_coverage(p['contig']) or '',
                                  panel_species=h['panel_species'], aa_sha256=sha256(p['aa']),
                                  members=';'.join(m['qseqid'] for m in L['members'])))
        strong_o = [c for c in cls_list[1:] if c[0] == 'strong']
        weak = [c for c in cls_list[1:] if c[0] == 'weak']
        weak_diff = [c for c in weak if c[1] is not None and c[2] >= CDS_MIN_ALN and c[1] < CDS_ID_DIFF]
        rec.update(n_strong_other=len(strong_o), n_weak=len(weak), n_weak_diff=len(weak_diff),
                   n_other_eligible=sum(1 for c in cls_list[1:] if c[0] == 'other_eligible'),
                   weak_cds_identity=';'.join('' if c[1] is None else f'{c[1]:.4f}/{c[2]:.2f}' for c in weak))
        rec['dup_counted'] = bool(strong_o or weak_diff)
        rec['v3_state'] = 'DETECTED_MULTI' if (strong_o or weak) else 'DETECTED_SINGLE'
        rec['v3_detected'] = True
        for q in [m['qseqid'] for L in loci for m in L['members']]:
            used.add(q)
        return rec, loci_rows
    idq = [h for h in qhits if h['id'] >= MIN_ID]
    if len(idq) >= 2 and frag_union(idq) >= FRAG_UNION:
        rec.update(v3_state='DETECTED_SPLIT', v3_detected=True, v3_split=True, n_loci=1,
                   best_protein=idq[0]['qseqid'], best_identity=round(idq[0]['id'], 4))
        for h in idq:
            used.add(h['qseqid'])
            p = prots[h['qseqid']]
            loci_rows.append(dict(gene=gene, source='panel', locus_rank=1, copy_class='fragment', protein=h['qseqid'],
                                  contig=p['contig'], start=p['start'], end=p['end'], strand=p['strand'],
                                  partial=p['partial'], aa_len=p['len'], identity=round(h['id'], 4),
                                  qcov=round(h['qcov'], 4), tcov=round(h['tcov'], 4), bitscore=h['bitscore']))
        return rec, loci_rows
    # rule 6: profile search
    cands = [c for c in hmm_by_gene.get(gene, []) if c['protein'] not in used]
    good = [c for c in cands if c['i_evalue'] <= HMM_E and c['model_cov'] >= HMM_COV
            and prots[c['protein']]['partial'] == '00' and lo is not None
            and lo <= prots[c['protein']]['len'] <= hi]
    good.sort(key=lambda c: c['score'], reverse=True)
    if good:
        b = good[0]
        p = prots[b['protein']]
        panel_note = 'no panel hit'
        ph = [h for h in qhits if h['qseqid'] == b['protein']]
        if ph:
            h = ph[0]
            if h['id'] < MIN_ID or h['qcov'] < MIN_QCOV:
                panel_note = 'below panel similarity criteria'
            elif not h['len_ok']:
                panel_note = 'outside panel length tolerance'
        rec.update(v3_state='DETECTED_PROFILE', v3_detected=True, v3_reason=panel_note, n_loci=1,
                   profile_loci=len(good), best_protein=b['protein'], best_len=p['len'], best_partial=p['partial'],
                   best_aa_sha256=sha256(p['aa']), best_cds_sha256=sha256(p['cds']),
                   profile_i_evalue=b['i_evalue'], profile_model_cov=round(b['model_cov'], 3),
                   best_identity=round(ph[0]['id'], 4) if ph else '')
        for c in good:
            used.add(c['protein'])
            q = prots[c['protein']]
            loci_rows.append(dict(gene=gene, source='profile', locus_rank=good.index(c) + 1,
                                  copy_class='profile_best' if c is b else 'profile_copy', protein=c['protein'],
                                  contig=q['contig'], start=q['start'], end=q['end'], strand=q['strand'],
                                  partial=q['partial'], aa_len=q['len'], profile_i_evalue=c['i_evalue'],
                                  profile_model_cov=round(c['model_cov'], 3), aa_sha256=sha256(q['aa'])))
        return rec, loci_rows
    # rule 9: reason for non-detection (profile candidates first, then panel hits that failed length)
    prof = sorted([c for c in cands if c['i_evalue'] <= HMM_E], key=lambda c: c['score'], reverse=True)
    pan = [h for h in qhits if h['id'] >= MIN_ID and h['qcov'] >= MIN_QCOV and h['qseqid'] not in used]
    pick = None
    if prof:
        pick = (prots[prof[0]['protein']], 'profile')
    elif pan:
        pick = (prots[pan[0]['qseqid']], 'panel')
    if pick:
        p, src = pick
        if p['partial'] != '00':
            st = 'UNDETECTED_TRUNCATED'
        elif lo is not None and p['len'] < lo:
            st = 'UNDETECTED_SHORT'
        elif hi is not None and p['len'] > hi:
            st = 'UNDETECTED_LONG'
        else:
            st = 'UNDETECTED_RESIDUAL'
        rec.update(v3_state=st, best_protein=p['id'], best_len=p['len'], best_partial=p['partial'],
                   v3_reason=f'candidate from {src}')
    elif qhits and (max(h['id'] for h in qhits) >= RESID or max(h['qcov'] for h in qhits) >= RESID):
        rec.update(v3_state='UNDETECTED_RESIDUAL', v3_reason='residual panel similarity')
    else:
        rec.update(v3_state='UNDETECTED_MISSING', v3_reason='no evidence among predicted proteins')
    return rec, loci_rows


def v3_sample(recs):
    c = defaultdict(int)
    for r in recs.values():
        c[r['v3_state']] += 1
    F = sum(1 for r in recs.values() if r['v3_split'])
    C = sum(1 for r in recs.values() if r['v3_detected'] and not r['v3_split'])
    U = 27 - C - F
    dup = sum(1 for r in recs.values() if r['dup_counted'])
    score = max(0.0, 100.0 * (C + 0.5 * F) / 27 - U)
    return dict(counts=dict(c), detected=C + F, C=C, F=F, U=U, status='PASS' if C + F >= N_DETECT else 'FAIL',
                duplicated_genes=dup, duplication_flag=dup >= N_DUPFLAG, score=round(score, 4))


# ---------------------------------------------------------------- DNA ORF diagnostic
STOPS = {'TAA', 'TAG', 'TGA'}


def six_frame_orfs(seqs, min_aa=ORF_SCAN_MIN_AA):
    from Bio.Seq import Seq
    for cid, s in seqs:
        rc = str(Seq(s).reverse_complement())
        for strand, seq in ((1, s), (-1, rc)):
            for fr in range(3):
                sub = seq[fr:]
                sub = sub[:len(sub) - len(sub) % 3]
                aa = str(Seq(sub).translate())
                pos = 0
                for piece in aa.split('*'):
                    if len(piece) >= min_aa:
                        st = fr + pos * 3
                        en = st + len(piece) * 3
                        if strand == 1:
                            a, b = st + 1, en
                        else:
                            a, b = len(s) - en + 1, len(s) - st
                        yield f'{cid}|{strand}|{a}|{b}', piece
                    pos += len(piece) + 1


def orf_diagnostic(seqs, genes, hmm, workdir, hmmsearch='hmmsearch', threads=1):
    """for genes undetected among predicted proteins: does any stop-to-stop ORF in the DNA match the profile?"""
    if not genes:
        return {}
    sub_hmm = os.path.join(workdir, 'sub.hmm')
    keep, buf = False, []
    with open(hmm) as fh, open(sub_hmm, 'w') as out:
        for line in fh:
            buf.append(line)
            if line.startswith('NAME'):
                keep = line.split()[1] in genes
            if line.startswith('//'):
                if keep:
                    out.write(''.join(buf))
                buf, keep = [], False
    orf = os.path.join(workdir, 'orfs.faa')
    with open(orf, 'w') as fh:
        for n, aa in six_frame_orfs(seqs):
            fh.write(f'>{n}\n{aa}\n')
    best = run_hmmsearch(orf, sub_hmm, workdir, threads, hmmsearch)
    out = {}
    for (t, g), d in best.items():
        if d['i_evalue'] <= HMM_E and (g not in out or d['score'] > out[g]['score']):
            out[g] = dict(d, orf=t)
    return out


# ---------------------------------------------------------------- main
def load_ranges(path):
    rng = {}
    for r in csv.DictReader(open(path), delimiter='\t'):
        rng[r['gene']] = (int(r['exp_min']), int(r['exp_max']))
    return rng


def analyse(sample, assembly, outdir, panel_db, hmm, ranges_path, threads=1, keep=False,
            prodigal='prodigal', diamond='diamond', hmmsearch='hmmsearch', orf_scan=True):
    t0 = time.time()
    os.makedirs(outdir, exist_ok=True)
    work = tempfile.mkdtemp(prefix=f'dryqc_{sample}_', dir=outdir)
    try:
        seqs = read_assembly(assembly)
        fasta = os.path.join(work, 'assembly.fna')
        with open(fasta, 'w') as fh:
            for cid, s in seqs:
                fh.write(f'>{cid}\n')
                for i in range(0, len(s), 80):
                    fh.write(s[i:i + 80] + '\n')
        faa, fna = run_prodigal(fasta, work, prodigal)
        prots, order = parse_prodigal(faa, fna)
        base = dict(sample=sample, assembly=os.path.abspath(assembly), input_mode='assembly',
                    contigs=len(seqs), total_bp=sum(len(s) for _, s in seqs))
        return _analyse_core(sample, prots, order, seqs if orf_scan else None, work, outdir, panel_db, hmm,
                             ranges_path, threads, prodigal, diamond, hmmsearch, base, t0)
    finally:
        if keep:
            log(f'kept work directory {work}')
        else:
            shutil.rmtree(work, ignore_errors=True)


def analyse_proteins(sample, faa, outdir, panel_db, hmm, ranges_path, threads=1, keep=False, cds_fna=None,
                     n_scaffolds=None, total_length=None, diamond='diamond', hmmsearch='hmmsearch'):
    """Stored-protein mode: the rules applied to the proteins of an earlier Prodigal run."""
    t0 = time.time()
    os.makedirs(outdir, exist_ok=True)
    work = tempfile.mkdtemp(prefix=f'dryqc_{sample}_', dir=outdir)
    try:
        prots, order = parse_stored_proteins(faa)
        n_cds = attach_cds(prots, cds_fna) if cds_fna else 0
        contigs_with_genes = len(set(x['contig'] for x in prots.values()))
        base = dict(sample=sample, assembly='', proteins_input=os.path.abspath(faa),
                    cds_input=os.path.abspath(cds_fna) if cds_fna else '', input_mode='stored_proteins',
                    contigs=n_scaffolds if n_scaffolds not in (None, '') else contigs_with_genes,
                    total_bp=total_length if total_length not in (None, '') else 0,
                    contigs_with_genes=contigs_with_genes, proteins_with_cds=n_cds)
        return _analyse_core(sample, prots, order, None, work, outdir, panel_db, hmm, ranges_path, threads,
                             None, diamond, hmmsearch, base, t0)
    finally:
        if keep:
            log(f'kept work directory {work}')
        else:
            shutil.rmtree(work, ignore_errors=True)


def _analyse_core(sample, prots, order, seqs, work, outdir, panel_db, hmm, ranges_path, threads, prodigal, diamond,
                  hmmsearch, base, t0):
    clean = os.path.join(work, 'proteins.clean.faa')
    write_proteins(prots, order, clean)
    hits = run_diamond(clean, panel_db, work, threads, diamond)
    qbest = best_per_query(hits)
    ranges = load_ranges(ranges_path)
    hmm_best = run_hmmsearch(clean, hmm, work, threads, hmmsearch)
    # each protein to its best-scoring profile only
    top = {}
    for (pname, g), d in hmm_best.items():
        if pname not in top or d['score'] > top[pname]['score']:
            top[pname] = d
    hmm_by_gene = defaultdict(list)
    for d in top.values():
        hmm_by_gene[d['gene']].append(d)
    # proteins eligible for any gene with the panel are reserved for that gene (v3 profile step)
    used = set(h['qseqid'] for g, hs in qbest.items() for h in hs if h['eligible'])
    v2, v2det, recs, loci = {}, {}, {}, []
    for g in GENES:
        qh = qbest.get(g, [])
        st, det = classify_v2(qh)
        v2[g], v2det[g] = st, det
        r, lr = classify_v3(g, qh, prots, hmm_by_gene, used, ranges)
        r['v2_state'] = st
        r['v2_ambiguous_type'] = det.get('ambiguous_type', '')
        r['v2_n_eligible'] = det.get('n_eligible', 0)
        r['v2_n_strong'] = det.get('n_strong', '')
        lo, hi = ranges.get(g, ('', ''))
        r['exp_min'], r['exp_max'] = lo, hi
        recs[g] = r
        loci.extend(lr)
    undet = [g for g in GENES if not recs[g]['v3_detected']]
    if seqs is not None and undet:
        diag = orf_diagnostic(seqs, undet, hmm, work, hmmsearch, threads)
        for g in undet:
            d = diag.get(g)
            recs[g]['dna_orf_profile_hit'] = '' if d is None else (
                f"{d['orf']};model_cov={d['model_cov']:.2f};iE={d['i_evalue']:.1e}")
    elif undet and base.get('input_mode') == 'stored_proteins':
        for g in undet:
            recs[g]['dna_orf_profile_hit'] = 'not possible (no assembly sequence)'
    s2, s3 = v2_sample(v2), v3_sample(recs)
    weak_no_cds = sum(1 for r in loci if r.get('copy_class') == 'weak' and r.get('cds_identity_to_best') in ('', None))
    versions = dict(diamond=tool_version([diamond, 'version']), hmmsearch=tool_version([hmmsearch, '-h']))
    versions['prodigal'] = tool_version([prodigal, '-v']) if prodigal else 'stored proteins (earlier run)'
    summary = dict(base, module=MODULE_VERSION, proteins=len(prots), v2=s2, v3=s3,
                   weak_copies_without_cds=weak_no_cds, runtime_s=round(time.time() - t0, 1), versions=versions)
    cols = ['sample', 'gene', 'v2_state', 'v2_ambiguous_type', 'v2_n_eligible', 'v2_n_strong', 'v3_state',
            'v3_detected', 'v3_split', 'v3_reason', 'dup_counted', 'n_loci', 'n_strong_other', 'n_weak',
            'n_weak_diff', 'weak_cds_identity', 'n_other_eligible', 'profile_loci', 'best_protein', 'best_len',
            'best_partial', 'best_identity', 'best_qcov', 'best_tcov', 'best_bitscore', 'best_panel_species',
            'best_panel_assembly', 'profile_i_evalue', 'profile_model_cov', 'exp_min', 'exp_max',
            'best_aa_sha256', 'best_cds_sha256', 'dna_orf_profile_hit']
    with open(os.path.join(outdir, f'{sample}.genes.tsv'), 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter='\t', extrasaction='ignore')
        w.writeheader()
        for g in GENES:
            row = dict(recs[g])
            row['sample'] = sample
            w.writerow(row)
    lcols = ['sample', 'gene', 'source', 'locus_rank', 'copy_class', 'protein', 'contig', 'start', 'end', 'strand',
             'partial', 'aa_len', 'identity', 'qcov', 'tcov', 'bitscore', 'bitscore_ratio', 'cds_identity_to_best',
             'cds_aligned_fraction', 'contig_coverage', 'panel_species', 'profile_i_evalue', 'profile_model_cov',
             'aa_sha256', 'members']
    with open(os.path.join(outdir, f'{sample}.loci.tsv'), 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=lcols, delimiter='\t', extrasaction='ignore')
        w.writeheader()
        for r in loci:
            r = dict(r)
            r['sample'] = sample
            w.writerow(r)
    # best-copy marker sequences for overlap checks (A4) and copy-origin checks (B5)
    with open(os.path.join(outdir, f'{sample}.markers.faa'), 'w') as fa, \
            open(os.path.join(outdir, f'{sample}.markers.fna'), 'w') as fn:
        for r in loci:
            if r.get('source') == 'panel' and r.get('copy_class') in ('best', 'strong', 'weak'):
                x = prots[r['protein']]
                tag = f"{sample}|{r['gene']}|{r['copy_class']}|{r['protein']}"
                fa.write(f'>{tag}\n{x["aa"]}\n')
                if x['cds']:
                    fn.write(f'>{tag}\n{x["cds"]}\n')
    with open(os.path.join(outdir, f'{sample}.sample.json'), 'w') as fh:
        json.dump(summary, fh, indent=1)
    return summary


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--sample', required=True)
    ap.add_argument('--assembly', help='FASTA or GenBank flat file (.gz allowed)')
    ap.add_argument('--proteins', help='stored-protein mode: Prodigal protein FASTA of the original run')
    ap.add_argument('--cds', help='stored-protein mode: Prodigal CDS FASTA used for the weak-copy test')
    ap.add_argument('--n-scaffolds', type=int, help='stored-protein mode: number of sequences of the original')
    ap.add_argument('--total-length', type=int, help='stored-protein mode: total length of the original')
    ap.add_argument('--out', required=True)
    ap.add_argument('--panel-db', required=True, help='DIAMOND database of data/panel.faa')
    ap.add_argument('--hmm', required=True, help='data/ribo27.hmm')
    ap.add_argument('--ranges', required=True, help='data/expected_length_ranges.tsv')
    ap.add_argument('--threads', type=int, default=1)
    ap.add_argument('--keep', action='store_true', help='keep the work directory (Prodigal, DIAMOND outputs)')
    ap.add_argument('--no-orf-scan', action='store_true')
    ap.add_argument('--prodigal', default='prodigal')
    ap.add_argument('--diamond', default='diamond')
    ap.add_argument('--hmmsearch', default='hmmsearch')
    a = ap.parse_args()
    if bool(a.assembly) == bool(a.proteins):
        ap.error('give either --assembly or --proteins')
    if a.proteins:
        s = analyse_proteins(a.sample, a.proteins, a.out, a.panel_db, a.hmm, a.ranges, a.threads, a.keep, a.cds,
                             a.n_scaffolds, a.total_length, a.diamond, a.hmmsearch)
    else:
        s = analyse(a.sample, a.assembly, a.out, a.panel_db, a.hmm, a.ranges, a.threads, a.keep,
                    a.prodigal, a.diamond, a.hmmsearch, not a.no_orf_scan)
    print(json.dumps(dict(sample=s['sample'], v2=s['v2']['status'], v2_detected=s['v2']['detected'],
                          v3=s['v3']['status'], v3_detected=s['v3']['detected'],
                          duplicated=s['v3']['duplicated_genes'], flag=s['v3']['duplication_flag'])))


if __name__ == '__main__':
    main()
