# DryQC

DryQC is a quality-control tool for bacterial whole-genome sequencing data. It checks whether 27 ribosomal protein
genes, which generally occur once in every bacterial genome, are recovered in paired-end Illumina reads and in
genome assemblies, and it flags assemblies in which nearly all of these genes are duplicated.

DryQC reports gene-level evidence (which genes were detected, split, duplicated or not detected, and why). It does
not estimate genome completeness, identify species or diagnose contamination.

## Components

| Script | Input | What it does |
|---|---|---|
| `scripts/dryqc_reads.py` | Paired-end FASTQ | Read check: reads are aligned to the reference panel with DIAMOND blastx; a gene is detected when at least five reads support it |
| `scripts/dryqc_assembly.py` | Assembly (FASTA or GenBank) | Assembly check: genes are predicted with Prodigal and assigned with DIAMOND blastp and profile hidden Markov models; copies are counted per locus |

A sample passes when at least 24 of the 27 genes are detected. An assembly is flagged for extensive duplication when
at least 25 genes are duplicated. The two outcomes are independent.

## Requirements

- Python 3.8 or later (3.11 used), Biopython
- DIAMOND, Prodigal 2.6.3, HMMER 3.4

Versions used in the study: DIAMOND v0.9.30.131 for the read check; DIAMOND v2.1.8, Prodigal v2.6.3 and HMMER 3.4
for the assembly check. A DIAMOND database must be built with the same DIAMOND version that searches it.

```bash
conda env create -f environment.yml          # assembly check and tests (DIAMOND 2.1.8)
conda env create -f environment-reads.yml    # read check (DIAMOND 0.9.30)
```

## Usage

Read check (Illumina paired-end reads):

```bash
conda activate dryqc-reads
python scripts/dryqc_reads.py --r1 sample_R1.fastq.gz --r2 sample_R2.fastq.gz \
    --panel data/panel.faa --outdir out_reads --sample sample --platform illumina --threads 8
```

The DIAMOND database of the panel is built next to `data/panel.faa` on the first run.
Outputs: `<sample>.ribo27.fastq.summary.csv` (sample result), `<sample>.ribo27.fastq.detail.csv` (one row per gene)
and `<sample>.ribo27.reads.m8` (alignments).

Assembly check:

```bash
conda activate dryqc
diamond makedb --in data/panel.faa -d panel
python scripts/dryqc_assembly.py --sample sample --assembly sample.fasta --out out_assembly \
    --panel-db panel.dmnd --hmm data/dryqc_profiles.hmm --ranges data/expected_length_ranges.tsv --threads 8
```

Outputs: `<sample>.genes.tsv` (one row per gene), `<sample>.loci.tsv` (one row per candidate locus) and
`<sample>.sample.json` (sample result and software versions). Each file reports the final rules (`v3_*` columns)
and, for comparison, the earlier rules (`v2_*` columns).

Only Illumina paired-end reads were evaluated in the study. The `--platform ont`, `pacbio_hifi` and `pacbio_clr`
settings of the read check have not been evaluated.

## Reference data

| File | Content | SHA-256 |
|---|---|---|
| `data/panel.faa` | 3,417 protein records of the 27 genes from 127 RefSeq assemblies | `745abb9fdbbe55785f57cc24307226ed365377c61532308b5ddd159ea2418ba1` |
| `data/panel_records.tsv` | Gene, source assembly and protein accession of each panel record | `192e88c8226302071d3f9b32e4bbd2a73a117ffef95b894daa0ca91fd85d5e90` |
| `data/dryqc_profiles.hmm` | Profile hidden Markov models of the 27 genes (HMMER 3.4) | `e1552b58e08d086306ecc278b41c097f444c6ee9be2ee0fdaeea441e6e200403` |
| `data/expected_length_ranges.tsv` | Expected protein length range of each gene from type strains | `38db98b749830f987b223ae4dcc49c613f063111aef8709032e58cfbfb15e506` |

## Tests

```bash
conda activate dryqc
python tests/run_tests.py test_work 4
```

Synthetic assemblies with known deletions, splits, copies, truncations, long forms and a second genome are checked;
all 11 checks should pass.

## Citation

Software (version 1.0.0):

On Y, Yoon EJ. DryQC (v1.0.0). Zenodo. 2026. https://doi.org/10.5281/zenodo.24452626

Manuscript:

On Y, Yoon EJ. DryQC: a biologically informed tool for quality assessment of bacterial whole-genome sequencing data.
Manuscript submitted.

## License

MIT License. See `LICENSE`.
