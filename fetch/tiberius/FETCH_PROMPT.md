# Fetch: Tiberius training/test/validation sets

## Task

Download the training, test, and validation datasets used by
[Tiberius](https://github.com/Gaius-Augustus/Tiberius), a deep learning
gene finder for eukaryotic genomes.

## What to fetch

Tiberius training data consists of genome sequences (FASTA) paired with
gene annotations (GTF/GFF3) for supervised training of a sequence-to-annotation
model. The datasets are split into training/test/validation by species.

1. **Check the Tiberius repository** at https://github.com/Gaius-Augustus/Tiberius
   for data download instructions, typically in README.md or a `data/` directory.

2. **Look for**:
   - Pre-built training sets (often hosted on a university server or Zenodo)
   - The species list for each split (train/test/val)
   - Genome FASTA files (softmasked, per-species)
   - Reference annotations (GTF/GFF3, e.g., from Ensembl/NCBI)
   - Any preprocessing scripts provided by the Tiberius authors

3. **Download to**:
   - `data/tiberius/training/` — training species genomes + annotations
   - `data/tiberius/test/` — held-out test species
   - `data/tiberius/validation/` — validation species

4. **Key references**:
   - Gabriel et al. (2024) "Tiberius: End-to-End Deep Learning with an HMM
     for Gene Prediction" — check supplementary for exact species/versions
   - AUGUSTUS training data may overlap: https://github.com/Gaius-Augustus/Augustus

## Idempotency

Follow the conventions in `fetch/common.py`:
- Skip files that already exist
- Never delete symlinks
- Default --outdir to `data/tiberius/`
