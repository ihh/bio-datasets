# Fetch: AnnEvo training/test/validation sets

## Task

Download the training, test, and validation datasets used by
[AnnEvo](https://github.com/Gaius-Augustus/AnnEvo), a genome annotation
tool that uses evolutionary sequence conservation (multiple alignments)
to improve gene prediction accuracy.

## What to fetch

AnnEvo extends Tiberius with cross-species alignment features. Its
training data includes both genome sequences and multiple sequence
alignments (typically Cactus or multiz alignments).

1. **Check the AnnEvo repository** at https://github.com/Gaius-Augustus/AnnEvo
   for data download instructions.

2. **Look for**:
   - Species-level genome FASTA + annotation (GTF/GFF3) for each split
   - Multiple alignment files (HAL, MAF, or pre-extracted alignment features)
   - The species phylogeny / tree used for alignment
   - Any Cactus alignment outputs or pre-computed conservation scores
   - Preprocessing/feature extraction scripts from the authors

3. **Download to**:
   - `data/annevo/training/` — training species genomes + alignments + annotations
   - `data/annevo/test/` — held-out test species
   - `data/annevo/validation/` — validation species

4. **Key references**:
   - Check the AnnEvo paper and supplementary materials for exact
     species lists, genome assembly versions, and alignment parameters
   - Training data may be hosted on the Hiller/Stanke lab servers,
     Zenodo, or linked from the GitHub README
   - Overlap with Tiberius training data is likely — reuse shared
     genomes/annotations where possible (symlink within data/)

## Idempotency

Follow the conventions in `fetch/common.py`:
- Skip files that already exist
- Never delete symlinks
- Default --outdir to `data/annevo/`
