# bio-datasets

Shared biological datasets. Data lives here (`~/bio-datasets/data/`), symlinked into project repos.

## Structure

```
data/           ← .gitignore'd, holds actual data
  rfam/         ← Rfam alignments + clan membership
  silva/raw/    ← SILVA NR99 rRNA alignments + trees
  crw/          ← CRW curated rRNA alignments + BPSEQ structures
  gtrnadb/      ← tRNA sequences (via Rfam RF00005)
  ucsc/assembly/hg38/msa/multiz100way/  ← UCSC 100-way MAFs
  pfam/         ← Pfam seed alignments (Stockholm)
  treefam/      ← TreeFam gene family alignments + trees
  balibase/     ← BAliBASE alignment benchmark
  zoonomia/     ← Zoonomia Cactus 241-mammal HAL alignment (~600 GB)
  tiberius/     ← Tiberius training genomes + annotations (37 mammals)
  annevo/       ← AnnEvo training data (genomes + Cactus alignment)

fetch/          ← parallel structure with download scripts
  common.py     ← shared utilities (idempotent download, symlink safety)
  rfam/fetch.py
  silva/fetch.py
  crw/fetch.py
  gtrnadb/fetch.py
  ucsc/assembly/hg38/msa/multiz100way/fetch.py
  pfam/fetch.py
  treefam/fetch.py
  balibase/fetch.py
  zoonomia/fetch.py
  tiberius/fetch.py
  annevo/fetch.py
```

## Usage

```bash
# Fetch a specific dataset
python fetch/rfam/fetch.py

# Fetch with custom output dir
python fetch/silva/fetch.py --outdir /scratch/silva

# Symlink into a project
ln -s ~/bio-datasets/data/rfam ~/my-project/rfam_data
```

## Contract

- `data/` is never committed (`.gitignore`'d)
- Fetch scripts are **idempotent**: skip existing files
- Fetch scripts **never delete symlinks** or existing data
- Each script defaults to `data/<dataset>/` relative to repo root

## Migration from ~/datasets/

If data already exists in `~/datasets/`, move it into `data/` here:
```bash
mv ~/datasets/pfam ~/bio-datasets/data/pfam
# Update symlinks in project repos to point to ~/bio-datasets/data/pfam
```
