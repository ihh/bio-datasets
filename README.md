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
  pfam/         ← (placeholder)

fetch/          ← parallel structure with download scripts
  common.py     ← shared utilities (idempotent download, symlink safety)
  rfam/fetch.py
  silva/fetch.py
  crw/fetch.py
  gtrnadb/fetch.py
  ucsc/assembly/hg38/msa/multiz100way/fetch.py
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
