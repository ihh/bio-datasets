# Agent Orientation: bio-datasets

## What this repo is

Centralized dataset management for biological data used across multiple projects.
Real data lives in `data/` (gitignored). Fetch scripts live in `fetch/` (committed).

## Key design patterns

1. **`fetch/` mirrors `data/`**: every `data/foo/` has a `fetch/foo/fetch.py` that can
   recreate it from scratch. The fetch scripts are the source of truth for how data
   was obtained.

2. **Idempotent, symlink-safe**: `common.py` provides `safe_download()` and `safe_outdir()`
   which NEVER delete symlinks or overwrite existing data. This is critical because
   project repos symlink into `data/`.

3. **Symlink convention**: Project repos symlink to `~/bio-datasets/data/<dataset>`,
   not to `~/datasets/` (legacy) or absolute paths.

## Adding a new dataset

1. Create `fetch/<dataset>/fetch.py` following existing examples (rfam, pfam, etc.)
2. Use `common.py` utilities for downloads
3. Test: run fetch, confirm it downloads; kill it; symlink in existing data; run again,
   confirm it skips the symlinks (does NOT delete them)
4. Commit the fetch script (never commit data)
5. Move actual data into `data/<dataset>/`
6. Update symlinks in project repos

## Migration procedure (existing data on a machine)

When data already exists somewhere (e.g. `~/datasets/pfam/` or inside a project repo):

```bash
# 1. Create fetch script if it doesn't exist
# 2. Push fetch script to github
# 3. Test fetch script starts correctly, then kill it
# 4. Symlink existing data into data/ to test symlink safety
ln -s ~/datasets/pfam ~/bio-datasets/data/pfam
python fetch/pfam/fetch.py  # should print SKIP (symlink) for all files
# 5. Move real data in, update symlinks from project repos
mv ~/datasets/pfam/* ~/bio-datasets/data/pfam/
# 6. Update project repo symlinks
ln -sf ~/bio-datasets/data/pfam ~/project/pfam_data
```

## Current datasets

| Dataset | Size | Fetch script | Notes |
|---------|------|-------------|-------|
| rfam | ~50 GB | fetch/rfam/fetch.py | Full Rfam alignments + CMs |
| silva | ~10 GB | fetch/silva/fetch.py | SILVA NR99 rRNA |
| crw | ~200 MB | fetch/crw/fetch.py | CRW curated rRNA structures |
| gtrnadb | ~50 MB | fetch/gtrnadb/fetch.py | tRNA sequences |
| ucsc | ~100 GB | fetch/ucsc/ | UCSC multiz100way MAFs |
| pfam | ~4 MB-27K files | fetch/pfam/fetch.py | Pfam seed alignments |
| treefam | ~16 GB | fetch/treefam/fetch.py | TreeFam gene families |
| balibase | ~124 MB | fetch/balibase/fetch.py | BAliBASE alignment benchmark |
