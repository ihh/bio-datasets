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

4. **Override convention**: The bio-datasets location (`~/bio-datasets` by default) is
   overridable via `BIO_DATASETS_HOME` env var or `--bio-datasets DIR` CLI flag.
   Every project script that uses bio-datasets must support both overrides.

## Migrating other repos' data pipelines

When a project repo has its own data-fetching code (download scripts, crawlers,
`--fetch` modes), migrate it to prefer bio-datasets:

### Process (confirm with user before starting each phase)

**Phase 1: Plan** — Identify the repo's data-fetching entry points and where data
currently lives. Confirm the plan with the user before writing any code.

**Phase 2: Execute** — After user approval:

1. Create `fetch/<dataset>/fetch.py` here if it doesn't exist
2. Add a `bio_datasets.py` utility to the project repo (or copy the pattern from
   `tkf-mixdom/python/tkfmixdom/jax/util/bio_datasets.py`) that:
   - Auto-detects `~/bio-datasets` (or `$BIO_DATASETS_HOME`)
   - Falls back to local directory when bio-datasets is absent
   - Creates symlinks from the project's local data dir to bio-datasets
   - Supports `--bio-datasets DIR` CLI override
3. Wire the project's fetch/download scripts to use `resolve_data_dir()`
4. Wire the project's training/analysis scripts to use `resolve_data_dir()`
5. Test: verify scripts work both with and without bio-datasets present

### Reference implementation

See `tkf-mixdom` for the canonical example:
- `python/tkfmixdom/jax/util/bio_datasets.py` — utility module
- `train_pfam.py` — uses `resolve_data_dir("pfam", "pfam/")` + `ensure_symlinks()`
- `maraschino.py fetch` — uses `resolve_data_dir("pfam", args.out_dir)`

Key API:
```python
from tkfmixdom.jax.util.bio_datasets import resolve_data_dir, ensure_symlinks

# Auto-detects ~/bio-datasets, falls back to local
data_dir = resolve_data_dir("pfam", local_fallback="pfam/")

# For mixed dirs (data + checkpoints), create per-file symlinks
ensure_symlinks(data_dir, Path("pfam/"), pattern="*.sto")
```

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
| zoonomia | ~600 GB | fetch/zoonomia/fetch.py | Cactus 241-mammal HAL alignment |
| tiberius | ~100-150 GB | fetch/tiberius/fetch.py | 37 mammalian genomes + RefSeq annotations |
| annevo | ~100-150 GB+ | fetch/annevo/fetch.py | AnnEvo training (genomes + Cactus alignment) |
