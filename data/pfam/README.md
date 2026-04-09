# Pfam

Pfam is a database of protein families, each represented by a multiple sequence alignment and a hidden Markov model. Pfam data is now maintained as part of InterPro.

## Data layout

```
data/pfam/
  Pfam-A.seed.gz              # Bulk seed alignments (Stockholm format)
  fasta/{PFID}.fasta           # Per-family aligned FASTA (gaps preserved)
  trees/{PFID}.nwk             # FastTree -lg approximately-ML trees (LG08 model)
  esm2/{PFID}.npz              # ESM2-650M per-position embeddings (float16)
```

### Trees

Trees are built with FastTree using the LG08 amino acid substitution model
(`FastTree -lg`). FastTree uses a heuristic minimum-evolution starting tree
refined by subtree-pruning-regrafting (SPR) moves, producing trees comparable
in quality to PhyML but much faster. Branch lengths are in expected
substitutions per site.

### ESM2 embeddings format

Each `.npz` file contains one array per sequence in the family:
- **Key**: sequence name (matching FASTA header)
- **Value**: `(L_ungapped, 1280)` float16 array — ESM2 last-layer representation
- Sequences capped at 64 per family, 1022 residues per sequence
- Gap positions have zero embeddings (map ungapped ESM2 output to aligned columns)

## Fetch and preprocess

```bash
# 1. Download bulk seed alignments
python fetch/pfam/fetch.py --bulk

# 2. Parse per-family FASTA and build trees
python fetch/pfam/preprocess.py --stage all

# 3. (Optional) Precompute ESM2 embeddings — requires torch + esm
#    This is model-specific; see downstream project for the script.
```

## Splits

Family train/val/test splits live in `fetch/pfam/splits/`:

- **`811-clan-resistant.json`**: 21,667 train / 2,430 val / 3,384 test families
  - 812 clans, 27,481 total families
  - Clan-aware splitting to prevent homology leakage
  - 59 benchmark families (held out from all splits)
  - Seed: 42, ratios: 80/10/10
