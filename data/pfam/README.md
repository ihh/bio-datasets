# Pfam

Pfam is a database of protein families, each represented by a multiple sequence alignment and a hidden Markov model. Pfam data is now maintained as part of InterPro.

## Data layout

```
data/pfam/
  Pfam-A.seed.gz              # Bulk seed alignments (Stockholm format)
  fasta/{PFID}.fasta           # Per-family aligned FASTA (gaps preserved)
  trees/{PFID}.nwk             # FastTree -lg (LG08 model) Newick trees
  esm2/{PFID}.npz              # ESM2-650M per-position embeddings (float16)
```

### ESM2 embeddings format

Each `.npz` file contains one array per sequence in the family:
- **Key**: sequence name (matching FASTA header)
- **Value**: `(L_ungapped, 1280)` float16 array — ESM2 last-layer representation
- Sequences capped at 64 per family, 1022 residues per sequence
- Gap positions have zero embeddings (map ungapped ESM2 output to aligned columns)

## Fetch

```bash
# 1. Download bulk seed alignments
python fetch/pfam/fetch.py --bulk

# 2. Parse, build trees, compute ESM2 embeddings
# (run from carabs repo with pfam-split.json)
python experiments/preprocess_pfam.py --stage all
```

## Splits

Family train/val/test splits live in `splits/`:

- **`811-clan-resistant.json`**: 21,667 train / 2,430 val / 3,384 test families
  - 812 clans, 27,481 total families
  - Clan-aware splitting to prevent homology leakage
  - 59 benchmark families (held out from all splits)
  - Seed: 42, ratios: 80/10/10
