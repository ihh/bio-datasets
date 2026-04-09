#!/usr/bin/env python3
"""Generate a train/val/test split for PANTHER families.

PANTHER doesn't have clan structure like Pfam, so we use a simple
deterministic random 8:1:1 split seeded for reproducibility.

Usage:
    python fetch/panther/make_split.py
    python fetch/panther/make_split.py --msa-dir data/panther/PANTHER19.0_data/Tree_MSF
"""

import argparse
import hashlib
import json
import os
import sys
from pathlib import Path


def make_split(families, ratios=(0.8, 0.1, 0.1), seed=42):
    """Create a deterministic random split.

    Uses hash-based assignment for stability: adding/removing families
    doesn't change the split of other families.
    """
    train_r, val_r, test_r = ratios

    # Hash-based deterministic assignment
    splits = {'train': [], 'val': [], 'test': []}
    for fam in families:
        h = hashlib.sha256(f"{seed}:{fam}".encode()).hexdigest()
        v = int(h[:8], 16) / 0xFFFFFFFF  # uniform [0, 1]
        if v < train_r:
            splits['train'].append(fam)
        elif v < train_r + val_r:
            splits['val'].append(fam)
        else:
            splits['test'].append(fam)

    for s in splits:
        splits[s].sort()

    return splits


def main():
    parser = argparse.ArgumentParser(
        description='Generate PANTHER train/val/test split')
    parser.add_argument('--msa-dir', type=str, default=None,
                        help='Directory with *.AN.fasta files')
    parser.add_argument('--out', type=str, default=None,
                        help='Output JSON file (default: fetch/panther/splits/811.json)')
    parser.add_argument('--ratios', type=str, default='0.8,0.1,0.1',
                        help='Train/val/test ratios (default: 0.8,0.1,0.1)')
    parser.add_argument('--seed', type=int, default=42)
    args = parser.parse_args()

    # Find MSA directory
    if args.msa_dir:
        msa_dir = Path(args.msa_dir)
    else:
        data_dir = Path(__file__).resolve().parent.parent.parent / "data"
        msa_dir = data_dir / "panther" / "PANTHER19.0_data" / "Tree_MSF"

    if not msa_dir.is_dir():
        print(f"ERROR: MSA directory not found: {msa_dir}")
        print("Run fetch.py and prepare.py first.")
        sys.exit(1)

    # Find all families
    families = sorted(
        f.stem.replace('.AN', '')
        for f in msa_dir.glob('*.AN.fasta')
    )
    print(f"Found {len(families)} families in {msa_dir}")

    # Make split
    ratios = tuple(float(r) for r in args.ratios.split(','))
    splits = make_split(families, ratios, args.seed)

    counts = {s: len(splits[s]) for s in ['train', 'val', 'test']}
    print(f"\nSplit results (seed={args.seed}, ratios={ratios}):")
    for s in ['train', 'val', 'test']:
        print(f"  {s}: {counts[s]} families "
              f"({100 * counts[s] / len(families):.1f}%)")

    # Output path
    if args.out:
        out_path = args.out
    else:
        out_dir = Path(__file__).resolve().parent / "splits"
        out_dir.mkdir(parents=True, exist_ok=True)
        out_path = str(out_dir / "811.json")

    # Save
    output = {
        'version': 1,
        'seed': args.seed,
        'ratios': list(ratios),
        'n_families': len(families),
        'counts': counts,
        'train': splits['train'],
        'val': splits['val'],
        'test': splits['test'],
    }
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, 'w') as f:
        json.dump(output, f, indent=2)
    print(f"\nSaved: {out_path}")


if __name__ == '__main__':
    main()
