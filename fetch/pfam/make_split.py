#!/usr/bin/env python3
"""Generate a canonical clan-aware train/val/test split for Pfam.

Key properties:
  - No clan spans multiple splits (leak-safe)
  - Benchmark families (BAliBASE etc.) forced into test split
  - Stratified by clan size so splits are balanced in family count
  - Deterministic from seed
  - Saved as JSON for reuse by all scripts

Usage:
    python fetch/pfam/make_split.py --pfam-seed data/pfam-seed/ --out data/pfam-seed/splits/v1.json
    python fetch/pfam/make_split.py --pfam-seed data/pfam-seed/ --out data/pfam-seed/splits/v1.json --benchmark-families BAliBASE.txt
"""

import argparse
import gzip
import hashlib
import json
import os
import sys
from collections import defaultdict
from pathlib import Path


def load_clan_membership(clan_file):
    """Load Pfam-A.clans.tsv.gz → {family: clan_id}. Singletons get pseudo-clan."""
    clans = {}
    opener = gzip.open if str(clan_file).endswith('.gz') else open
    with opener(clan_file, 'rt') as f:
        for line in f:
            parts = line.strip().split('\t')
            fam = parts[0]
            clan = parts[1] if len(parts) >= 2 and parts[1] else None
            clans[fam] = clan
    return clans


def make_split(families, clan_membership, benchmark_families=None,
               ratios=(0.8, 0.1, 0.1), seed=42):
    """Create a clan-aware split.

    Algorithm:
    1. Force benchmark families' clans into the test split
    2. Group remaining families by clan (singletons = individual groups)
    3. Sort clans by size (largest first) for balanced packing
    4. Assign clans to train/val/test to approximate target ratios
       using a greedy algorithm that assigns each clan to the split
       that is furthest below its target count

    Returns dict with keys: train, val, test, metadata
    """
    benchmark_families = set(benchmark_families or [])
    train_r, val_r, test_r = ratios
    n_total = len(families)

    # Group families by clan
    clan_to_fams = defaultdict(list)
    for fam in families:
        clan = clan_membership.get(fam)
        if clan:
            clan_to_fams[clan].append(fam)
        else:
            # Singleton: its own group
            clan_to_fams[f'_singleton_{fam}'].append(fam)

    # Identify clans that must go to test (contain benchmark families)
    forced_test_clans = set()
    for fam in benchmark_families:
        clan = clan_membership.get(fam)
        if clan:
            forced_test_clans.add(clan)
        elif fam in families:
            forced_test_clans.add(f'_singleton_{fam}')

    # Separate forced-test clans from available clans
    forced_test_fams = []
    available_clans = []
    for clan_id, fams in clan_to_fams.items():
        if clan_id in forced_test_clans:
            forced_test_fams.extend(fams)
        else:
            available_clans.append((clan_id, fams))

    # Sort available clans by size (largest first) for balanced packing
    # Use hash for deterministic tiebreaking
    def sort_key(item):
        clan_id, fams = item
        h = hashlib.sha256(f"{seed}:{clan_id}".encode()).hexdigest()
        return (-len(fams), h)
    available_clans.sort(key=sort_key)

    # Target counts (excluding forced-test families)
    n_available = n_total - len(forced_test_fams)
    # Adjust ratios: test already has forced families, so reduce test target
    target_test_from_available = max(0, int(n_total * test_r) - len(forced_test_fams))
    target_train = int(n_available * train_r / (train_r + val_r + test_r * (target_test_from_available / max(n_available * test_r, 1))))
    target_val = int(n_available * val_r / (train_r + val_r))

    # Greedy assignment: assign each clan to the split furthest below target
    splits = {'train': [], 'val': [], 'test': list(forced_test_fams)}
    counts = {'train': 0, 'val': 0, 'test': len(forced_test_fams)}
    targets = {
        'train': int(n_total * train_r),
        'val': int(n_total * val_r),
        'test': int(n_total * test_r),
    }

    for clan_id, fams in available_clans:
        # Which split is furthest below its target?
        deficits = {s: targets[s] - counts[s] for s in ['train', 'val', 'test']}
        best = max(deficits, key=deficits.get)
        splits[best].extend(fams)
        counts[best] += len(fams)

    # Sort each split for reproducibility
    for s in splits:
        splits[s].sort()

    return splits, counts, targets


def main():
    parser = argparse.ArgumentParser(description='Generate Pfam train/val/test split')
    parser.add_argument('--pfam-seed', type=str, required=True,
                        help='Path to pfam-seed directory with .sto files')
    parser.add_argument('--out', type=str, required=True,
                        help='Output JSON file for the split')
    parser.add_argument('--clan-file', type=str, default=None,
                        help='Pfam-A.clans.tsv.gz (default: auto-detect in pfam-seed/)')
    parser.add_argument('--benchmark-families', type=str, default=None,
                        help='File with benchmark Pfam accessions (one per line)')
    parser.add_argument('--benchmark-dirs', type=str, nargs='*', default=None,
                        help='Benchmark directories to scan for PF* families')
    parser.add_argument('--ratios', type=str, default='0.8,0.1,0.1',
                        help='Train/val/test ratios (default: 0.8,0.1,0.1)')
    parser.add_argument('--seed', type=int, default=42)
    args = parser.parse_args()

    pfam_dir = Path(args.pfam_seed)

    # Find all families with MSAs
    families = sorted(
        f.stem for f in pfam_dir.glob('*.sto'))
    if not families:
        families = sorted(
            f.stem.replace('.sto', '') for f in pfam_dir.glob('*.sto.gz'))
    print(f"Found {len(families)} families in {pfam_dir}")

    # Load clan membership
    clan_file = args.clan_file or str(pfam_dir / 'Pfam-A.clans.tsv.gz')
    if not os.path.exists(clan_file):
        print(f"ERROR: clan file not found: {clan_file}")
        sys.exit(1)
    clans = load_clan_membership(clan_file)
    print(f"Loaded clan membership: {len(clans)} families, "
          f"{len(set(c for c in clans.values() if c))} clans")

    # Collect benchmark families
    bench_fams = set()
    if args.benchmark_families:
        with open(args.benchmark_families) as f:
            bench_fams = set(l.strip() for l in f if l.strip().startswith('PF'))
    if args.benchmark_dirs:
        import glob as glob_mod
        for d in args.benchmark_dirs:
            for pattern in ['**/PF*.???', '**/PF*.????']:
                for p in glob_mod.glob(os.path.join(d, pattern), recursive=True):
                    fam = os.path.basename(p).split('.')[0]
                    if fam.startswith('PF'):
                        bench_fams.add(fam)
    # Only keep benchmark families that we actually have
    bench_fams = bench_fams & set(families)
    print(f"Benchmark families: {len(bench_fams)} (forced to test split)")

    # Make split
    ratios = tuple(float(r) for r in args.ratios.split(','))
    splits, counts, targets = make_split(
        families, clans, bench_fams, ratios, args.seed)

    print(f"\nSplit results:")
    for s in ['train', 'val', 'test']:
        print(f"  {s}: {counts[s]} families (target: {targets[s]})")

    # Verify no clan leaks
    split_of = {}
    for s, fams in splits.items():
        for fam in fams:
            clan = clans.get(fam)
            if clan:
                if clan in split_of and split_of[clan] != s:
                    print(f"  ERROR: clan {clan} in both {split_of[clan]} and {s}!")
                split_of[clan] = s
    print(f"  Clan leak check: PASS ({len(split_of)} clans verified)")

    # Save
    output = {
        'version': 1,
        'seed': args.seed,
        'ratios': list(ratios),
        'n_families': len(families),
        'n_clans': len(set(c for c in clans.values() if c)),
        'benchmark_families': sorted(bench_fams),
        'counts': counts,
        'train': splits['train'],
        'val': splits['val'],
        'test': splits['test'],
    }
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open(args.out, 'w') as f:
        json.dump(output, f, indent=2)
    print(f"\nSaved: {args.out}")


if __name__ == '__main__':
    main()
