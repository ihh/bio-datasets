#!/usr/bin/env python3
"""Preprocess PFam seed alignments: parse bulk file and build trees.

Pipeline:
  1. Parse bulk Pfam-A.seed.gz → per-family aligned FASTA files
  2. Build approximately-ML trees using FastTree with the LG08 model

Usage:
  python fetch/pfam/preprocess.py --stage parse   # split bulk file into per-family FASTA
  python fetch/pfam/preprocess.py --stage trees   # build FastTree LG trees
  python fetch/pfam/preprocess.py --stage all     # run both stages
"""

import argparse
import gzip
import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import default_outdir, log


# ── Amino acid alphabet ──────────────────────────────────────────────────
AA_ALPHABET = "ACDEFGHIKLMNPQRSTVWY"
AA_TO_IDX = {aa: i + 1 for i, aa in enumerate(AA_ALPHABET)}  # 1-20
GAP_IDX = 21  # gap token
PAD_IDX = 0
VOCAB_SIZE = 22  # 0=pad, 1-20=AA, 21=gap

def tokenize_protein(seq: str) -> np.ndarray:
    """Tokenize a gapped protein sequence. Unknown residues → gap."""
    tokens = []
    for c in seq.upper():
        if c in AA_TO_IDX:
            tokens.append(AA_TO_IDX[c])
        elif c in ('-', '.'):
            tokens.append(GAP_IDX)
        else:
            tokens.append(GAP_IDX)  # X, B, Z, etc. → gap
    return np.array(tokens, dtype=np.int32)


# ── Stage 1: Parse bulk PFam ─────────────────────────────────────────────

def parse_pfam_bulk(bulk_path: Path, split_path: Path, out_dir: Path):
    """Parse Pfam-A.seed.gz and write per-family FASTA files for families in split."""
    with open(split_path) as f:
        split = json.load(f)

    wanted = set(split['train'] + split['val'] + split['test'])
    print(f"Wanted families: {len(wanted)}")

    out_dir.mkdir(parents=True, exist_ok=True)
    fasta_dir = out_dir / "fasta"
    fasta_dir.mkdir(exist_ok=True)

    found = 0
    skipped_small = 0

    with gzip.open(bulk_path, 'rt') as f:
        current_acc = None
        current_seqs = []  # list of (name, aligned_seq)
        in_block = False

        for line in f:
            line = line.rstrip('\n')

            if line.startswith('# STOCKHOLM 1.0'):
                current_acc = None
                current_seqs = []
                in_block = True
                continue

            if line.startswith('#=GF AC'):
                acc = line.split()[-1].split('.')[0]
                current_acc = acc
                continue

            if line.startswith('//'):
                # End of alignment block
                if current_acc and current_acc in wanted:
                    if len(current_seqs) >= 3:  # need ≥3 for meaningful tree
                        _write_fasta(fasta_dir / f"{current_acc}.fasta", current_seqs)
                        found += 1
                    else:
                        skipped_small += 1
                    if found % 2000 == 0 and found > 0:
                        print(f"  wrote {found} families...")
                in_block = False
                continue

            if line.startswith('#') or not in_block:
                continue

            # Sequence line: "name  aligned_sequence"
            parts = line.split()
            if len(parts) == 2:
                current_seqs.append((parts[0], parts[1]))

    print(f"Wrote {found} families ({skipped_small} skipped with <3 seqs)")
    return found


def _write_fasta(path: Path, seqs: list):
    """Write aligned sequences as FASTA (keeping gaps for FastTree)."""
    with open(path, 'w') as f:
        for name, seq in seqs:
            f.write(f">{name}\n{seq}\n")


# ── Stage 2: Build FastTree LG trees ─────────────────────────────────────

def build_trees(fasta_dir: Path, tree_dir: Path):
    """Build approximately-ML trees using FastTree with the LG08 model.

    FastTree uses a heuristic minimum-evolution starting tree refined by
    subtree-pruning-regrafting (SPR) moves under the LG08 amino acid model,
    producing trees comparable in quality to PhyML but much faster.
    """
    tree_dir.mkdir(parents=True, exist_ok=True)

    fasta_files = sorted(fasta_dir.glob("*.fasta"))
    print(f"Building trees for {len(fasta_files)} families...")

    done = 0
    errors = 0
    for fasta_path in fasta_files:
        acc = fasta_path.stem
        tree_path = tree_dir / f"{acc}.nwk"
        if tree_path.exists():
            done += 1
            continue

        try:
            result = subprocess.run(
                ["FastTree", "-lg", "-quiet", "-nosupport", "-nopr"],
                stdin=open(fasta_path),
                capture_output=True, text=True, timeout=60
            )
            if result.returncode == 0 and result.stdout.strip():
                tree_path.write_text(result.stdout.strip() + "\n")
                done += 1
            else:
                errors += 1
                if errors <= 5:
                    print(f"  FastTree error for {acc}: {result.stderr[:200]}")
        except subprocess.TimeoutExpired:
            errors += 1
            if errors <= 5:
                print(f"  FastTree timeout for {acc}")

        if done % 2000 == 0 and done > 0:
            print(f"  built {done} trees...")

    print(f"Trees: {done} built, {errors} errors")


# ── Main ──────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Preprocess PFam: parse bulk file and build trees")
    parser.add_argument("--stage", choices=["parse", "trees", "all"], required=True)
    parser.add_argument("--pfam-gz", type=Path,
                        default=default_outdir("pfam") / "Pfam-A.seed.gz")
    parser.add_argument("--split", type=Path,
                        default=Path(__file__).resolve().parent / "splits" / "811-clan-resistant.json")
    parser.add_argument("--out-dir", type=Path,
                        default=default_outdir("pfam"))
    args = parser.parse_args()

    fasta_dir = args.out_dir / "fasta"
    tree_dir = args.out_dir / "trees"

    if args.stage in ("parse", "all"):
        print("=== Stage 1: Parse bulk PFam ===")
        parse_pfam_bulk(args.pfam_gz, args.split, args.out_dir)

    if args.stage in ("trees", "all"):
        print("=== Stage 2: Build FastTree LG trees ===")
        build_trees(fasta_dir, tree_dir)


if __name__ == "__main__":
    main()
