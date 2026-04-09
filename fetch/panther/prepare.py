#!/usr/bin/env python3
"""
Prepare PANTHER 19.0 data from TreeGrafter download.

Assumes data/panther/PANTHER19.0_data.tar.gz has been downloaded by fetch.py.
Unpacks if needed and reports statistics about the MSAs.

Each PANTHER family has:
  - PTHR#####.AN.fasta    — aligned ancestral sequences (FASTA with gaps)
  - PTHR#####.newick       — phylogenetic tree (Newick, AN labels)
  - PTHR#####.bifurcate.newick — bifurcating version of tree
  - PTHR#####.tree         — NHX-annotated tree with species/event info

Usage:
    python fetch/panther/prepare.py
    python fetch/panther/prepare.py --data-dir /path/to/data
"""

import argparse
import os
import sys
import tarfile
from collections import Counter
from pathlib import Path


AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")


def parse_fasta_counts(filepath):
    """Parse aligned FASTA, return (n_seqs, alignment_length, n_valid_aa).

    Counts only standard amino acids (not gaps or ambiguous chars).
    """
    n_seqs = 0
    aln_len = 0
    total_aa = 0
    with open(filepath) as f:
        seq = ""
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if seq:
                    aln_len = max(aln_len, len(seq))
                    total_aa += sum(1 for c in seq if c in AMINO_ACIDS)
                    seq = ""
                n_seqs += 1
            else:
                seq += line
        if seq:
            aln_len = max(aln_len, len(seq))
            total_aa += sum(1 for c in seq if c in AMINO_ACIDS)
    return n_seqs, aln_len, total_aa


def unpack_tarball(data_dir):
    """Unpack PANTHER19.0_data.tar.gz if not already done."""
    tarball = data_dir / "panther" / "PANTHER19.0_data.tar.gz"
    tree_msf = data_dir / "panther" / "PANTHER19.0_data" / "Tree_MSF"

    if tree_msf.is_dir() and any(tree_msf.glob("*.AN.fasta")):
        n = len(list(tree_msf.glob("*.AN.fasta")))
        print(f"  Already unpacked: {tree_msf} ({n} families)")
        return tree_msf

    if not tarball.exists():
        print(f"ERROR: tarball not found: {tarball}")
        print(f"  Download it first: python fetch/panther/fetch.py")
        sys.exit(1)

    print(f"  Unpacking {tarball}...")
    with tarfile.open(tarball, "r:gz") as tar:
        tar.extractall(path=data_dir / "panther")
    print(f"  Unpacked to {tree_msf}")
    return tree_msf


def report_stats(tree_msf):
    """Report statistics about PANTHER families."""
    fasta_files = sorted(tree_msf.glob("*.AN.fasta"))
    n_families = len(fasta_files)
    print(f"\n  Total families: {n_families}")

    seq_counts = []
    aln_lengths = []
    aa_counts = []
    small_families = 0

    for i, f in enumerate(fasta_files):
        n_seqs, aln_len, n_aa = parse_fasta_counts(f)
        seq_counts.append(n_seqs)
        aln_lengths.append(aln_len)
        aa_counts.append(n_aa)
        if n_seqs < 2:
            small_families += 1
        if (i + 1) % 5000 == 0:
            print(f"    Scanned {i + 1}/{n_families}...")

    import numpy as np
    seq_arr = np.array(seq_counts)
    aln_arr = np.array(aln_lengths)

    print(f"\n  Families with <2 sequences: {small_families}")
    print(f"  Families with >=2 sequences: {n_families - small_families}")
    print(f"\n  Sequences per family:")
    print(f"    min={seq_arr.min()}, median={int(np.median(seq_arr))}, "
          f"mean={seq_arr.mean():.1f}, max={seq_arr.max()}")
    for pct in [10, 25, 50, 75, 90, 95, 99]:
        print(f"    p{pct}: {int(np.percentile(seq_arr, pct))}")

    print(f"\n  Alignment length:")
    print(f"    min={aln_arr.min()}, median={int(np.median(aln_arr))}, "
          f"mean={aln_arr.mean():.1f}, max={aln_arr.max()}")

    return n_families, seq_counts, aln_lengths


def main():
    parser = argparse.ArgumentParser(
        description="Prepare PANTHER data from TreeGrafter download")
    parser.add_argument("--data-dir", type=Path, default=None,
                        help="Path to data/ directory (default: auto-detect)")
    parser.add_argument("--stats-only", action="store_true",
                        help="Only report stats, skip unpacking")
    args = parser.parse_args()

    if args.data_dir:
        data_dir = args.data_dir
    else:
        data_dir = Path(__file__).resolve().parent.parent.parent / "data"

    if not data_dir.is_dir():
        print(f"ERROR: data directory not found: {data_dir}")
        sys.exit(1)

    print(f"Data directory: {data_dir}")

    # Step 1: Unpack
    if not args.stats_only:
        print("\n[1/2] Unpacking tarball...")
        tree_msf = unpack_tarball(data_dir)
    else:
        tree_msf = data_dir / "panther" / "PANTHER19.0_data" / "Tree_MSF"

    # Step 2: Report stats
    print("\n[2/2] Scanning families...")
    report_stats(tree_msf)

    print("\nDone.")


if __name__ == "__main__":
    main()
