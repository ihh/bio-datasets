#!/usr/bin/env python3
"""
Prepare Pfam seed data from bulk download.

Assumes data/pfam/Pfam-A.seed.gz already exists (downloaded by fetch.py --bulk).
Splits the bulk concatenated Stockholm file into individual .sto files in
data/pfam/seed/, then creates convenience symlinks.

Usage:
    python fetch/pfam/prepare.py
    python fetch/pfam/prepare.py --data-dir /path/to/data
"""

import argparse
import gzip
import os
import sys
from pathlib import Path


def split_bulk_stockholm(bulk_gz: Path, seed_dir: Path) -> int:
    """Split bulk Pfam-A.seed.gz into individual PFxxxxx.sto files.

    Each family in the bulk file starts with '# STOCKHOLM 1.0' and ends
    with '//'. The accession is extracted from the '#=GF AC' line.

    Returns the number of families written.
    """
    # Skip if seed_dir already has .sto files
    existing = list(seed_dir.glob("*.sto"))
    if existing:
        print(f"  {seed_dir} already has {len(existing)} .sto files, skipping split")
        return len(existing)

    if not bulk_gz.exists():
        print(f"ERROR: bulk file not found: {bulk_gz}")
        print(f"  Download it first: python fetch/pfam/fetch.py --bulk")
        sys.exit(1)

    seed_dir.mkdir(parents=True, exist_ok=True)

    count = 0
    lines = []
    accession = None

    print(f"  Splitting {bulk_gz} into {seed_dir}/")
    with gzip.open(bulk_gz, "rt", encoding="latin-1") as f:
        for line in f:
            lines.append(line)

            # Extract accession from #=GF AC line
            if line.startswith("#=GF AC"):
                # Format: #=GF AC   PF00001.23
                parts = line.strip().split()
                if len(parts) >= 3:
                    # Take accession without version: PF00001.23 -> PF00001
                    accession = parts[2].split(".")[0]

            # End of family block
            if line.strip() == "//":
                if accession:
                    outpath = seed_dir / f"{accession}.sto"
                    with open(outpath, "w") as out:
                        out.writelines(lines)
                    count += 1
                    if count % 5000 == 0:
                        print(f"    {count} families written...")
                else:
                    print(f"  WARNING: skipping block without accession "
                          f"(lines {len(lines)})")
                lines = []
                accession = None

    print(f"  Wrote {count} families to {seed_dir}/")
    return count


def create_symlinks(data_dir: Path):
    """Create convenience symlinks.

    - data/pfam/seed/splits/v1.json -> ../../../../fetch/pfam/splits/811-clan-resistant.json
    - data/pfam-seed -> pfam/seed  (backward compatibility)
    """
    # Splits symlink
    splits_dir = data_dir / "pfam" / "seed" / "splits"
    splits_dir.mkdir(parents=True, exist_ok=True)
    v1_link = splits_dir / "v1.json"

    # Relative path from data/pfam/seed/splits/ to fetch/pfam/splits/
    splits_target = Path("../../../../fetch/pfam/splits/811-clan-resistant.json")

    if v1_link.is_symlink() or v1_link.exists():
        if v1_link.is_symlink() and os.readlink(v1_link) == str(splits_target):
            print(f"  {v1_link} already points to {splits_target}")
        else:
            v1_link.unlink()
            v1_link.symlink_to(splits_target)
            print(f"  Created {v1_link} -> {splits_target}")
    else:
        v1_link.symlink_to(splits_target)
        print(f"  Created {v1_link} -> {splits_target}")

    # Backward-compat symlink: data/pfam-seed -> pfam/seed
    compat_link = data_dir / "pfam-seed"
    compat_target = Path("pfam/seed")

    if compat_link.is_symlink():
        if os.readlink(compat_link) == str(compat_target):
            print(f"  {compat_link} already points to {compat_target}")
        else:
            compat_link.unlink()
            compat_link.symlink_to(compat_target)
            print(f"  Created {compat_link} -> {compat_target}")
    elif compat_link.exists():
        print(f"  WARNING: {compat_link} exists and is not a symlink, skipping")
    else:
        compat_link.symlink_to(compat_target)
        print(f"  Created {compat_link} -> {compat_target}")


def main():
    parser = argparse.ArgumentParser(
        description="Prepare Pfam seed data from bulk download")
    parser.add_argument("--data-dir", type=Path, default=None,
                        help="Path to data/ directory (default: auto-detect)")
    args = parser.parse_args()

    # Auto-detect data directory
    if args.data_dir:
        data_dir = args.data_dir
    else:
        # Assume running from repo root
        data_dir = Path(__file__).resolve().parent.parent.parent / "data"

    if not data_dir.is_dir():
        print(f"ERROR: data directory not found: {data_dir}")
        sys.exit(1)

    print(f"Data directory: {data_dir}")

    # Step 1: Split bulk file into individual .sto files
    bulk_gz = data_dir / "pfam" / "Pfam-A.seed.gz"
    # Also check the seed subdirectory (where we may have moved it)
    if not bulk_gz.exists():
        alt = data_dir / "pfam" / "seed" / "Pfam-A.seed.gz"
        if alt.exists():
            bulk_gz = alt

    seed_dir = data_dir / "pfam" / "seed"
    print("\n[1/2] Splitting bulk Stockholm file...")
    split_bulk_stockholm(bulk_gz, seed_dir)

    # Step 2: Create symlinks
    print("\n[2/2] Creating symlinks...")
    create_symlinks(data_dir)

    print("\nDone.")


if __name__ == "__main__":
    main()
