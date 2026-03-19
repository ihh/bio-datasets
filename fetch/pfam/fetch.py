#!/usr/bin/env python3
"""
Fetch Pfam seed alignments (Stockholm format).

Downloads seed alignments from Pfam (InterPro). Supports fetching
a random subset or specific families.

Usage:
    python fetch/pfam/fetch.py                          # fetch all (~20k families)
    python fetch/pfam/fetch.py --random 100 --seed 42   # random 100 families
    python fetch/pfam/fetch.py --families PF00001,PF00002
"""

import argparse
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import default_outdir, safe_outdir, safe_download, log

# Pfam is now part of InterPro; seed alignments available via EBI FTP
BASE_URL = "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/"
SEED_URL = BASE_URL + "Pfam-A.seed.gz"


def fetch_individual(family: str, outdir: Path) -> bool:
    """Fetch a single family's seed alignment via InterPro API."""
    url = f"https://www.ebi.ac.uk/interpro/api/entry/pfam/{family}/?annotation=alignment:seed&download"
    dest = outdir / f"{family}.sto"
    return safe_download(url, dest, min_size=100)


def main():
    parser = argparse.ArgumentParser(description="Fetch Pfam seed alignments")
    parser.add_argument("--outdir", type=Path, default=default_outdir("pfam"))
    parser.add_argument("--families", type=str, default=None,
                        help="Comma-separated Pfam accessions")
    parser.add_argument("--random", type=int, default=0,
                        help="Fetch N random families (0 = all or --families)")
    parser.add_argument("--seed", type=int, default=42,
                        help="Random seed for --random selection")
    parser.add_argument("--bulk", action="store_true",
                        help="Download bulk Pfam-A.seed.gz instead of individual families")
    args = parser.parse_args()

    outdir = safe_outdir(args.outdir)

    if args.bulk:
        log.info("Downloading bulk Pfam-A.seed.gz")
        safe_download(SEED_URL, outdir / "Pfam-A.seed.gz", min_size=1000)
        return

    if args.families:
        families = [f.strip() for f in args.families.split(",")]
    elif args.random > 0:
        # Use maraschino.py's fetch logic if available, else use a fixed list
        import random
        rng = random.Random(args.seed)
        # Pfam families are PF00001-PF20000+ ; sample from known range
        all_fams = [f"PF{i:05d}" for i in range(1, 20001)]
        families = rng.sample(all_fams, min(args.random, len(all_fams)))
    else:
        log.error("Specify --families, --random N, or --bulk")
        sys.exit(1)

    log.info("Fetching %d Pfam families to %s", len(families), outdir)
    ok = 0
    for fam in families:
        if fetch_individual(fam, outdir):
            ok += 1
    log.info("Pfam: %d/%d families fetched", ok, len(families))


if __name__ == "__main__":
    main()
