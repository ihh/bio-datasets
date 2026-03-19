#!/usr/bin/env python3
"""
Fetch Zoonomia Cactus 241-species whole-genome alignment (HAL format).

This is the progressive Cactus alignment of 241 placental mammal genomes
produced by the Zoonomia Project. It is used by ClaMSA/Tiberius de novo mode
for cross-species conservation features.

Source: https://cglgenomics.ucsc.edu/data/cactus/
Paper: Zoonomia Consortium (2020), Nature 587, 240-245

Expected files:
  - 241-mammalian-2020v2.hal  (~600 GB, the full alignment in HAL format)
  - 241-mammalian-2020v2b.hal (~600 GB, updated version if available)

The HAL file can be queried with hal2maf, halStats, etc. from the
HAL tools suite (https://github.com/ComparativeGenomicsToolkit/hal).

WARNING: This is a very large download (~600 GB). Make sure you have
sufficient disk space before running.

Usage:
    python fetch/zoonomia/fetch.py [--outdir data/zoonomia] [--dry-run]
"""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import default_outdir, safe_outdir, safe_download, log

# Primary source: UCSC CGL Genomics
CACTUS_BASE = "https://cglgenomics.ucsc.edu/data/cactus/"
HAL_URL = "https://cglgenomics.ucsc.edu/data/cactus/241-mammalian-2020v2b.hal"
# Alternative mirror (Zoonomia project page)
HAL_URL_ALT = "https://zoonomiaproject.org/the-data/"

# Species tree used for the alignment
TREE_URL = "https://cglgenomics.ucsc.edu/data/cactus/241-mammalian-2020v2.nh"

EXPECTED_SIZES = {
    "241-mammalian-2020v2b.hal": "~600 GB",
    "241-mammalian-2020v2.nh": "~10 KB",
}


def main():
    parser = argparse.ArgumentParser(
        description="Fetch Zoonomia Cactus 241-species alignment (HAL)")
    parser.add_argument("--outdir", type=Path, default=default_outdir("zoonomia"))
    parser.add_argument("--dry-run", action="store_true",
                        help="Print what would be downloaded without downloading")
    args = parser.parse_args()

    if args.dry_run:
        log.info("DRY RUN — would download:")
        log.info("  %s (%s)", HAL_URL, EXPECTED_SIZES["241-mammalian-2020v2b.hal"])
        log.info("  %s (%s)", TREE_URL, EXPECTED_SIZES["241-mammalian-2020v2.nh"])
        log.info("To: %s", args.outdir)
        log.info("")
        log.info("Alternative source page: %s", HAL_URL_ALT)
        log.info("NOTE: The HAL file is ~600 GB. Ensure sufficient disk space.")
        return

    outdir = safe_outdir(args.outdir)

    # Species tree (small)
    safe_download(TREE_URL, outdir / "241-mammalian-2020v2.nh")

    # HAL alignment (very large)
    log.info("Downloading HAL alignment (~600 GB). This will take a long time.")
    safe_download(HAL_URL, outdir / "241-mammalian-2020v2b.hal",
                  min_size=1_000_000)


if __name__ == "__main__":
    main()
