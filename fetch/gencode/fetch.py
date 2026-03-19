#!/usr/bin/env python3
"""
Fetch GENCODE human gene annotations for hg38.

Downloads GENCODE v45 GTF annotation files for use in gene-structure
prediction tasks (e.g. aa-msa-mamba Task 3).

Files downloaded:
  - gencode.v45.annotation.gtf.gz          (~50 MB)  comprehensive annotation
  - gencode.v45.basic.annotation.gtf.gz     (~25 MB)  basic annotation subset

Expected total size: ~75 MB compressed

Directory structure after fetch:
  data/gencode/human/
    gencode.v45.annotation.gtf.gz
    gencode.v45.basic.annotation.gtf.gz

Usage:
    python fetch/gencode/fetch.py [--outdir data/gencode/human] [--dry-run]
"""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import default_outdir, safe_outdir, safe_download, log

GENCODE_VERSION = "45"
GENCODE_BASE = (
    f"https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human"
    f"/release_{GENCODE_VERSION}"
)

# Files to download: (filename, description, approximate size)
FILES = [
    (
        f"gencode.v{GENCODE_VERSION}.annotation.gtf.gz",
        "Comprehensive gene annotation (GTF)",
        "~50 MB",
    ),
    (
        f"gencode.v{GENCODE_VERSION}.basic.annotation.gtf.gz",
        "Basic gene annotation (GTF)",
        "~25 MB",
    ),
]


def main():
    parser = argparse.ArgumentParser(
        description="Fetch GENCODE human gene annotations")
    parser.add_argument(
        "--outdir", type=Path,
        default=default_outdir("gencode") / "human")
    parser.add_argument(
        "--dry-run", action="store_true",
        help="Print what would be downloaded without downloading")
    parser.add_argument(
        "--version", type=str, default=GENCODE_VERSION,
        help=f"GENCODE release version (default: {GENCODE_VERSION})")
    args = parser.parse_args()

    # Allow overriding version
    version = args.version
    base_url = (
        f"https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human"
        f"/release_{version}"
    )
    files = [
        (f"gencode.v{version}.annotation.gtf.gz",
         "Comprehensive gene annotation (GTF)", "~50 MB"),
        (f"gencode.v{version}.basic.annotation.gtf.gz",
         "Basic gene annotation (GTF)", "~25 MB"),
    ]

    if args.dry_run:
        log.info("DRY RUN — GENCODE v%s human annotations", version)
        log.info("Target directory: %s", args.outdir)
        log.info("")
        for fname, desc, size in files:
            url = f"{base_url}/{fname}"
            log.info("  %-50s  %s", fname, size)
            log.info("    %s", url)
        log.info("")
        log.info("Estimated total: ~75 MB compressed")
        return

    outdir = safe_outdir(args.outdir)

    ok = 0
    for fname, desc, size in files:
        url = f"{base_url}/{fname}"
        dest = outdir / fname
        log.info("Fetching %s ...", desc)
        if safe_download(url, dest, min_size=1000):
            ok += 1

    log.info("GENCODE: %d/%d files fetched", ok, len(files))


if __name__ == "__main__":
    main()
