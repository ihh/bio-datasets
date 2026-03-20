#!/usr/bin/env python3
"""
Fetch ProteinGym DMS benchmark data.

Downloads DMS indel and substitution variant datasets, reference files,
and optionally MSA files for per-protein model training.

Usage:
    python fetch/proteingym/fetch.py                    # indel data only (small)
    python fetch/proteingym/fetch.py --substitutions    # also substitution data
    python fetch/proteingym/fetch.py --msas             # also MSA files (5.2 GB)
"""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import default_outdir, safe_outdir, safe_download, log

BASE_URL = "https://marks.hms.harvard.edu/proteingym/"
GITHUB_RAW = "https://raw.githubusercontent.com/OATML-Markslab/ProteinGym/main/"

INDEL_FILES = [
    ("DMS_indels.zip", BASE_URL + "DMS_indels.zip", 200_000_000),
]

SUBSTITUTION_FILES = [
    ("DMS_substitutions.zip", BASE_URL + "DMS_substitutions.zip", 500_000_000),
]

REFERENCE_FILES = [
    ("DMS_indels.csv", GITHUB_RAW + "reference_files/DMS_indels.csv", 1000),
    ("DMS_substitutions.csv", GITHUB_RAW + "reference_files/DMS_substitutions.csv", 1000),
]

MSA_FILES = [
    ("DMS_MSAs.zip", BASE_URL + "DMS_MSAs.zip", 5_000_000_000),
]


def main():
    parser = argparse.ArgumentParser(description="Fetch ProteinGym data")
    parser.add_argument("--outdir", type=Path, default=default_outdir("proteingym"))
    parser.add_argument("--substitutions", action="store_true",
                        help="Also download DMS substitution data (~500 MB)")
    parser.add_argument("--msas", action="store_true",
                        help="Also download MSA files (~5.2 GB)")
    args = parser.parse_args()

    outdir = safe_outdir(args.outdir)

    # Reference files (always)
    log.info("Downloading reference files...")
    for name, url, min_size in REFERENCE_FILES:
        safe_download(url, outdir / name, min_size=min_size)

    # Indel data (always)
    log.info("Downloading DMS indel data...")
    for name, url, min_size in INDEL_FILES:
        safe_download(url, outdir / name, min_size=min_size)

    # Substitution data (optional)
    if args.substitutions:
        log.info("Downloading DMS substitution data...")
        for name, url, min_size in SUBSTITUTION_FILES:
            safe_download(url, outdir / name, min_size=min_size)

    # MSA files (optional, large)
    if args.msas:
        log.info("Downloading MSA files (5.2 GB)...")
        for name, url, min_size in MSA_FILES:
            safe_download(url, outdir / name, min_size=min_size)

    log.info("Done. Unzip with: cd %s && unzip DMS_indels.zip", outdir)


if __name__ == "__main__":
    main()
