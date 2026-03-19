#!/usr/bin/env python3
"""
Fetch TreeFam gene family alignments and trees.

Downloads TreeFam data from the EBI archive.

Usage:
    python fetch/treefam/fetch.py
"""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import default_outdir, safe_outdir, safe_download, log

# TreeFam FTP (archived)
BASE_URL = "http://www.treefam.org/static/download/"
FILES = [
    "TreeFam9.tar.gz",
]


def main():
    parser = argparse.ArgumentParser(description="Fetch TreeFam data")
    parser.add_argument("--outdir", type=Path, default=default_outdir("treefam"))
    args = parser.parse_args()

    outdir = safe_outdir(args.outdir)

    for fname in FILES:
        safe_download(BASE_URL + fname, outdir / fname, min_size=1000)

    log.info("TreeFam: download complete. Extract with: tar xzf %s/TreeFam9.tar.gz -C %s/",
             outdir, outdir)


if __name__ == "__main__":
    main()
