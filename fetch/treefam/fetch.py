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

# TreeFam download URL
BASE_URL = "https://www.treefam.org/static/download/"
FAMILY_DATA = "treefam_family_data.tar.gz"


def main():
    parser = argparse.ArgumentParser(description="Fetch TreeFam data")
    parser.add_argument("--outdir", type=Path, default=default_outdir("treefam"))
    args = parser.parse_args()

    outdir = safe_outdir(args.outdir)

    safe_download(BASE_URL + FAMILY_DATA, outdir / FAMILY_DATA, min_size=1000)

    # Extract
    import subprocess
    tarball = outdir / FAMILY_DATA
    if not (outdir / "treefam_family_data").exists():
        log.info("Extracting %s...", tarball)
        subprocess.run(["tar", "xzf", str(tarball), "-C", str(outdir)], check=True)

    log.info("TreeFam: %s families available in %s/treefam_family_data/",
             len(list((outdir / "treefam_family_data").glob("*.aa.fasta"))),
             outdir)


if __name__ == "__main__":
    main()
