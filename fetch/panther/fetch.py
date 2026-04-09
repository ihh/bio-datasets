#!/usr/bin/env python3
"""
Fetch PANTHER 19.0 TreeGrafter data (MSAs + trees for ~15k protein families).

Downloads the TreeGrafter data package from pantherdb.org FTP.

Usage:
    python fetch/panther/fetch.py
    python fetch/panther/fetch.py --outdir /path/to/data/panther
"""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import default_outdir, safe_outdir, safe_download, log

# PANTHER 19.0 TreeGrafter data (~1-2 GB compressed)
PANTHER_URL = "https://data.pantherdb.org/ftp/downloads/TreeGrafter/PANTHER19.0_data.tar.gz"


def main():
    parser = argparse.ArgumentParser(description="Fetch PANTHER TreeGrafter data")
    parser.add_argument("--outdir", type=Path, default=default_outdir("panther"))
    args = parser.parse_args()

    outdir = safe_outdir(args.outdir)
    dest = outdir / "PANTHER19.0_data.tar.gz"

    log.info("Downloading PANTHER 19.0 TreeGrafter data to %s", outdir)
    ok = safe_download(PANTHER_URL, dest, min_size=100_000_000)  # expect >100MB

    if ok:
        log.info("Download complete: %s (%.1f MB)",
                 dest, dest.stat().st_size / 1e6)
    else:
        log.error("Download failed")
        sys.exit(1)


if __name__ == "__main__":
    main()
