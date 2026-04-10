#!/usr/bin/env python3
"""
Fetch Pfam FULL alignments (Stockholm format, all members).

Downloads the bulk Pfam-A.full.gz from EBI FTP. This file is ~15-20 GB
compressed and contains all Pfam families with their full (non-seed)
member alignments.

Usage:
    python fetch/pfam/fetch_full.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import safe_download, log

BASE_URL = "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/"
FULL_URL = BASE_URL + "Pfam-A.full.gz"


def main():
    data_dir = Path(__file__).resolve().parent.parent.parent / "data"
    outdir = data_dir / "pfam" / "full"
    outdir.mkdir(parents=True, exist_ok=True)

    dest = outdir / "Pfam-A.full.gz"
    log.info("Downloading Pfam-A.full.gz (~15-20 GB) to %s", dest)
    ok = safe_download(FULL_URL, dest, min_size=1_000_000_000)  # expect >1 GB
    if ok:
        size_gb = dest.stat().st_size / (1024**3)
        log.info("Download complete: %.2f GB", size_gb)
    else:
        log.error("Download failed")
        sys.exit(1)


if __name__ == "__main__":
    main()
