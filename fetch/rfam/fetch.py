#!/usr/bin/env python3
"""
Fetch Rfam full alignments (Stockholm) and clan memberships.

Downloads all RF*.sto.gz from Rfam CURRENT, plus clan_membership.txt.gz.
Idempotent: skips files that already exist.

Usage:
    python fetch/rfam/fetch.py [--outdir data/rfam] [--max-families 0] [--jobs 4]
"""

import argparse
import re
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import requests

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import default_outdir, safe_outdir, safe_download, log

BASE_URL = "https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/"
ALIGN_URL = BASE_URL + "full_alignments/"
CLAN_URL = BASE_URL + "database_files/clan_membership.txt.gz"
CM_URL = BASE_URL + "Rfam.cm.gz"


def list_alignment_files(url: str = ALIGN_URL) -> list[str]:
    """Scrape directory listing for RF*.sto.gz filenames."""
    log.info("Fetching alignment listing from %s", url)
    resp = requests.get(url, timeout=60)
    resp.raise_for_status()
    pattern = re.compile(r'href="(RF\d+\.sto(?:\.gz)?)"')
    names = sorted(set(pattern.findall(resp.text)))
    log.info("Found %d alignment files", len(names))
    return names


def main():
    parser = argparse.ArgumentParser(description="Fetch Rfam data")
    parser.add_argument("--outdir", type=Path, default=default_outdir("rfam"))
    parser.add_argument("--max-families", type=int, default=0,
                        help="Limit families (0 = all)")
    parser.add_argument("--jobs", type=int, default=4)
    args = parser.parse_args()

    outdir = safe_outdir(args.outdir)
    align_dir = safe_outdir(outdir / "alignments")

    # Clan membership
    safe_download(CLAN_URL, outdir / "clan_membership.txt.gz")

    # Covariance models
    safe_download(CM_URL, outdir / "Rfam.cm.gz")

    # Alignments
    names = list_alignment_files()
    if args.max_families > 0:
        names = names[:args.max_families]

    def fetch_one(name):
        url = ALIGN_URL + name
        dest = align_dir / name
        return safe_download(url, dest)

    ok = 0
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futures = {pool.submit(fetch_one, n): n for n in names}
        for f in as_completed(futures):
            if f.result():
                ok += 1

    log.info("Rfam: %d/%d alignments fetched", ok, len(names))


if __name__ == "__main__":
    main()
