#!/usr/bin/env python3
"""
Fetch UCSC hg38 multiz100way multiple alignment MAF files.

Downloads per-chromosome MAF files from the UCSC genome browser.

Usage:
    python fetch/ucsc/assembly/hg38/msa/multiz100way/fetch.py \
        [--outdir data/ucsc/assembly/hg38/msa/multiz100way]
"""

import argparse
import re
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import requests

sys.path.insert(0, str(Path(__file__).resolve().parents[5]))
from common import default_outdir, safe_outdir, safe_download, log, repo_root

BASE_URL = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/multiz100way/maf/"


def list_maf_files() -> list[str]:
    """Scrape directory listing for *.maf.gz files."""
    log.info("Fetching MAF listing from %s", BASE_URL)
    resp = requests.get(BASE_URL, timeout=60)
    resp.raise_for_status()
    pattern = re.compile(r'href="(chr[^"]+\.maf\.gz)"')
    names = sorted(set(pattern.findall(resp.text)))
    log.info("Found %d MAF files", len(names))
    return names


def main():
    parser = argparse.ArgumentParser(description="Fetch hg38 multiz100way MAFs")
    parser.add_argument("--outdir", type=Path,
                        default=repo_root() / "data/ucsc/assembly/hg38/msa/multiz100way")
    parser.add_argument("--jobs", type=int, default=4)
    parser.add_argument("--chroms", nargs="*", default=None,
                        help="Limit to specific chromosomes (e.g., chr1 chr22)")
    args = parser.parse_args()

    outdir = safe_outdir(args.outdir)
    names = list_maf_files()

    if args.chroms:
        prefixes = tuple(c + "." for c in args.chroms)
        names = [n for n in names if n.startswith(prefixes)]
        log.info("Filtered to %d files for chroms: %s", len(names), args.chroms)

    def fetch_one(name):
        return safe_download(BASE_URL + name, outdir / name)

    ok = 0
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futures = {pool.submit(fetch_one, n): n for n in names}
        for f in as_completed(futures):
            if f.result():
                ok += 1

    log.info("multiz100way: %d/%d MAFs fetched", ok, len(names))


if __name__ == "__main__":
    main()
