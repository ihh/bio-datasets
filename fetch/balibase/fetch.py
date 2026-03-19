#!/usr/bin/env python3
"""
Fetch BAliBASE alignment benchmark.

BAliBASE (Benchmark Alignment dataBASE) provides reference alignments
for evaluating multiple sequence alignment methods.

Usage:
    python fetch/balibase/fetch.py
"""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import default_outdir, safe_outdir, safe_download, log

# BAliBASE v3 from the LBB
BASE_URL = "https://lbgi.fr/balibase/"
FILES = [
    "BAliBASE_R1-5.tar.gz",
]


def main():
    parser = argparse.ArgumentParser(description="Fetch BAliBASE benchmark")
    parser.add_argument("--outdir", type=Path, default=default_outdir("balibase"))
    args = parser.parse_args()

    outdir = safe_outdir(args.outdir)

    for fname in FILES:
        safe_download(BASE_URL + fname, outdir / fname, min_size=1000)

    log.info("BAliBASE: download complete. Extract with: tar xzf %s/%s -C %s/",
             outdir, FILES[0], outdir)


if __name__ == "__main__":
    main()
