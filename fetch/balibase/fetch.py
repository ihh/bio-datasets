#!/usr/bin/env python3
"""
Fetch BAliBASE alignment benchmark.

BAliBASE (Benchmark Alignment dataBASE) provides reference alignments
for evaluating multiple sequence alignment methods. We use the drive5
mirror, which bundles BAliBASE v3 alongside PREFAB v4, OXBENCH, and
SABRE in a single 21 MB tarball, all converted to standard FASTA.

The original IGBMC FTP distribution (ftp-igbmc.u-strasbg.fr/pub/BAliBASE3)
and the LBGI HTTP distribution (lbgi.fr/balibase/BAliBASE_R1-5.tar.gz)
have both been unreliable in our experience; the drive5 bundle is what
we use in practice.

After download, the tarball is extracted in place and yields:
    bench/bali3/in/   FASTA inputs (one .tfa per benchmark)
    bench/bali3/ref/  reference alignments
    bench/prefab4/...
    bench/oxbench/...
    bench/sabre/...

Usage:
    python fetch/balibase/fetch.py
"""

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import default_outdir, safe_outdir, safe_download, log


URL = "https://drive5.com/bench/bench.tar.gz"
ARCHIVE_NAME = "bench.tar.gz"
EXTRACT_MARKER = "bench1.0"   # top-level dir inside the drive5 tarball


def main():
    parser = argparse.ArgumentParser(description="Fetch BAliBASE benchmark")
    parser.add_argument("--outdir", type=Path,
                            default=default_outdir("balibase"))
    parser.add_argument("--no-extract", action="store_true",
                            help="Skip extraction; leave the tarball as-is.")
    args = parser.parse_args()

    outdir = safe_outdir(args.outdir)
    archive = outdir / ARCHIVE_NAME
    safe_download(URL, archive, min_size=1_000_000)

    extracted = outdir / EXTRACT_MARKER
    if args.no_extract:
        log.info("Skipping extraction (--no-extract set).")
    elif extracted.exists():
        log.info("Already extracted at %s; skipping.", extracted)
    else:
        log.info("Extracting %s ...", archive)
        subprocess.run(["tar", "xzf", str(archive), "-C", str(outdir)],
                          check=True)
        log.info("Extracted to %s", extracted)

    log.info("BAliBASE: ready at %s", outdir)


if __name__ == "__main__":
    main()
