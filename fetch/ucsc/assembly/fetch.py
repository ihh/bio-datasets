#!/usr/bin/env python3
"""
Fetch UCSC reference genome assembly FASTA files.

Downloads chromFa.tar.gz or {genome}.fa.gz from the UCSC genome browser
and produces a single {genome}.fa file.

Usage:
    python fetch/ucsc/assembly/fetch.py --genome ce11
    python fetch/ucsc/assembly/fetch.py --genome hg19
    python fetch/ucsc/assembly/fetch.py --genome hg38
    python fetch/ucsc/assembly/fetch.py --genome ce11 hg19 hg38
"""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from common import safe_outdir, safe_download, safe_decompress_gz, log, repo_root

BIGZIPS_URL = "https://hgdownload.soe.ucsc.edu/goldenPath/{genome}/bigZips/{genome}.fa.gz"


def fetch_genome(genome: str, base_outdir: Path) -> bool:
    """Fetch a single reference genome FASTA."""
    outdir = safe_outdir(base_outdir / genome)
    fa_path = outdir / f"{genome}.fa"
    gz_path = outdir / f"{genome}.fa.gz"

    # If the final .fa already exists (or is a symlink), skip entirely
    if fa_path.is_symlink():
        log.info("SKIP (symlink): %s", fa_path)
        return True
    if fa_path.exists() and fa_path.stat().st_size > 0:
        log.info("SKIP (exists): %s (%d bytes)", fa_path.name, fa_path.stat().st_size)
        return True

    url = BIGZIPS_URL.format(genome=genome)
    if not safe_download(url, gz_path):
        return False

    return safe_decompress_gz(gz_path, fa_path)


def main():
    parser = argparse.ArgumentParser(
        description="Fetch UCSC reference genome FASTA files"
    )
    parser.add_argument(
        "--genome", nargs="+", required=True,
        help="UCSC genome identifier(s), e.g. ce11 hg19 hg38"
    )
    parser.add_argument(
        "--outdir", type=Path,
        default=repo_root() / "data" / "ucsc" / "assembly",
    )
    args = parser.parse_args()

    outdir = safe_outdir(args.outdir)
    ok = 0
    for genome in args.genome:
        if fetch_genome(genome, outdir):
            ok += 1

    log.info("UCSC assemblies: %d/%d fetched", ok, len(args.genome))


if __name__ == "__main__":
    main()
