#!/usr/bin/env python3
"""
Fetch tRNA data via Rfam RF00005 (tRNA family).

GtRNAdb bulk downloads are no longer available (404 as of 2026-03).
Fallback: Rfam RF00005 full alignment contains ~700k tRNA sequences
pre-aligned against the tRNA covariance model with SS_cons annotation.

Usage:
    python fetch/gtrnadb/fetch.py [--outdir data/gtrnadb]
"""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import default_outdir, safe_outdir, safe_download, safe_decompress_gz, log

RFAM_TRNA_FASTA = "https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/fasta_files/RF00005.fa.gz"
RFAM_TRNA_CM = "https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz"
RFAM_TRNA_SEED = "https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/full_alignments/RF00005.sto.gz"


def main():
    parser = argparse.ArgumentParser(description="Fetch tRNA data")
    parser.add_argument("--outdir", type=Path, default=default_outdir("gtrnadb"))
    args = parser.parse_args()

    outdir = safe_outdir(args.outdir)
    align_dir = safe_outdir(outdir / "alignments")

    # RF00005 FASTA
    gz = align_dir / "RF00005.fa.gz"
    if safe_download(RFAM_TRNA_FASTA, gz):
        safe_decompress_gz(gz, align_dir / "RF00005.fa")

    # RF00005 Stockholm alignment
    sto_gz = align_dir / "RF00005.sto.gz"
    safe_download(RFAM_TRNA_SEED, sto_gz)

    log.info("GtRNAdb/tRNA fetch complete: %s", outdir)


if __name__ == "__main__":
    main()
