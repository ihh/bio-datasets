#!/usr/bin/env python3
"""
Fetch SILVA NR99 rRNA alignments and guide trees.

Downloads SSU+LSU NR99 FASTA alignments and taxonomy trees.
Decompresses .gz files after download.

Usage:
    python fetch/silva/fetch.py [--outdir data/silva] [--version 138.2]
"""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import default_outdir, safe_outdir, safe_download, safe_decompress_gz, log

SILVA_FTP = "https://ftp.arb-silva.de/release_{version}"

FILES = {
    "ssu_fasta": "Exports/SILVA_{version}_SSURef_NR99_tax_silva_full_align_trunc.fasta.gz",
    "lsu_fasta": "Exports/SILVA_{version}_LSURef_NR99_tax_silva_full_align_trunc.fasta.gz",
    "ssu_tree": "Exports/taxonomy/tax_slv_ssu_{version}.tre.gz",
    "lsu_tree": "Exports/taxonomy/tax_slv_lsu_{version}.tre.gz",
}


def main():
    parser = argparse.ArgumentParser(description="Fetch SILVA data")
    parser.add_argument("--outdir", type=Path, default=default_outdir("silva/raw"))
    parser.add_argument("--version", default="138.2")
    args = parser.parse_args()

    outdir = safe_outdir(args.outdir)
    base_url = SILVA_FTP.format(version=args.version.replace(".", "_"))

    for key, template in FILES.items():
        rel_path = template.format(version=args.version)
        url = f"{base_url}/{rel_path}"
        gz_name = Path(rel_path).name
        gz_dest = outdir / gz_name

        if safe_download(url, gz_dest):
            # Decompress
            plain_name = gz_name.removesuffix(".gz")
            if plain_name != gz_name:
                safe_decompress_gz(gz_dest, outdir / plain_name)

    log.info("SILVA fetch complete: %s", outdir)


if __name__ == "__main__":
    main()
