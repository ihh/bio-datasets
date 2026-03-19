#!/usr/bin/env python3
"""
Fetch genome assemblies from NCBI GenBank by accession.

Downloads genomic FASTA files from the NCBI FTP server given GCA/GCF accessions.

Usage:
    python fetch/ncbi/assembly/fetch.py --accessions GCA_016989095.1 GCA_016989105.1
    python fetch/ncbi/assembly/fetch.py --accessions-file accessions.txt --jobs 4
"""

import argparse
import re
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import requests

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from common import safe_outdir, safe_download, safe_decompress_gz, log, repo_root

NCBI_FTP = "https://ftp.ncbi.nlm.nih.gov"
DATASETS_SUMMARY = "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{accession}"


def resolve_ftp_path(accession: str) -> str | None:
    """Resolve a GCA/GCF accession to its NCBI FTP directory URL.

    Uses the standard path convention:
      ftp.ncbi.nlm.nih.gov/genomes/all/GCA/016/989/095/GCA_016989095.1_*/
    """
    m = re.match(r"(GC[AF])_(\d{3})(\d{3})(\d{3})\.(\d+)", accession)
    if not m:
        log.error("Invalid accession format: %s", accession)
        return None
    prefix, d1, d2, d3, version = m.groups()
    base = f"{NCBI_FTP}/genomes/all/{prefix}/{d1}/{d2}/{d3}/"

    # List directory to find the full assembly name
    log.info("Resolving FTP path for %s", accession)
    try:
        resp = requests.get(base, timeout=60)
        resp.raise_for_status()
    except requests.RequestException as e:
        log.error("Failed to list %s: %s", base, e)
        return None

    # Find directory matching this accession
    pattern = re.compile(rf'href="({re.escape(accession)}_[^"]+)/"')
    matches = pattern.findall(resp.text)
    if not matches:
        log.error("No matching directory for %s at %s", accession, base)
        return None

    asm_name = matches[0]
    return f"{base}{asm_name}/{asm_name}_genomic.fna.gz"


def fetch_assembly(accession: str, outdir: Path) -> bool:
    """Fetch a single assembly FASTA by accession."""
    asm_dir = safe_outdir(outdir / accession)
    gz_name = f"{accession}_genomic.fna.gz"
    fna_name = f"{accession}_genomic.fna"
    gz_path = asm_dir / gz_name
    fna_path = asm_dir / fna_name

    # If decompressed file exists, skip
    if fna_path.is_symlink():
        log.info("SKIP (symlink): %s", fna_path)
        return True
    if fna_path.exists() and fna_path.stat().st_size > 0:
        log.info("SKIP (exists): %s (%d bytes)", fna_path.name, fna_path.stat().st_size)
        return True

    url = resolve_ftp_path(accession)
    if not url:
        return False

    if not safe_download(url, gz_path):
        return False

    return safe_decompress_gz(gz_path, fna_path)


def main():
    parser = argparse.ArgumentParser(
        description="Fetch NCBI GenBank/RefSeq genome assemblies"
    )
    parser.add_argument(
        "--accessions", nargs="*", default=[],
        help="GCA/GCF accession(s), e.g. GCA_016989095.1"
    )
    parser.add_argument(
        "--accessions-file", type=Path, default=None,
        help="File with one accession per line"
    )
    parser.add_argument(
        "--outdir", type=Path,
        default=repo_root() / "data" / "ncbi" / "assembly",
    )
    parser.add_argument("--jobs", type=int, default=2)
    args = parser.parse_args()

    accessions = list(args.accessions)
    if args.accessions_file:
        with open(args.accessions_file) as f:
            accessions.extend(
                line.strip() for line in f if line.strip() and not line.startswith("#")
            )

    if not accessions:
        parser.error("No accessions specified (use --accessions or --accessions-file)")

    outdir = safe_outdir(args.outdir)
    ok = 0

    def fetch_one(acc):
        return fetch_assembly(acc, outdir)

    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futures = {pool.submit(fetch_one, a): a for a in accessions}
        for f in as_completed(futures):
            if f.result():
                ok += 1

    log.info("NCBI assemblies: %d/%d fetched", ok, len(accessions))


if __name__ == "__main__":
    main()
