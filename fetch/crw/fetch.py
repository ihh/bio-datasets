#!/usr/bin/env python3
"""
Fetch CRW (Comparative RNA Web) alignments and BPSEQ structure files.

CRW provides hand-curated rRNA alignments with experimentally verified
secondary structures.

Usage:
    python fetch/crw/fetch.py [--outdir data/crw]
"""

import argparse
import shutil
import sys
import zipfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import default_outdir, safe_outdir, safe_download, log

CRW_BASE = "https://crw-site.chemistry.gatech.edu/DAT/3C"

ALIGNMENT_FILES = {
    "16S_A": "Alignment/Files/16S/16S.A.ALL.alnfasta.zip",
    "16S_B": "Alignment/Files/16S/16S.B.ALL.alnfasta.zip",
    "16S_E": "Alignment/Files/16S/16S.E.ALL.alnfasta.zip",
    "16S_C": "Alignment/Files/16S/16S.C.ALL.alnfasta.zip",
    "16S_M": "Alignment/Files/16S/16S.M.ALL.alnfasta.zip",
    "23S_A": "Alignment/Files/23S/23S.A.ALL.alnfasta.zip",
    "23S_B": "Alignment/Files/23S/23S.B.ALL.alnfasta.zip",
    "23S_E": "Alignment/Files/23S/23S.E.ALL.alnfasta.zip",
}

BPSEQ_FILES = {
    "16S_bpseq": "Structures/Whole_Structures/16S/bpseq.zip",
    "23S_bpseq": "Structures/Whole_Structures/23S/bpseq.zip",
}


def safe_unzip(zip_path: Path, dest_dir: Path):
    """Extract zip if dest_dir is empty or doesn't have expected files."""
    if dest_dir.is_symlink():
        log.info("SKIP unzip (symlink): %s", dest_dir)
        return
    if any(dest_dir.iterdir()):
        log.info("SKIP unzip (non-empty): %s", dest_dir)
        return
    log.info("Extracting: %s → %s", zip_path.name, dest_dir)
    with zipfile.ZipFile(zip_path) as zf:
        zf.extractall(dest_dir)


def main():
    parser = argparse.ArgumentParser(description="Fetch CRW data")
    parser.add_argument("--outdir", type=Path, default=default_outdir("crw"))
    args = parser.parse_args()

    outdir = safe_outdir(args.outdir)
    align_dir = safe_outdir(outdir / "alignments")
    bpseq_dir = safe_outdir(outdir / "bpseq")

    # Download and extract alignments
    for key, rel_path in ALIGNMENT_FILES.items():
        url = f"{CRW_BASE}/{rel_path}"
        zip_dest = outdir / f"{key}.zip"
        if safe_download(url, zip_dest):
            extract_dir = align_dir / key
            extract_dir.mkdir(exist_ok=True)
            safe_unzip(zip_dest, extract_dir)

    # Download and extract BPSEQ structures
    for key, rel_path in BPSEQ_FILES.items():
        url = f"{CRW_BASE}/{rel_path}"
        zip_dest = outdir / f"{key}.zip"
        if safe_download(url, zip_dest):
            extract_dir = bpseq_dir / key
            extract_dir.mkdir(exist_ok=True)
            safe_unzip(zip_dest, extract_dir)

    log.info("CRW fetch complete: %s", outdir)


if __name__ == "__main__":
    main()
