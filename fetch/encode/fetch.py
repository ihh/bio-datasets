#!/usr/bin/env python3
"""
Fetch ENCODE BigWig tracks for Borzoi-style multi-track prediction.

Downloads a curated subset of ENCODE BigWig files covering key functional
genomics assays, starting with K562 and GM12878 cell lines. These tracks
are used for aa-msa-mamba Task 4 (Borzoi-style track prediction).

Track types:
  - CAGE:         gene expression / TSS identification
  - RNA-seq:      transcript quantification (total RNA, polyA)
  - ATAC-seq:     chromatin accessibility
  - H3K4me3:      active promoter mark (ChIP-seq)
  - CTCF:         insulator binding (ChIP-seq)

Expected total size: ~15-25 GB

Directory structure after fetch:
  data/encode/bigwig/
    cage/
    rnaseq/
    atacseq/
    h3k4me3/
    ctcf/

Usage:
    python fetch/encode/fetch.py [--outdir data/encode/bigwig] [--dry-run]
    python fetch/encode/fetch.py --track cage     # fetch only CAGE tracks
"""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import default_outdir, safe_outdir, safe_download, log

ENCODE_BASE = "https://www.encodeproject.org/files"

# Curated BigWig tracks from ENCODE.
# Format: (subdir, filename, ENCODE accession, cell_line, assay, description, approx_size)
#
# These are signal p-value or fold-change-over-control BigWig files from
# the ENCODE portal for well-characterized cell lines (K562, GM12878).
# Each URL is: https://www.encodeproject.org/files/{accession}/@@download/{accession}.bigWig

TRACKS = [
    # CAGE (ENCODE CAGE / RAMPAGE as proxy)
    ("cage", "ENCFF542GMN.bigWig", "ENCFF542GMN", "K562",
     "CAGE", "CAGE signal, K562, plus strand", "~2 GB"),
    ("cage", "ENCFF347MKI.bigWig", "ENCFF347MKI", "GM12878",
     "CAGE", "CAGE signal, GM12878, plus strand", "~2 GB"),

    # RNA-seq (total RNA)
    ("rnaseq", "ENCFF232GIC.bigWig", "ENCFF232GIC", "K562",
     "RNA-seq", "RNA-seq signal, K562, total RNA, plus strand", "~2 GB"),
    ("rnaseq", "ENCFF014OGT.bigWig", "ENCFF014OGT", "GM12878",
     "RNA-seq", "RNA-seq signal, GM12878, total RNA, plus strand", "~2 GB"),

    # ATAC-seq
    ("atacseq", "ENCFF667MDI.bigWig", "ENCFF667MDI", "K562",
     "ATAC-seq", "ATAC-seq signal, K562, fold change over control", "~2 GB"),
    ("atacseq", "ENCFF832UPC.bigWig", "ENCFF832UPC", "GM12878",
     "ATAC-seq", "ATAC-seq signal, GM12878, fold change over control", "~2 GB"),

    # H3K4me3 ChIP-seq
    ("h3k4me3", "ENCFF473YHH.bigWig", "ENCFF473YHH", "K562",
     "H3K4me3", "H3K4me3 ChIP-seq signal, K562, fold change over control", "~2 GB"),
    ("h3k4me3", "ENCFF133BES.bigWig", "ENCFF133BES", "GM12878",
     "H3K4me3", "H3K4me3 ChIP-seq signal, GM12878, fold change over control", "~2 GB"),

    # CTCF ChIP-seq
    ("ctcf", "ENCFF534DCJ.bigWig", "ENCFF534DCJ", "K562",
     "CTCF", "CTCF ChIP-seq signal, K562, fold change over control", "~2 GB"),
    ("ctcf", "ENCFF213VOZ.bigWig", "ENCFF213VOZ", "GM12878",
     "CTCF", "CTCF ChIP-seq signal, GM12878, fold change over control", "~2 GB"),
]


def track_url(accession: str) -> str:
    """Construct ENCODE download URL for a BigWig file."""
    return f"{ENCODE_BASE}/{accession}/@@download/{accession}.bigWig"


def main():
    parser = argparse.ArgumentParser(
        description="Fetch ENCODE BigWig tracks for Borzoi-style prediction")
    parser.add_argument(
        "--outdir", type=Path,
        default=default_outdir("encode") / "bigwig")
    parser.add_argument(
        "--dry-run", action="store_true",
        help="Print what would be downloaded without downloading")
    parser.add_argument(
        "--track", type=str, default=None,
        choices=sorted(set(t[0] for t in TRACKS)),
        help="Fetch only this track type")
    args = parser.parse_args()

    tracks = TRACKS
    if args.track:
        tracks = [t for t in tracks if t[0] == args.track]

    if args.dry_run:
        log.info("DRY RUN — ENCODE BigWig tracks (%d files):", len(tracks))
        log.info("Target directory: %s", args.outdir)
        log.info("")
        for subdir, fname, acc, cell, assay, desc, size in tracks:
            url = track_url(acc)
            log.info("  %s/%-25s  %s  %s", subdir, fname, size, desc)
            log.info("    %s", url)
        log.info("")
        log.info("Estimated total: ~15-25 GB (10 BigWig files)")
        log.info("")
        log.info("Note: ENCODE accessions above are representative examples.")
        log.info("Verify accessions at https://www.encodeproject.org/ for")
        log.info("your specific analysis requirements.")
        return

    outdir = safe_outdir(args.outdir)

    ok = 0
    for subdir, fname, acc, cell, assay, desc, size in tracks:
        track_dir = safe_outdir(outdir / subdir)
        url = track_url(acc)
        dest = track_dir / fname
        log.info("Fetching %s ...", desc)
        if safe_download(url, dest, min_size=1000):
            ok += 1

    log.info("ENCODE: %d/%d BigWig tracks fetched", ok, len(tracks))


if __name__ == "__main__":
    main()
