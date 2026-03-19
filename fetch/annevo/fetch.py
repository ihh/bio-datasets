#!/usr/bin/env python3
"""
Fetch AnnEvo training data: genomes, annotations, and Cactus alignments.

AnnEvo (Hoff et al.) extends Tiberius by incorporating cross-species
evolutionary conservation from Cactus whole-genome alignments. Its training
data includes:

1. Softmasked genome assemblies (FASTA) — same species as Tiberius
2. Gene annotations (GTF) — RefSeq annotations for each species
3. Cactus progressive alignment (HAL) — the Zoonomia 241-mammal alignment,
   from which per-species MAF slices are extracted

The Cactus HAL file is shared with the zoonomia dataset (fetch/zoonomia/).
This script fetches only the genome/annotation components. For the alignment,
run fetch/zoonomia/fetch.py separately, then use hal2maf to extract per-species
alignment windows.

Source: https://github.com/Gaius-Augustus/AnnEvo
Paper: Check the AnnEvo repository for the latest publication reference.

AnnEvo uses the same 37 mammalian species as Tiberius for training, plus
the Cactus alignment provides conservation features across 241 mammals.

Expected total size (genomes + annotations only): ~100-150 GB
  (The Cactus HAL adds ~600 GB — managed by fetch/zoonomia/)

Directory structure after fetch:
  data/annevo/
    genomes/{species}/genome.fa.gz
    annotations/{species}/annotation.gtf.gz
    -> data/zoonomia/241-mammalian-2020v2b.hal  (symlink to shared alignment)

Usage:
    python fetch/annevo/fetch.py [--outdir data/annevo] [--dry-run]
"""

import argparse
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import default_outdir, safe_outdir, log, repo_root

# AnnEvo reuses the same species as Tiberius; import its species list
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "tiberius"))
from fetch import TIBERIUS_SPECIES, ncbi_genome_url, safe_download


def main():
    parser = argparse.ArgumentParser(
        description="Fetch AnnEvo training data (genomes + annotations)")
    parser.add_argument("--outdir", type=Path, default=default_outdir("annevo"))
    parser.add_argument("--dry-run", action="store_true",
                        help="Print what would be downloaded without downloading")
    parser.add_argument("--species", type=str, default=None,
                        help="Fetch only this species (common name)")
    parser.add_argument("--symlink-zoonomia", action="store_true",
                        help="Create symlink to zoonomia HAL alignment")
    args = parser.parse_args()

    species_list = TIBERIUS_SPECIES
    if args.species:
        species_list = [s for s in species_list if s[0] == args.species]
        if not species_list:
            log.error("Unknown species: %s", args.species)
            sys.exit(1)

    if args.dry_run:
        log.info("DRY RUN — AnnEvo training data:")
        log.info("")
        log.info("1. Genomes + annotations: %d species (same as Tiberius)", len(species_list))
        for common, species, acc, asm in species_list:
            log.info("  %-25s %s (%s)", common, acc, asm)
        log.info("")
        log.info("2. Cactus alignment: Zoonomia 241-mammal HAL (~600 GB)")
        log.info("   -> Run: python fetch/zoonomia/fetch.py")
        log.info("")
        log.info("3. Per-species MAF extraction (post-download):")
        log.info("   -> Use hal2maf from HAL tools to extract species-specific alignments")
        log.info("")
        log.info("Estimated total: ~100-150 GB genomes + ~600 GB HAL")
        log.info("To: %s", args.outdir)
        return

    outdir = safe_outdir(args.outdir)
    genomes_dir = safe_outdir(outdir / "genomes")
    annot_dir = safe_outdir(outdir / "annotations")

    # Optionally symlink the zoonomia HAL
    if args.symlink_zoonomia:
        hal_src = repo_root() / "data" / "zoonomia" / "241-mammalian-2020v2b.hal"
        hal_link = outdir / "241-mammalian-2020v2b.hal"
        if hal_src.exists() and not hal_link.exists():
            os.symlink(hal_src, hal_link)
            log.info("Symlinked Zoonomia HAL: %s -> %s", hal_link, hal_src)
        elif not hal_src.exists():
            log.warning("Zoonomia HAL not found at %s — run fetch/zoonomia/fetch.py first",
                        hal_src)

    # Download genomes and annotations (same as Tiberius)
    ok_genomes, ok_annots = 0, 0
    for common, species, acc, asm in species_list:
        base = ncbi_genome_url(acc, asm)
        prefix = f"{acc}_{asm}"

        sp_genome = safe_outdir(genomes_dir / common)
        genome_url = f"{base}/{prefix}_genomic.fna.gz"
        if safe_download(genome_url, sp_genome / "genome.fa.gz", min_size=1000):
            ok_genomes += 1

        sp_annot = safe_outdir(annot_dir / common)
        gtf_url = f"{base}/{prefix}_genomic.gtf.gz"
        if safe_download(gtf_url, sp_annot / "annotation.gtf.gz", min_size=1000):
            ok_annots += 1

    log.info("AnnEvo: %d/%d genomes, %d/%d annotations fetched",
             ok_genomes, len(species_list), ok_annots, len(species_list))
    log.info("NOTE: For the Cactus alignment, run: python fetch/zoonomia/fetch.py")


if __name__ == "__main__":
    main()
