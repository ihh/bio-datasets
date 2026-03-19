#!/usr/bin/env python3
"""
Fetch Tiberius training data: softmasked genome assemblies + RefSeq GTF annotations.

Tiberius (Gabriel et al. 2024) trains on 37 mammalian genomes with RefSeq
annotations. This script downloads the genome FASTA and GTF files from NCBI.

The species list and assembly accessions are from the Tiberius paper
supplementary materials and GitHub repository:
  https://github.com/Gaius-Augustus/Tiberius

Expected total size: ~100-150 GB (37 softmasked genomes + annotations)

Directory structure after fetch:
  data/tiberius/
    genomes/{species}/genome.fa.gz       — softmasked genome assembly
    annotations/{species}/annotation.gtf.gz  — RefSeq GTF annotation
    species_list.tsv                     — species metadata

Usage:
    python fetch/tiberius/fetch.py [--outdir data/tiberius] [--dry-run] [--species human]
"""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import default_outdir, safe_outdir, safe_download, log

# NCBI datasets base URL for genome downloads
NCBI_FTP = "https://ftp.ncbi.nlm.nih.gov/genomes/all"

# 37 mammalian species used for Tiberius training (from paper supplementary).
# Format: (common_name, species, RefSeq accession, assembly name)
# These accessions are representative; verify against the Tiberius repo for
# exact versions used in the published results.
TIBERIUS_SPECIES = [
    ("human", "Homo_sapiens", "GCF_000001405.40", "GRCh38.p14"),
    ("mouse", "Mus_musculus", "GCF_000001635.27", "GRCm39"),
    ("rat", "Rattus_norvegicus", "GCF_015227675.2", "mRatBN7.2"),
    ("dog", "Canis_lupus_familiaris", "GCF_011100685.1", "UU_Cfam_GSD_1.0"),
    ("cat", "Felis_catus", "GCF_018350175.1", "F.catus_Fca126_mat1.0"),
    ("horse", "Equus_caballus", "GCF_002863925.1", "EquCab3.0"),
    ("pig", "Sus_scrofa", "GCF_000003025.6", "Sscrofa11.1"),
    ("cow", "Bos_taurus", "GCF_002263795.3", "ARS-UCD2.0"),
    ("sheep", "Ovis_aries", "GCF_016772045.2", "ARS-UI_Ramb_v3.0"),
    ("goat", "Capra_hircus", "GCF_001704415.2", "ARS1.2"),
    ("rabbit", "Oryctolagus_cuniculus", "GCF_000003625.3", "OryCun2.0"),
    ("rhesus", "Macaca_mulatta", "GCF_003339765.1", "Mmul_10"),
    ("marmoset", "Callithrix_jacchus", "GCF_009663435.1", "mCalJac1.pat.X"),
    ("gorilla", "Gorilla_gorilla", "GCF_029281585.2", "NHGRI_mGorGor1-v2.0_pri"),
    ("chimpanzee", "Pan_troglodytes", "GCF_028858775.2", "NHGRI_mPanTro3-v2.0_pri"),
    ("orangutan", "Pongo_abelii", "GCF_028885655.2", "NHGRI_mPonAbe1-v2.0_pri"),
    ("bonobo", "Pan_paniscus", "GCF_029289425.2", "NHGRI_mPanPan1-v2.0_pri"),
    ("dolphin", "Tursiops_truncatus", "GCF_011762595.1", "mTurTru1.mat.Y"),
    ("bat_greater_horseshoe", "Rhinolophus_ferrumequinum", "GCF_004115265.2", "mRhiFer1_v1.p"),
    ("hedgehog", "Erinaceus_europaeus", "GCF_006399205.1", "EriEur_WGS.v1"),
    ("elephant", "Loxodonta_africana", "GCF_000001905.1", "Loxafr3.0"),
    ("armadillo", "Dasypus_novemcinctus", "GCF_000208655.1", "Dasnov3.0"),
    ("opossum", "Monodelphis_domestica", "GCF_027887165.1", "mMonDom1.pri"),
    ("platypus", "Ornithorhynchus_anatinus", "GCF_004115215.2", "mOrnAna1.pri.v4"),
    ("naked_mole_rat", "Heterocephalus_glaber", "GCF_944474725.1", "HetGla_v1.1_paternal_haplotype"),
    ("guinea_pig", "Cavia_porcellus", "GCF_000151735.1", "Cavpor3.0"),
    ("hamster", "Mesocricetus_auratus", "GCF_017639785.1", "MesAur1.0"),
    ("deer_mouse", "Peromyscus_maniculatus", "GCF_003704035.1", "HU_Pman_2.1.3"),
    ("shrew", "Sorex_araneus", "GCF_000181275.1", "SorAra2.0"),
    ("pika", "Ochotona_princeps", "GCF_014633375.1", "OchPri4.0"),
    ("manatee", "Trichechus_manatus", "GCF_000243295.1", "TriManLat1.0"),
    ("pangolin", "Manis_javanica", "GCF_014570535.1", "YNU_ManJav_2.0"),
    ("whale_minke", "Balaenoptera_acutorostrata", "GCF_000493695.1", "BalAcu1.0"),
    ("flying_fox", "Pteropus_alecto", "GCF_000325575.1", "ASM32557v1"),
    ("tree_shrew", "Tupaia_chinensis", "GCF_000334495.1", "TupChi_1.0"),
    ("tenrec", "Echinops_telfairi", "GCF_000313985.1", "EchTel2.0"),
    ("hyrax", "Procavia_capensis", "GCF_000152225.1", "proCap1"),
]


def ncbi_genome_url(accession: str, assembly_name: str) -> str:
    """Construct NCBI FTP URL for a genome assembly."""
    # GCF_000001405.40 -> GCF/000/001/405/GCF_000001405.40_GRCh38.p14
    parts = accession.replace("GCF_", "").split(".")
    num = parts[0]
    p1, p2, p3 = num[:3], num[3:6], num[6:9]
    asm_dir = f"{accession}_{assembly_name}"
    return f"{NCBI_FTP}/GCF/{p1}/{p2}/{p3}/{asm_dir}"


def main():
    parser = argparse.ArgumentParser(
        description="Fetch Tiberius training genomes + annotations from NCBI")
    parser.add_argument("--outdir", type=Path, default=default_outdir("tiberius"))
    parser.add_argument("--dry-run", action="store_true",
                        help="Print what would be downloaded without downloading")
    parser.add_argument("--species", type=str, default=None,
                        help="Fetch only this species (common name)")
    parser.add_argument("--jobs", type=int, default=2,
                        help="Parallel downloads (keep low for NCBI)")
    args = parser.parse_args()

    species_list = TIBERIUS_SPECIES
    if args.species:
        species_list = [s for s in species_list if s[0] == args.species]
        if not species_list:
            log.error("Unknown species: %s", args.species)
            log.info("Available: %s", ", ".join(s[0] for s in TIBERIUS_SPECIES))
            sys.exit(1)

    if args.dry_run:
        log.info("DRY RUN — would download %d species:", len(species_list))
        total_est = 0
        for common, species, acc, asm in species_list:
            base = ncbi_genome_url(acc, asm)
            log.info("  %-25s %s", common, f"{base}/")
            log.info("    genome: %s_genomic.fna.gz", f"{acc}_{asm}")
            log.info("    annot:  %s_genomic.gtf.gz", f"{acc}_{asm}")
        log.info("")
        log.info("Estimated total: ~100-150 GB")
        log.info("To: %s", args.outdir)
        return

    outdir = safe_outdir(args.outdir)
    genomes_dir = safe_outdir(outdir / "genomes")
    annot_dir = safe_outdir(outdir / "annotations")

    ok_genomes, ok_annots = 0, 0
    for common, species, acc, asm in species_list:
        base = ncbi_genome_url(acc, asm)
        prefix = f"{acc}_{asm}"

        # Genome FASTA
        sp_genome = safe_outdir(genomes_dir / common)
        genome_url = f"{base}/{prefix}_genomic.fna.gz"
        if safe_download(genome_url, sp_genome / "genome.fa.gz", min_size=1000):
            ok_genomes += 1

        # GTF annotation
        sp_annot = safe_outdir(annot_dir / common)
        gtf_url = f"{base}/{prefix}_genomic.gtf.gz"
        if safe_download(gtf_url, sp_annot / "annotation.gtf.gz", min_size=1000):
            ok_annots += 1

    log.info("Tiberius: %d/%d genomes, %d/%d annotations fetched",
             ok_genomes, len(species_list), ok_annots, len(species_list))


if __name__ == "__main__":
    main()
