"""Shared configuration and helpers for the MVA Track 1 pipeline."""
from __future__ import annotations

import os
from pathlib import Path

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
ROOT_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = ROOT_DIR / "data"
VCF = DATA_DIR / "WGS_EX2312012_HGWCNDSX7.vcf.gz"
TBI = DATA_DIR / "WGS_EX2312012_HGWCNDSX7.vcf.gz.tbi"
FASTQ_PATTERN = DATA_DIR / "WGS_EX2312012_HGWCNDSX7_S16_L00{lanes}_R{read}_001.fastq.gz"
NOTES_DIR = ROOT_DIR / "notes"
OUT_DIR = NOTES_DIR  # analysis outputs land here
PROBAND_ID = "WGS_EX2312012"   # from VCF #CHROM sample column — the schema's proband_id

os.makedirs(OUT_DIR, exist_ok=True)

# ---------------------------------------------------------------------------
# MVA / spindle-assembly-checkpoint (SAC) gene panel, GRCh38 (no-chr naming)
# ---------------------------------------------------------------------------
GENE_PANEL: dict[str, tuple[str, int, int]] = {
    "BUB1B":  ("15", 40161204, 40221180),
    "CEP57":  ("11", 95577010, 95622499),
    "TRIP13": ("5",  1353562,  1378191),
    "BUB1":   ("2",  110637579, 110678361),
    "BUB3":   ("10", 124922098, 124933422),
    "MAD2L1": ("4",  121076933, 121096285),
    "MAD2L2": ("1",  11826811,  11837222),
    "CDC20":  ("1",  43358981,  43363203),
    "PLK1":   ("16", 23677073,  23690318),
    "AURKB":  ("17", 8207644,   8225940),
    "SGO1":   ("3",  20108421,  20133361),
    "CENPE":  ("4",  103105200, 103198623),
    "NUF2":   ("1",  163266089, 163304940),
}

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def gene_for_pos(chrom: str, pos: int) -> list[str]:
    """Return panel genes overlapping a genomic position."""
    hits = []
    for gene, (c, start, end) in GENE_PANEL.items():
        if chrom == c and start <= pos <= end:
            hits.append(gene)
    return hits


def parse_alt_fraction(ad: tuple[int, ...] | None) -> float | None:
    """Alt allele fraction from FORMAT/AD (ref, alt, ...)."""
    if not ad or len(ad) < 2:
        return None
    total = sum(ad)
    if total == 0:
        return None
    return round(ad[1] / total, 4)


def parse_gt(gt: tuple[int | None, ...] | None) -> str:
    if gt is None:
        return ""
    sep = "|" if "|" in str(gt) else "/"
    parts = [str(g) if g is not None else "." for g in gt]
    return sep.join(parts)