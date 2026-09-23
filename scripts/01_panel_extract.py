"""Phase 1: Extract all variants (PASS AND filtered) in the MVA/SAC gene panel.

Outputs a TSV with genotype, allelic depths, alt fraction, and quality metrics.
Why include filtered calls? MVA is a mosaic condition - the causal allele may
sit at low allele fraction and look "bad" by standard WGS filters, so we must
not throw away LowQual/QD2/MQ40 calls in the disease genes before inspection.

Usage:
    .venv/bin/python 01_panel_extract.py [--out notes/phase1_mva_genes.tsv]
"""
from __future__ import annotations

import argparse
import csv
import sys

import pysam

from common import GENE_PANEL, NOTES_DIR, PROBAND_ID, VCF, parse_gt, parse_alt_fraction

# build chrom -> [(start, end, gene)]
INDEX: dict[str, list[tuple[int, int, str]]] = {}
for gene, (chrom, start, end) in GENE_PANEL.items():
    INDEX.setdefault(chrom, []).append((start, end, gene))

OUT_COLS = [
    "chrom", "pos", "id", "ref", "alt", "qual", "filter",
    "AC", "AF", "AN", "infoDP", "QD", "FS", "MQ", "SOR",
    "dbSNP", "GT", "AD", "fmtDP", "GQ", "alt_frac", "gene",
]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=NOTES_DIR / "phase1_mva_genes.tsv")
    ap.add_argument("--vcf", default=str(VCF))
    args = ap.parse_args()

    rows: list[list[str]] = []
    with pysam.VariantFile(args.vcf) as vcf:
        for rec in vcf.fetch():
            if rec.contig not in INDEX:
                continue
            genes = [g for (s, e, g) in INDEX[rec.contig] if s <= rec.pos <= e]
            if not genes:
                continue
            sample = rec.samples[PROBAND_ID]
            info = rec.info
            db = "DB" if "DB" in info else ""
            ad = sample.get("AD")
            rows.append([
                rec.contig, str(rec.pos), rec.id or ".", rec.ref,
                ",".join(str(a) for a in (rec.alts or [])),
                str(rec.qual) if rec.qual is not None else ".",
                rec.filter.keys()[0] if rec.filter else "PASS",
                str(info.get("AC", "")), str(info.get("AF", "")), str(info.get("AN", "")),
                str(info.get("DP", "")), str(info.get("QD", "")), str(info.get("FS", "")),
                str(info.get("MQ", "")), str(info.get("SOR", "")),
                db, parse_gt(sample["GT"]),
                ",".join(str(x) for x in ad) if ad else "",
                str(sample.get("DP", "")), str(sample.get("GQ", "")),
                str(parse_alt_fraction(ad) or ""),
                ",".join(genes),
            ])

    with open(args.out, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(OUT_COLS)
        w.writerows(rows)

    n = len(rows)
    ndb = sum(1 for r in rows if r[OUT_COLS.index("dbSNP")] == "")
    print(f"panel variants: {n}  (non-dbSNP: {ndb})")
    print(f"written: {args.out}")


if __name__ == "__main__":
    sys.exit(main())