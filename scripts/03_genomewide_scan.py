"""Phase 2: Genome-wide rare functional candidate scan + compound-het + mosaic.

Streams the VCF once and keeps variants that look potentially causal in a
single affected proband with an autosomal-recessive / mosaic disease:

  * PASS (or a small allowlist of informative filtered tags) AND
  * not in dbSNP (and where stage allows, not common) AND
  * heterozygous/homozygous with reasonable depth, AND
  * allele fraction in a range consistent with germline-het (0.2-0.8)
    OR flagged separately as a putative LOW-VAF MOSAIC call (0.04-0.30).

It also detects compound-heterozygous *pairs*: two rare variants in the same
gene. Gene membership for CH detection uses the panel + a cached gene-BED if
provided (--gene-bed).

Outputs:
  notes/phase2_candidates.tsv       - candidate variants
  notes/phase2_compound_het.tsv     - candidate gene-level CH pairs
  notes/phase2_mosaic.tsv           - low-VAF mosaic candidates

Usage:
    .venv/bin/python 03_genomewide_scan.py [--min-gq 30] [--min-dp 8]
"""
from __future__ import annotations

import argparse
import csv
import sys
from collections import defaultdict

import pysam

from common import NOTES_DIR, PROBAND_ID, VCF, parse_gt, parse_alt_fraction

# genes WITHOUT an explicit panel entry are ignored for CH pairing unless a
# gene BED is provided; panel genes are always considered.
PANEL_GENES = {
    "BUB1B", "CEP57", "TRIP13", "BUB1", "BUB3", "MAD2L1", "MAD2L2",
    "CDC20", "PLK1", "AURKB", "SGO1", "CENPE", "NUF2",
}

MOSAIC_MIN, MOSAIC_MAX = 0.04, 0.30   # low-VAF band
HET_MIN, HET_MAX = 0.20, 0.80         # germline het band


def load_gene_bed(path: str | None) -> dict[tuple[str, int], str]:
    """Load chr<TAB>start<TAB>end<TAB>gene -> map pos->gene (any overlap)."""
    out: dict[tuple[str, int], str] = {}
    if not path:
        return out
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 4:
                continue
            chrom, start, end, gene = p[0], int(p[1]), int(p[2]), p[3]
            for pos in range(start, end + 1):
                out[(chrom, pos)] = gene
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf", default=str(VCF))
    ap.add_argument("--gene-bed", default=None,
                    help="chr start end gene BED to map genome-wide variants to genes (for CH)")
    ap.add_argument("--min-gq", type=int, default=30)
    ap.add_argument("--min-dp", type=int, default=8)
    ap.add_argument("--max-mosaic-dp", type=int, default=12,
                    help="mosaic band only kept if DP <= this (low-VAF at deep DP is often noise)")
    ap.add_argument("--out-dir", default=NOTES_DIR)
    args = ap.parse_args()

    gene_bed = load_gene_bed(args.gene_bed)
    cand_cols = ["chrom", "pos", "id", "ref", "alt", "qual", "filter",
                 "GT", "AD", "fmtDP", "GQ", "alt_frac", "gene", "note"]
    mosaic_cols = cand_cols
    ch_cols = ["gene", "chrom_1", "pos_1", "alt_1", "altfrac_1",
               "chrom_2", "pos_2", "alt_2", "altfrac_2"]

    candidates: list[list[str]] = []
    mosaic: list[list[str]] = []
    rare_by_gene: defaultdict[str, list[tuple]] = defaultdict(list)

    with pysam.VariantFile(args.vcf) as vcf:
        for rec in vcf.fetch():
            if rec.contig.startswith(("GL", "KI", "KN", "Un_")):
                continue
            sample = rec.samples[PROBAND_ID]
            gq = sample.get("GQ") or 0
            dp = sample.get("DP") or 0
            ad = sample.get("AD")
            af = parse_alt_fraction(ad)
            chrom, pos = rec.contig, rec.pos
            gene = (gene_bed or {}).get((chrom, pos), "")
            is_panel = False
            # panel membership from known coords (overrides BED when present)
            if not gene:
                from common import GENE_PANEL
                for gname, (c, s, e) in GENE_PANEL.items():
                    if chrom == c and s <= pos <= e:
                        gene, is_panel = gname, True
                        break

            in_db = "DB" in rec.info
            if in_db:
                continue

            # quality gate
            if gq < args.min_gq and not is_panel:
                continue
            if dp < args.min_dp:
                continue
            if af is None:
                continue

            gt = parse_gt(sample["GT"])
            alt = ",".join(str(a) for a in (rec.alts or []))
            row = [chrom, str(pos), rec.id or ".", rec.ref, alt,
                   str(rec.qual if rec.qual is not None else "."),
                   rec.filter.keys()[0] if rec.filter else "PASS",
                   gt, ",".join(str(x) for x in ad) if ad else "",
                   str(dp), str(gq), str(af), gene, ""]

            if MOSAIC_MIN <= af <= MOSAIC_MAX and dp <= args.max_mosaic_dp:
                row[-1] = "low-VAF mosaic"
                mosaic.append(row)
            elif HET_MIN <= af <= HET_MAX or (af > 0.8 and rec.alts and len(rec.alts) == 1):
                candidates.append(row)
                if gene:
                    rare_by_gene[gene].append((chrom, pos, alt, af))

    # compound-heterozygous pairs: two rare het variants in same gene
    ch_pairs: list[list[str]] = []
    for gene, variants in rare_by_gene.items():
        if len(variants) < 2:
            continue
        variants.sort(key=lambda x: x[1])
        for i in range(len(variants)):
            for j in range(i + 1, len(variants)):
                v1, v2 = variants[i], variants[j]
                # require both het (alt fraction in het band or mosaic band)
                if (HET_MIN <= v1[3] <= HET_MAX or MOSAIC_MIN <= v1[3] <= MOSAIC_MAX) and \
                   (HET_MIN <= v2[3] <= HET_MAX or MOSAIC_MIN <= v2[3] <= MOSAIC_MAX):
                    ch_pairs.append([gene, v1[0], str(v1[1]), v1[2], str(v1[3]),
                                     v2[0], str(v2[1]), v2[2], str(v2[3])])

    def write_tsv(path, header, rows):
        with open(path, "w", newline="") as f:
            w = csv.writer(f, delimiter="\t")
            w.writerow(header)
            w.writerows(rows)

    write_tsv(args.out_dir / "phase2_candidates.tsv", cand_cols, candidates)
    write_tsv(args.out_dir / "phase2_mosaic.tsv", mosaic_cols, mosaic)
    write_tsv(args.out_dir / "phase2_compound_het.tsv", ch_cols, ch_pairs)

    print(f"candidates (PASS-ish, non-dbSNP, good depth): {len(candidates)}")
    print(f"low-VAF mosaic candidates: {len(mosaic)}")
    print(f"compound-het gene pairs: {len(ch_pairs)}")
    print(f"outputs in {args.out_dir}")


if __name__ == "__main__":
    sys.exit(main())