"""Phase 2b: Genome-wide compound-heterozygosity pairs using the full gene BED.

MVA1 (and MVA in general) is autosomal recessive; the causal configuration is
very often TWO rare variants in the SAME gene, in trans. The panel screen (03)
only pairs panel genes; this extends pairing to all 61k Gencode genes using a
sorted-interval lookup over notes/gencode.v44.genes.bed.

Input: phase2_candidates.tsv (non-dbSNP, PASS, good depth het/hom variants).
Output: notes/phase3_compound_het_genomewide.tsv — one row per gene with >=2
rare alleles, listing all pairs (for manual/clinical triage).

Usage:
    .venv/bin/python 03b_ch_genomewide.py [--bed notes/gencode.v44.genes.bed]
"""
from __future__ import annotations

import argparse
import bisect
import csv
import sys
from collections import defaultdict

from common import NOTES_DIR

HEMIMIN, HEMIMAX = 0.20, 0.80   # germline het band
HOM_MIN = 0.85                  # near-homozygous alt


def load_intervals(path: str) -> dict[str, list]:
    """chrom -> sorted list of (start, end, gene)."""
    per_chrom: dict[str, list] = defaultdict(list)
    with open(path) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 4:
                continue
            chrom = p[0]
            if chrom.startswith("chr"):
                chrom = chrom[3:]
            per_chrom[chrom].append((int(p[1]), int(p[2]), p[3]))
    for chrom in per_chrom:
        per_chrom[chrom].sort(key=lambda x: x[0])
    return per_chrom


def genes_at(intervals: dict[str, list], chrom: str, pos: int) -> list[str]:
    """Genes overlapping pos, via bisect on sorted starts."""
    lst = intervals.get(chrom)
    if not lst:
        return []
    starts = [x[0] for x in lst]
    i = bisect.bisect_right(starts, pos) - 1
    out = []
    # scan backwards for any interval containing pos
    j = i
    while j >= 0 and lst[j][1] >= pos:
        out.append(lst[j][2])
        j -= 1
    # also scan forward for intervals starting at same pos
    k = i + 1
    while k < len(lst) and lst[k][0] <= pos <= lst[k][1]:
        out.append(lst[k][2])
        k += 1
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="inp", default=str(NOTES_DIR / "phase2_candidates.tsv"))
    ap.add_argument("--bed", default=str(NOTES_DIR / "gencode.v44.genes.bed"))
    ap.add_argument("--out", default=str(NOTES_DIR / "phase3_compound_het_genomewide.tsv"))
    args = ap.parse_args()

    print("loading gene intervals...", flush=True)
    intervals = load_intervals(args.bed)

    # chrom,pos -> candidate row
    with open(args.inp) as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    print(f"{len(rows)} candidates", flush=True)

    # assign genes; keep only het/hom with alt_frac in plausible diploid bands
    by_gene: defaultdict[str, list[dict]] = defaultdict(list)
    for r in rows:
        try:
            af = float(r["alt_frac"])
        except (TypeError, ValueError):
            continue
        if not (HEMIMIN <= af <= HEMIMAX or af >= HOM_MIN):
            continue
        genes = genes_at(intervals, r["chrom"], int(r["pos"]))
        for g in genes:
            by_gene[g].append(r)

    pairs = []
    for gene, vs in by_gene.items():
        if len(vs) < 2:
            continue
        vs.sort(key=lambda r: (r["chrom"], int(r["pos"])))
        for i in range(len(vs)):
            for j in range(i + 1, len(vs)):
                a, b = vs[i], vs[j]
                # panel-related or not: report all; triage later
                pairs.append([
                    gene,
                    a["chrom"], a["pos"], a["ref"], a["alt"], a["alt_frac"], a["GQ"],
                    b["chrom"], b["pos"], b["ref"], b["alt"], b["alt_frac"], b["GQ"],
                ])

    with open(args.out, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["gene", "chrom_a", "pos_a", "ref_a", "alt_a", "af_a", "gq_a",
                    "chrom_b", "pos_b", "ref_b", "alt_b", "af_b", "gq_b"])
        w.writerows(pairs)

    print(f"genes with >=2 rare alleles: {sum(1 for g, v in by_gene.items() if len(v) >= 2)}")
    print(f"total pairs: {len(pairs)} -> {args.out}")

    # highlight MVA-relevant genes if present
    mvagenes = {"BUB1B", "CEP57", "TRIP13", "BUB1", "BUB3", "MAD2L1", "MAD2L2",
                "CDC20", "PLK1", "AURKB", "SGO1", "CENPE", "NUF2"}
    for p in pairs:
        if p[0] in mvagenes:
            print("  ** MVA-gene CH pair:", p)


if __name__ == "__main__":
    sys.exit(main())