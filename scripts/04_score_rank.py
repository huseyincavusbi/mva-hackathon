"""Phase 3: Score & rank candidates into an EPCR-style ranked list.

Merges evidence per variant:
  - phenotype fit (gene in MVA/SAC panel or HPO-prioritized gene list)  [0..1]
  - consequence severity (LOF > missense > splice-adjacent > intronic)  [0..1]
  - rarity (gnomAD AF; lower = more likely causal)                      [0..1]
  - genotype/read evidence (alt fraction in het/mosaic band, depth)     [0..1]
  - in-silico pathogenicity flags (AlphaMissense/CADD/REVEL where known) [0..1]

Final score = weighted combination; EPCR (estimated probability of causal
relationship) is the normalized score in (0,1]. Rows are then ranked desc.

Input: annotated TSV(s) from 02/03. Output: notes/phase4_ranked.tsv (sorted,
with epcr + finding_type + notes columns matching submission schema fields).

Usage:
    .venv/bin/python 04_score_rank.py --in notes/phase2_candidates.tsv
"""
from __future__ import annotations

import argparse
import csv
import sys

from common import GENE_PANEL, NOTES_DIR

LOF_TERMS = {"stop_gained", "frameshift_variant", "splice_acceptor_variant",
             "splice_donor_variant", "start_lost", "stop_lost", "exon_loss_variant"}
MISSENSE_TERMS = {"missense_variant", "inframe_deletion", "inframe_insertion",
                   "protein_altering_variant", "splice_region_variant"}
SYNONYMOUS_TERMS = {"synonymous_variant"}
FUNCTIONAL_TERMS = LOF_TERMS | MISSENSE_TERMS


def is_functional(cons: str) -> bool:
    """Functional gate: only protein-altering or splice-region consequences
    can be causal candidates. Pure intronic/ncRNA/downstream/upstream calls
    (e.g. repeat artifacts) are excluded before scoring so they can no
    longer inflate the ranked list."""
    terms = {t.strip() for t in cons.split(",") if t.strip()}
    return bool(terms & FUNCTIONAL_TERMS)

# Weights for evidence axes (heuristic; calibrate after leaderboard feedback)
W = {"phenotype": 0.30, "consequence": 0.25, "rarity": 0.20,
     "reads": 0.15, "insilico": 0.10}


def consequence_score(cons: str) -> float:
    terms = {t.strip() for t in cons.split(",") if t.strip()}
    if any(t in terms for t in LOF_TERMS):
        return 1.0
    if any(t in terms for t in MISSENSE_TERMS):
        return 0.7
    if any(t in terms for t in SYNONYMOUS_TERMS):
        return 0.1
    return 0.2  # intronic/ncRNA default


def rarity_score(gnomad_af: str) -> float:
    if not gnomad_af:
        return 0.5  # unknown, not evidence of rarity (blank means
                    # unannotated, not absent from gnomAD)
    try:
        af = float(gnomad_af)
    except ValueError:
        return 0.5
    if af == 0:
        return 1.0
    if af < 1e-4:
        return 0.85
    if af < 1e-3:
        return 0.6
    if af < 0.01:
        return 0.3
    return 0.05


def phenotype_score(gene: str) -> float:
    if gene in GENE_PANEL:
        return 1.0  # MVA/SAC panel gene = maximal phenotype fit
    if gene in {"CENPE", "BUB1", "BUB3", "MAD2L1", "MAD2L2", "CDC20", "PLK1",
                "AURKB", "SGO1", "NUF2"}:
        return 0.9
    return 0.1


def reads_score(af: str, gq: str, dp: str) -> float:
    try:
        afv = float(af)
    except (TypeError, ValueError):
        return 0.1
    # clean het or clear hom, and decent depth
    if 0.3 <= afv <= 0.7 or afv >= 0.9:
        base = 0.8
    elif afv >= 0.2:
        base = 0.5
    else:
        base = 0.2
    try:
        dpv = int(dp)
        if dpv >= 20:
            base = min(1.0, base + 0.15)
        elif dpv < 8:
            base *= 0.5
    except (TypeError, ValueError):
        pass
    return base


def insilico_score(gnomad_af: str, cons: str, clinvar: str) -> float:
    # no full AlphaMissense/CADD columns yet -> rely on consequence + ClinVar
    s = consequence_score(cons) * 0.6
    if clinvar and any(k in clinvar.upper() for k in ("PATHOGENIC", "LIKELY_PATHOGENIC")):
        s = max(s, 0.9)
    return s


def score_variant(row: dict[str, str]) -> float:
    s = (
        W["phenotype"] * phenotype_score(row.get("gene", ""))
        + W["consequence"] * consequence_score(row.get("vep_consequence", ""))
        + W["rarity"] * rarity_score(row.get("gnomad_AF", ""))
        + W["reads"] * reads_score(row.get("alt_frac", ""), row.get("GQ", ""), row.get("fmtDP", ""))
        + W["insilico"] * insilico_score(row.get("gnomad_AF", ""),
                                         row.get("vep_consequence", ""),
                                         row.get("clinvar", ""))
    )
    return round(s, 4)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="inp", required=True, help="annotated candidate TSV")
    ap.add_argument("--out", default=NOTES_DIR / "phase4_ranked.tsv")
    ap.add_argument("--top", type=int, default=20, help="rows to keep in ranked list")
    args = ap.parse_args()

    with open(args.inp) as f:
        reader = csv.DictReader(f, delimiter="\t")
        rows = [dict(r) for r in reader]

    rows = [r for r in rows if is_functional(r.get("vep_consequence", ""))]

    for r in rows:
        r["score"] = score_variant(r)

    # stable sort: score desc, then phenotype/rarity as tiebreak
    rows.sort(key=lambda r: (r["score"],
                             phenotype_score(r.get("gene", "")),
                             rarity_score(r.get("gnomad_AF", ""))),
              reverse=True)
    rows = rows[: args.top]

    # normalize score to EPCR in (0,1]
    max_s = max(r["score"] for r in rows) if rows else 1.0
    for i, r in enumerate(rows, start=1):
        r["rank"] = i
        epcr = round(max(0.01, r["score"] / max_s), 4) if max_s > 0 else 0.01
        r["epcr"] = epcr
        r["finding_type"] = "primary" if i <= 10 else "secondary"

    cols = ["rank", "chrom", "pos", "ref", "alt", "gene", "vep_consequence",
            "gnomad_AF", "alt_frac", "score", "epcr", "finding_type", "clinvar", "notes"]
    with open(args.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)

    print(f"ranked {len(rows)} rows -> {args.out}")
    for r in rows[:10]:
        print(f"  #{r['rank']} {r['chrom']}:{r['pos']} {r['ref']}>{r['alt']} "
              f"gene={r.get('gene','-')} cons={r.get('vep_consequence','-')[:30]} epcr={r['epcr']}")


if __name__ == "__main__":
    sys.exit(main())