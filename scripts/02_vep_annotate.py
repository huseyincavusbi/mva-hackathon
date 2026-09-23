"""Phase 1b: Annotate variants via Ensembl VEP REST (consequence, gnomAD, ClinVar).

Reads a TSV (from 01_panel_extract.py or any chr/pos/ref/alt table) and queries
the Ensembl VEP REST region endpoint per variant. Produces a merged out TSV with
VEP consequence terms, affected gene/transcript, and ColocatedVariants-derived
gnomAD AF + ClinVar significance.

Ensembl best practice: max ~15 requests/sec, batch when possible. We use the
region endpoint with the allele in the path, properly URL-encoded via requests.

Usage:
    .venv/bin/python 02_vep_annotate.py --in notes/phase1_mva_genes.tsv \
        [--out notes/phase1_mva_genes.vep.tsv] [--rate 5] [--only-nondbsnp]
"""
from __future__ import annotations

import argparse
import csv
import sys
import time
from typing import Any

import requests

from common import NOTES_DIR

ENSEMBL = "https://rest.ensembl.org"
HEADERS = {"Content-Type": "application/json", "Accept": "application/json",
           "User-Agent": "Mozilla/5.0"}


def vep_region(chrom: str, pos: int, ref: str, alt: str, session: requests.Session,
               retries: int = 3) -> dict[str, Any] | None:
    """Annotate one variant. VEP region endpoint expects ONLY the ALT allele.

    Correct format: {chrom}:{start}-{end}/{ALT} where the region spans all REF
    bases (for indels, end = pos + len(ref) - 1). Passing ref/alt yields 400
    ("allele matches reference") for SNPs, so never do that.
    """
    alt = alt.split(",")[0]              # annotate primary alt only
    end = pos + len(ref) - 1              # span all reference bases
    url = f"{ENSEMBL}/vep/human/region/{chrom}:{pos}-{end}/{alt}?"
    for attempt in range(retries):
        try:
            r = session.get(url, headers=HEADERS, timeout=90)
            if r.status_code == 200:
                return r.json()[0]
            if r.status_code in (429, 503):
                time.sleep(2 * (attempt + 1))  # backoff
                continue
            # 400 often means allele doesn't match ref genome at pos; mark it
            return {"_http_error": r.status_code, "input": f"{chrom}:{pos} {ref}>{alt}"}
        except requests.RequestException:
            time.sleep(2 * (attempt + 1))
    return {"_http_error": "timeout", "input": f"{chrom}:{pos} {ref}>{alt}"}


def summarize(vep: dict[str, Any] | None) -> dict[str, str]:
    """Extract consequence/gene/gnomAD/ClinVar summary from a VEP response."""
    if not vep:
        return {}
    if vep.get("_http_error"):
        return {"vep_status": f"error:{vep['_http_error']}"}
    cons: set[str] = set()
    genes: set[str] = set()
    tx: set[str] = set()
    gnomad_afs: list[float] = []
    clinvar: set[str] = set()
    for t in vep.get("transcript_consequences", []):
        cons.update(t.get("consequence_terms", []))
        if t.get("gene_symbol"):
            genes.add(t["gene_symbol"])
        if t.get("transcript_id"):
            tx.add(t["transcript_id"])
        for cv in t.get("colocated_variants", []):
            freqs = cv.get("frequencies", {})
            for src in ("gnomAD", "gnomad"):
                if src in freqs and "AF" in freqs[src]:
                    gnomad_afs.append(float(freqs[src]["AF"]))
            cl = cv.get("clin_sig")
            if cl:
                clinvar.update(cl)
    return {
        "vep_consequence": ",".join(sorted(cons)),
        "vep_gene": ",".join(sorted(genes)),
        "vep_transcript": ",".join(sorted(tx))[:120],
        "gnomad_AF": min(gnomad_afs) if gnomad_afs else "",
        "clinvar": ",".join(sorted(clinvar)),
        "vep_status": "ok",
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="inp", default=str(NOTES_DIR / "phase1_mva_genes.tsv"))
    ap.add_argument("--out", default=None)  # computed below
    ap.add_argument("--rate", type=float, default=5.0, help="requests per second (<=15)")
    ap.add_argument("--only-nondbsnp", action="store_true",
                    help="only annotate variants without dbSNP id (rarer)")
    args = ap.parse_args()

    if not args.out:
        args.out = args.inp.replace(".tsv", ".vep.tsv")

    with open(args.inp) as f:
        reader = csv.DictReader(f, delimiter="\t")
        cols = reader.fieldnames or []
        vari = [dict(r) for r in reader]

    if args.only_nondbsnp:
        vari = [v for v in vari if not v.get("dbSNP")]

    print(f"annotating {len(vari)} variants -> {args.out}", flush=True)
    session = requests.Session()
    out_cols = cols + ["vep_consequence", "vep_gene", "vep_transcript",
                       "gnomad_AF", "clinvar", "vep_status"]
    n_ok = 0
    with open(args.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=out_cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for v in vari:
            vep = vep_region(v["chrom"], int(v["pos"]), v["ref"], v["alt"], session)
            v.update(summarize(vep))
            w.writerow(v)
            if v.get("vep_status") == "ok":
                n_ok += 1
            time.sleep(1.0 / args.rate)
    print(f"done: {n_ok}/{len(vari)} annotated OK")


if __name__ == "__main__":
    sys.exit(main())