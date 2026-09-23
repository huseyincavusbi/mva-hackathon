"""Phase 5: Generate the Track 1 submission CSV in the exact hackathon schema.

Schema (from the live app / notes/02_submission_schema.md):
    proband_id, chrom_1, pos_1, ref_1, alt_1, chrom_2, pos_2, ref_2, alt_2,
    epcr, finding_type, notes

Rules enforced:
  - max 10 rows
  - epcr in (0,1], monotone decreasing with rank
  - GRCh38 coordinates; chrom as "chr15"
  - finding_type in {primary, secondary}; secondary rows don't hurt automated score
  - compound-het pairs use the _2 columns (chrom_2 etc.)

Usage:
    .venv/bin/python 06_make_submission.py --in notes/phase4_ranked.tsv --team NAME
"""
from __future__ import annotations

import argparse
import csv
import sys

from common import NOTES_DIR, PROBAND_ID


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="inp", required=True, help="ranked TSV from 04")
    ap.add_argument("--out", default=NOTES_DIR / "submission_track1.csv")
    ap.add_argument("--team", default="team", help="team/username for filename + notes")
    ap.add_argument("--pair", default=None,
                    help="optional: emit first two rows as ONE compound-het row using "
                         "chrom_2/pos_2/ref_2/alt_2 (schema supports compound-het pairs; "
                         "half-credit if only one allele recovered). Format: chrom2,pos2,ref2,alt2")
    args = ap.parse_args()

    with open(args.inp) as f:
        rows = list(csv.DictReader(f, delimiter="\t"))

    if len(rows) > 11:
        rows = rows[:11]

    pair = None
    if args.pair:
        c2, p2, r2, a2 = args.pair.split(",")
        row1 = rows.pop(0)
        row2 = rows.pop(0)
        pair = (row1, row2, (c2, p2, r2, a2))

    out_cols = ["proband_id", "chrom_1", "pos_1", "ref_1", "alt_1",
                "chrom_2", "pos_2", "ref_2", "alt_2",
                "epcr", "finding_type", "notes"]
    written = 0
    with open(args.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=out_cols)
        w.writeheader()

        def write_row(r, chrom2="", pos2="", ref2="", alt2="", epcr_=None):
            nonlocal written
            try:
                epcr = float(epcr_ if epcr_ is not None else r.get("epcr", 0))
            except ValueError:
                epcr = 0.01
            epcr = max(0.01, min(1.0, round(epcr, 4)))
            gene = r.get("gene", "")
            cons = r.get("vep_consequence", "")
            notes = f"{gene}: {cons}" if gene else cons
            w.writerow({
                "proband_id": PROBAND_ID,
                "chrom_1": f"chr{r['chrom']}" if not r["chrom"].startswith("chr") else r["chrom"],
                "pos_1": r["pos"],
                "ref_1": r["ref"],
                "alt_1": r["alt"].split(",")[0],
                "chrom_2": chrom2,
                "pos_2": pos2,
                "ref_2": ref2,
                "alt_2": alt2,
                "epcr": f"{epcr:.4f}",
                "finding_type": r.get("finding_type", "primary"),
                "notes": notes,
            })
            written += 1

        if pair:
            (r1, r2, (c2, p2, r2a, a2a)) = pair
            write_row(r1, chrom2=f"chr{c2}" if not c2.startswith("chr") else c2,
                      pos2=p2, ref2=r2a, alt2=a2a, epcr_=1.0)
        for r in rows:
            write_row(r)

    print(f"submission written: {args.out} ({written} rows)")
    print(f"  -> rename to {args.team}_<modelname>.csv before uploading")
    with open(args.out) as f:
        print(f.read())


if __name__ == "__main__":
    sys.exit(main())