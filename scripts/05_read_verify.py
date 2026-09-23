"""Phase 4: Read-level verification of top candidate variants from FASTQ.

C-speed approach: for each candidate we fetch the GRCh38 reference window,
build the 21-mer centered on the variant (REF allele) and the corresponding
ALT 21-mer, then count exact occurrences of each (plus reverse complements)
across whole FASTQ files using `gzip -cd | grep -o` subprocesses. Each read
that covers the locus contributes one occurrence of exactly one of the
REF/ALT k-mers, so the ratio is a true read-level estimate.

Usage:
    .venv/bin/python 05_read_verify.py --in notes/phase4_ranked.tsv \
        --fastq data/WGS_EX2312012_HGWCNDSX7_S16_L001_R1_001.fastq.gz [--n 10]
"""
from __future__ import annotations

import argparse
import csv
import subprocess
import sys

import requests

from common import NOTES_DIR

ENSEMBL = "https://rest.ensembl.org"
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
HALF = 10                      # 21-mer, variant base at center


def rc(s: str) -> str:
    return s.translate(COMP)[::-1]


def fetch_window(chrom: str, pos: int, pad: int = 60) -> str | None:
    url = f"{ENSEMBL}/sequence/region/human/{chrom}:{pos - pad}-{pos + pad + 1}:1?"
    r = requests.get(url, headers={"Accept": "application/json",
                                   "User-Agent": "Mozilla/5.0"}, timeout=60)
    if r.status_code != 200:
        return None
    return r.json().get("seq", "")


def count_occ(kmer: str, fastq_path: str) -> int:
    """Exact occurrence count of kmer in a gzipped FASTQ (C-speed)."""
    gz = subprocess.Popen(["gzip", "-cd", fastq_path],
                          stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    gr = subprocess.Popen(["grep", "-o", "-F", "--", kmer],
                          stdin=gz.stdout, stdout=subprocess.PIPE,
                          stderr=subprocess.DEVNULL)
    gz.stdout.close()
    out, _ = gr.communicate()
    gz.wait()
    return out.decode(errors="replace").count("\n")


def total_reads(fastq_path: str) -> int:
    gz = subprocess.Popen(["gzip", "-cd", fastq_path],
                          stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    wc = subprocess.Popen(["wc", "-l"], stdin=gz.stdout,
                          stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    gz.stdout.close()
    out, _ = wc.communicate()
    gz.wait()
    return int(out.decode().strip()) // 4


def verify_one(chrom: str, pos: int, ref: str, alt: str,
               fastq_paths: list[str]) -> dict[str, int | str]:
    win = fetch_window(chrom, pos)
    if not win:
        return {"error": "no ref window"}
    # window is [pos-pad, pos+pad]: variant base at index PAD (=center)
    center = len(win) // 2
    refk = win[center - HALF: center + HALF + 1]
    if refk[HALF] != ref:
        return {"error": f"ref mismatch: window center {refk[HALF]} != {ref}"}
    # Count reads carrying each possible base AT the center position:
    # read supports base b if it contains flank-left + b + flank-right (= refk
    # with center substituted), or its reverse complement. This is robust to
    # the observed position in the read because the 21-mer anchors the site.
    def count_center_base(b: str, fp: str) -> int:
        k = refk[:HALF] + b + refk[HALF + 1:]
        return count_occ(k, fp) + count_occ(rc(k), fp)

    counts: dict[str, int] = {}
    for b in "ACGT":
        counts[b] = sum(count_center_base(b, fp) for fp in fastq_paths)
    refc = counts.get(ref, 0)
    altc = counts.get(alt, 0)
    return {"ref": refc, "alt": altc, "error": "",
            "bases": counts,
            "total_reads": sum(total_reads(fp) for fp in fastq_paths)}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="inp", required=True)
    ap.add_argument("--fastq", action="append", required=True)
    ap.add_argument("--n", type=int, default=10)
    ap.add_argument("--out", default=str(NOTES_DIR / "phase4_readverify.tsv"))
    args = ap.parse_args()

    with open(args.inp) as f:
        rows = list(csv.DictReader(f, delimiter="\t"))[: args.n]

    out_cols = list(rows[0].keys()) + ["total_reads", "reads_ref", "reads_alt",
                                       "obs_alt_frac", "verdict"]
    with open(args.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=out_cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            v = verify_one(r["chrom"], int(r["pos"]), r["ref"], r["alt"].split(",")[0],
                           args.fastq)
            r["total_reads"] = v.get("total_reads", "")
            r["reads_ref"] = v.get("ref", 0)
            r["reads_alt"] = v.get("alt", 0)
            tot = int(v.get("ref", 0)) + int(v.get("alt", 0))
            r["obs_alt_frac"] = round(int(v.get("alt", 0)) / tot, 3) if tot else ""
            r["verdict"] = ("CONFIRMED" if tot and int(v.get("alt", 0)) >= 3 else
                            "MISSING" if tot == 0 else "WEAK")
            w.writerow(r)
            print(f"  {r['chrom']}:{r['pos']} {r['ref']}>{r['alt']} "
                  f"ref={v.get('ref')} alt={v.get('alt')} frac={r['obs_alt_frac']} "
                  f"reads={v.get('total_reads')} -> {r['verdict']}", flush=True)
    print(f"written: {args.out}")


if __name__ == "__main__":
    sys.exit(main())