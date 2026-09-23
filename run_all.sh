#!/usr/bin/env bash
# Run the full Track 1 pipeline end-to-end.
set -euo pipefail
cd "$(dirname "$0")"
PY=scripts/.venv/bin/python
[ -x "$PY" ] || { echo "venv missing — run: uv venv scripts/.venv && uv pip install pysam pandas numpy requests"; exit 1; }

"$PY" scripts/01_panel_extract.py --out notes/phase1_mva_genes.tsv
"$PY" scripts/02_vep_annotate.py --in notes/phase1_mva_genes.tsv --out notes/phase1_mva_genes.vep.tsv --rate 5
"$PY" scripts/03_genomewide_scan.py
"$PY" scripts/04_score_rank.py --in notes/phase1_mva_genes.vep.tsv --out notes/phase4_ranked.tsv
"$PY" scripts/05_read_verify.py --in notes/phase4_ranked.tsv --fastq data/WGS_EX2312012_HGWCNDSX7_S16_L001_R1_001.fastq.gz --out notes/phase4_readverify.tsv
"$PY" scripts/06_make_submission.py --in notes/phase4_ranked.tsv --out notes/submission_track1.csv --team myteam
echo "DONE — see notes/submission_track1.csv"