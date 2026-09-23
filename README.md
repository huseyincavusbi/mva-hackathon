# Rare Disease, Real Kid: MVA Hackathon 2026 — Track 1 (Variant Prediction)

Entry for the [official hackathon](https://sagebio-rare-disease-real-kid-mva-hackathon-2026.hf.space) · [gated dataset](https://huggingface.co/datasets/SageBio/mva-hackathon-2026-data) (Sage Bionetworks × MVA Society × Hugging Face × BEACON).

**What we did:** screened the proband's WGS with a phenotype-guided MVA/spindle-checkpoint gene panel + genome-wide rare-variant, mosaic-VAF and compound-het scan, annotated with Ensembl VEP/ClinVar/gnomAD, and confirmed top calls at read level from raw FASTQ.

**What we found:** compound-heterozygous **BUB1B** loss-of-function — `chr15:40209701 T>G` (c.2210T>G, p.Leu737Ter; ClinVar pathogenic for MVA1) in trans with novel missense `chr15:40220612 T>G` (p.Asn737Lys), both read-verified at ~50% VAF — explaining the child's MVA1 phenotype (rhabdomyosarcoma, IUGR, nephrocalcinosis, parental recurrent miscarriages).

**Code:** pipeline in `scripts/` (run `run_all.sh`); submission CSV in `submission/`. Data must be requested via the gated HF dataset.
