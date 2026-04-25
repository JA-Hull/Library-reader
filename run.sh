#!/usr/bin/env bash
# Demo: simulate reads, then count unique inserts with and without Q filter.
set -euo pipefail
cd "$(dirname "$0")"

LEFT=GGTACCATG
RIGHT=TAATTAAGC
mkdir -p out

python3 simulate_reads.py --out-dir data --left "$LEFT" --right "$RIGHT"

echo "--- FASTA, no Q ---"
python3 count_between_barcodes.py data/reads.fasta \
    --left "$LEFT" --right "$RIGHT" \
    --out      out/counts_fasta.tsv \
    --metadata out/counts_fasta.meta.txt

echo "--- FASTQ, mean Q >= 20 ---"
python3 count_between_barcodes.py data/reads.fastq \
    --left "$LEFT" --right "$RIGHT" --min-mean-q 20 \
    --out      out/counts_meanQ20.tsv \
    --metadata out/counts_meanQ20.meta.txt

echo "--- FASTQ, per-position Q >= 15 ---"
python3 count_between_barcodes.py data/reads.fastq \
    --left "$LEFT" --right "$RIGHT" --min-pos-q 15 \
    --out      out/counts_posQ15.tsv \
    --metadata out/counts_posQ15.meta.txt
