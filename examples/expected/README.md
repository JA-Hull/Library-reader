# Expected outputs

These TSVs are the result of running the demo with the default seed:

```bash
./run.sh
```

(barcodes `GGTACCATG` / `TAATTAAGC`, 2000 reads, seed 0).

Files:

- `counts_fasta.tsv`        — FASTA input, no Q filter
- `counts_fasta.meta.txt`   — diagnostics for the above
- `counts_meanQ20.tsv`      — FASTQ input, `--min-mean-q 20`
- `counts_meanQ20.meta.txt` — diagnostics for the above
- `counts_posQ15.tsv`       — FASTQ input, `--min-pos-q 15`
- `counts_posQ15.meta.txt`  — diagnostics for the above

To regenerate:

```bash
rm -rf data out
./run.sh
cp out/*.tsv out/*.meta.txt examples/expected/
```
