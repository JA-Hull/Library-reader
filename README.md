# barcode_counter

A small reimplementation of the custom Python script described by **Hull
et al. (2025, JVI)** for tallying unique sub-sequences flanked by a
fixed left and right barcode/adapter in FASTA or FASTQ reads.

## Background

The original script was written by me (**Joshua Alexander Hull**) to
prepare **count files of viral libraries** so that fold-changes between
different FASTQs could be computed downstream — count files are far
more memory-efficient to handle for FC analysis than re-parsing raw
reads each time.

I no longer hold the rights to that original code. **This repository is
a clean-room reimplementation written from the ground up, fully
independent of the original script** as a test of Claude Opus 4.6's capabilities. It is provided as-is under the
**MIT License** (see [LICENSE](LICENSE)), with no encumbrance from the
earlier proprietary version.

## Citations

If you use this tool, please cite the original methodological context:

- **Hull, J. A.**, Fusco, R. M., Tan, J., Ochoa, M. A., Hall, A., Wan,
  X., Loeb, E., & Asokan, A. (2025).
  *Functional orthogonality of parvoviral phospholipase A2 domains in
  adeno-associated virus transduction.* **J. Virol.** 99(9):e00799-25.
  [doi:10.1128/jvi.00799-25](https://doi.org/10.1128/jvi.00799-25)
  (PMID: 40793780, PMCID: PMC12455919)

- Rosales, A., Blondel, L. O., **Hull, J.**, Gao, Q., Aykun, N., Peek,
  J. L., Vargas, A., Fergione, S., Song, M., Wilson, M. H., Barbas,
  A. S., & Asokan, A. (2025). *Evolving adeno-associated viruses for
  gene transfer to the kidney via cross-species cycling of capsid
  libraries.* **Nature Biomedical Engineering** 9:1086–1100.
  [doi:10.1038/s41551-024-01341-0](https://doi.org/10.1038/s41551-024-01341-0)

- Wolfson, D. W., **Hull, J. A.**, Li, Y., Gonzalez, T. J., Jayaram,
  M. D., Devlin, G. W., Cigliola, V., Oonk, K. A., Rosales, A., Poss,
  K. D., et al. (2025). *Spatial and longitudinal tracking of
  enhancer-AAV vectors that target transgene expression to injured
  mouse myocardium.* **eLife** (reviewed preprint).
  [doi:10.7554/eLife.107148](https://doi.org/10.7554/eLife.107148)

## Maintainer

Joshua Alexander Hull — [LinkedIn](https://www.linkedin.com/in/ja-hull/)

## Features

- **Standard library only** — no `pip install` required (Python 3.9+).
- Accepts **FASTA or FASTQ**, auto-detected from the first record.
- **Outermost-barcode rule:** leftmost left + rightmost right (after the
  left), so spurious internal barcode copies are skipped.
- Optional Phred quality filters (FASTQ): **mean-Q** or **per-position
  Q**.
- Useful **diagnostics on stderr**: read counts in/out, drop reasons,
  mean Q of kept reads, insert length stats and ASCII histogram.
- A companion **synthetic data generator** that injects a wide variety
  of edge cases (empty records, missing or duplicated barcodes, very
  low quality, etc.) to exercise the counter.

## Layout

```
barcode_counter/
├── count_between_barcodes.py   # the tool
├── simulate_reads.py           # synthetic FASTA + FASTQ generator
├── run.sh                      # end-to-end demo
├── examples/expected/          # committed golden outputs (see below)
├── README.md
├── LICENSE                     # MIT
└── .gitignore
```

## Quick start

```bash
git clone <url> && cd Library-reader
./run.sh
```

That will:

1. Generate `data/reads.fasta` and `data/reads.fastq` (2000 reads each).
2. Run the counter three ways and write TSVs to `out/`.
3. Print diagnostics to the terminal.

## Standalone usage

```bash
python3 count_between_barcodes.py reads.fastq \
    --left  GGTACCATG \
    --right TAATTAAGC \
    --min-mean-q 20 \
    --out counts.tsv
```

`counts.tsv` is sorted by descending count:

```
insert	count
ACGT...	123
...
```

Diagnostics go to stderr:

```
# === diagnostics ===
#   reads_in                 2000
#   passing_q                1687
#   dropped_low_q            313
#   dropped_no_left_bc       100
#   dropped_no_right_bc      64
#   dropped_empty_insert     19
#   reads_out                1504
#   unique_inserts           50
#   mean_Q_kept_reads        33.49
#   insert_len min/median/max 3/21/63
#   insert length histogram (per length):
#          3 |     36 #
#         21 |   1313 ########################################
#         35 |     33 #
#         51 |     67 ##
#         63 |     40 #
```

## CLI options

### `count_between_barcodes.py`

| flag             | meaning                                                  |
| ---------------- | -------------------------------------------------------- |
| `input`          | FASTA or FASTQ file (required, positional)               |
| `--left STR`     | Left (5') barcode (required)                             |
| `--right STR`    | Right (3') barcode (required)                            |
| `-o`, `--out`    | Output TSV path; `-` (default) means stdout              |
| `--metadata PATH`| Also write the diagnostics block to this `.txt` file     |
| `--min-mean-q F` | Drop reads whose mean Phred Q is below `F` (FASTQ only)  |
| `--min-pos-q N`  | Drop reads with ANY base Phred Q below `N` (FASTQ only)  |

### `simulate_reads.py`

| flag             | default        | meaning                                |
| ---------------- | -------------- | -------------------------------------- |
| `--out-dir`      | `data`         | Where to write `reads.fasta/.fastq`    |
| `--n-reads`      | `2000`         | Total reads                            |
| `--n-variants`   | `50`           | Library size                           |
| `--insert-len`   | `21`           | Length of each library variant         |
| `--left`         | `GGTACCATG`    | Left barcode embedded in clean reads   |
| `--right`        | `TAATTAAGC`    | Right barcode embedded in clean reads  |
| `--seed`         | `0`            | RNG seed                               |

The simulator's edge-case mix is hardcoded in `EDGE_MIX` at the top of
`simulate_reads.py` and currently covers:

- `lowq`, `lowq_with_barcodes` — exercises the Q filters.
- `empty_record`, `empty_insert` — exercises empty-handling.
- `no_left`, `no_right`, `no_barcodes` — drop-reason coverage.
- `double_left`, `double_right`, `right_inside` — outermost-barcode rule.
- `long_insert`, `short_insert` — populates the length histogram tails.

## Expected output

Golden TSVs from a fresh `./run.sh` (with `--seed 0`) are committed
under [`examples/expected/`](examples/expected/) so you can diff against
your own run.

## Authorship

This repository was generated by **Claude Opus 4** (Anthropic), driven
by an interactive pair-programming session in the Cursor CLI. All code,
documentation, and example outputs in this repo were produced by that
model.

## License

MIT — see [LICENSE](LICENSE).
