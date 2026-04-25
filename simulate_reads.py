"""Generate fake FASTA + FASTQ for count_between_barcodes.py.

Most reads are clean: [left barcode][insert from a small library][right barcode].
A configurable mix of edge cases (see EDGE_MIX) is injected so that the
counter's drop paths and outermost-barcode rule are all exercised.
The FASTA writer also sprinkles lowercase bases and blank lines to
verify the parser tolerates them.
"""

from __future__ import annotations

import argparse
import os
import random
from collections import Counter

# In FASTA output, every Nth record is written lowercase / followed by a
# blank line, to verify the counter's parser tolerates both.
LOWERCASE_EVERY = 50
BLANK_LINE_EVERY = 100


def rand_dna(rng, n):
    return "".join(rng.choices("ACGT", k=n))


def qual_string(rng, n, mean_q):
    return "".join(chr(max(2, min(40, int(rng.gauss(mean_q, 3)))) + 33)
                   for _ in range(n))


def make_read(rng, kind, library, left, right, insert_len):
    """Return (seq, mean_q_target) for a given edge-case kind."""
    pick = lambda: rng.choice(library)
    if kind == "clean":
        return left + pick() + right, 34
    if kind == "lowq":
        return left + pick() + right, 6
    if kind == "lowq_with_barcodes":
        return left + pick() + right, 5
    if kind == "empty_record":
        return "", 34
    if kind == "empty_insert":
        return left + right, 34
    if kind == "no_left":
        return rand_dna(rng, 5) + pick() + right, 34
    if kind == "no_right":
        return left + pick() + rand_dna(rng, 5), 34
    if kind == "no_barcodes":
        return rand_dna(rng, len(left) + insert_len + len(right)), 34
    if kind == "double_left":
        # Two left barcodes; outer (leftmost) is correct -> insert spans
        # [pick] + left + [pick].
        return left + pick() + left + pick() + right, 34
    if kind == "double_right":
        # Two right barcodes; outer (rightmost) is correct -> insert spans
        # [pick] + right + [pick].
        return left + pick() + right + pick() + right, 34
    if kind == "right_inside":
        # Spurious right barcode inside the insert; rightmost is the real one.
        return left + pick()[:5] + right + pick() + right, 34
    if kind == "long_insert":
        return left + rand_dna(rng, insert_len * 3) + right, 34
    if kind == "short_insert":
        return left + rand_dna(rng, 3) + right, 34
    raise ValueError(kind)


# Fractions sum to <1; remainder is "clean".
EDGE_MIX = {
    "lowq":               0.10,
    "lowq_with_barcodes": 0.05,
    "empty_record":       0.01,
    "empty_insert":       0.01,
    "no_left":            0.03,
    "no_right":           0.03,
    "no_barcodes":        0.02,
    "double_left":        0.02,
    "double_right":       0.02,
    "right_inside":       0.02,
    "long_insert":        0.02,
    "short_insert":       0.02,
}


def simulate(n_reads, library, left, right, insert_len, seed=0):
    rng = random.Random(seed)
    kinds = list(EDGE_MIX) + ["clean"]
    weights = list(EDGE_MIX.values()) + [1 - sum(EDGE_MIX.values())]
    reads = []
    for i in range(n_reads):
        kind = rng.choices(kinds, weights, k=1)[0]
        seq, mq = make_read(rng, kind, library, left, right, insert_len)
        reads.append((f"read_{i}_{kind}", seq, qual_string(rng, len(seq), mq), kind))
    return reads


def write_outputs(reads, out_dir):
    fa_path = os.path.join(out_dir, "reads.fasta")
    fq_path = os.path.join(out_dir, "reads.fastq")
    with open(fa_path, "w") as fa, open(fq_path, "w") as fq:
        for idx, (name, seq, qual, _kind) in enumerate(reads):
            fa_seq = seq.lower() if idx % LOWERCASE_EVERY == 0 else seq
            fa.write(f">{name}\n{fa_seq}\n")
            if idx % BLANK_LINE_EVERY == 0:
                fa.write("\n")
            fq.write(f"@{name}\n{seq}\n+\n{qual}\n")
    return fa_path, fq_path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", default="data")
    ap.add_argument("--n-reads", type=int, default=2000)
    ap.add_argument("--n-variants", type=int, default=50)
    ap.add_argument("--insert-len", type=int, default=21)
    ap.add_argument("--left", default="GGTACCATG")
    ap.add_argument("--right", default="TAATTAAGC")
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    rng = random.Random(args.seed)
    library = list({rand_dna(rng, args.insert_len)
                    for _ in range(args.n_variants * 3)})[:args.n_variants]

    reads = simulate(args.n_reads, library, args.left, args.right,
                     args.insert_len, seed=args.seed)
    write_outputs(reads, args.out_dir)

    print(f"Wrote {len(reads)} reads to {args.out_dir}/reads.fast{{a,q}}")
    for k, v in sorted(Counter(k for *_x, k in reads).items(), key=lambda kv: -kv[1]):
        print(f"  {k:22s} {v}")


if __name__ == "__main__":
    main()
