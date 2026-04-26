"""Count unique sub-sequences flanked by a left and right barcode in a
FASTA/FASTQ file (Hull et al. 2025 JVI reimplementation).

Outermost-barcode rule: leftmost left barcode + rightmost right barcode
after it, so spurious internal copies are skipped.

Optional Phred quality filters (FASTQ only):
  --min-mean-q F   drop reads whose mean Phred Q is below F
  --min-pos-q  N   drop reads with ANY base Phred Q below N

Diagnostics to stderr: reads in/out, drop reasons, mean Q of kept reads,
insert length min/median/max + ASCII histogram.
"""

from __future__ import annotations

import argparse
import sys
from collections import Counter, namedtuple
from contextlib import nullcontext
from statistics import median


Stats = namedtuple("Stats", "counter diag lengths mean_qs")


def iter_records(path):
    """Yield (sequence, qual_or_None) from FASTA or FASTQ.

    FASTQ uses fixed 4-line records (qual lines may legally start with '@').
    Blank lines and lowercase bases are tolerated.
    """
    with open(path, encoding="utf-8") as fh:
        first = fh.readline()
        while first and not first.strip():
            first = fh.readline()
        if not first:
            return

        if first.startswith(">"):
            parts = []
            for line in fh:
                if line.startswith(">"):
                    yield "".join(parts).upper(), None
                    parts = []
                elif line.strip():
                    parts.append(line.strip())
            yield "".join(parts).upper(), None

        elif first.startswith("@"):
            hdr = first
            while hdr:
                seq = fh.readline().rstrip("\n").upper()
                fh.readline()  # '+'
                qual = fh.readline().rstrip("\n")
                yield seq, qual
                hdr = fh.readline()
                while hdr and not hdr.strip():
                    hdr = fh.readline()
        else:
            raise ValueError(f"Not FASTA/FASTQ: {first!r}")


def quality_ok(qual, min_mean, min_pos):
    """True if qual passes thresholds. None qual (FASTA) always passes."""
    if qual is None or (min_mean is None and min_pos is None):
        return True
    if not qual:
        return False
    qs = [ord(c) - 33 for c in qual]
    if min_mean is not None and sum(qs) / len(qs) < min_mean:
        return False
    return min_pos is None or min(qs) >= min_pos


def count_inserts(path, left, right, min_mean_q=None, min_pos_q=None):
    counter, diag = Counter(), Counter()
    lengths, mean_qs = [], []

    for seq, qual in iter_records(path):
        diag["reads_in"] += 1
        if not quality_ok(qual, min_mean_q, min_pos_q):
            diag["dropped_low_q"] += 1
            continue
        diag["passing_q"] += 1

        i = seq.find(left)
        if i < 0:
            diag["dropped_no_left_bc"] += 1
            continue
        start = i + len(left)
        j = seq.rfind(right, start)
        if j < 0:
            diag["dropped_no_right_bc"] += 1
            continue

        insert = seq[start:j]
        if not insert:
            diag["dropped_empty_insert"] += 1
            continue

        counter[insert] += 1
        lengths.append(len(insert))
        if qual:
            mean_qs.append(sum(ord(c) - 33 for c in qual) / len(qual))

    diag["reads_out"] = sum(counter.values())
    diag["unique_inserts"] = len(counter)
    return Stats(counter, diag, lengths, mean_qs)


def length_histogram(lengths):
    """Return [(length, count)] for every observed insert length, sorted."""
    counts = Counter(lengths)
    return sorted(counts.items())


def write_counts(counter, out_path):
    cm = nullcontext(sys.stdout) if out_path == "-" else open(out_path, "w", encoding="utf-8")
    with cm as out:
        out.write("insert\tcount\n")
        for ins, c in counter.most_common():
            out.write(f"{ins}\t{c}\n")


def write_diagnostics(stats, fh=sys.stderr):
    fh.write("# === diagnostics ===\n")
    for k in ("reads_in", "passing_q", "dropped_low_q",
              "dropped_no_left_bc", "dropped_no_right_bc",
              "dropped_empty_insert", "reads_out", "unique_inserts"):
        fh.write(f"#   {k:24s} {stats.diag.get(k, 0)}\n")

    if stats.mean_qs:
        fh.write(f"#   mean_Q_kept_reads        {sum(stats.mean_qs)/len(stats.mean_qs):.2f}\n")

    if stats.lengths:
        fh.write(f"#   insert_len min/median/max "
                 f"{min(stats.lengths)}/{int(median(stats.lengths))}/{max(stats.lengths)}\n")
        fh.write("#   insert length histogram (per length):\n")
        hist = length_histogram(stats.lengths)
        peak = max(c for _, c in hist) or 1
        for L, c in hist:
            fh.write(f"#     {L:6d} | {c:6d} {'#' * int(40 * c / peak)}\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("input", help="FASTA or FASTQ file")
    ap.add_argument("--left", required=True, help="Left (5') barcode")
    ap.add_argument("--right", required=True, help="Right (3') barcode")
    ap.add_argument("-o", "--out", default="-", help="Output TSV (default stdout)")
    ap.add_argument("--metadata", help="Also write diagnostics to this .txt file")
    ap.add_argument("--min-mean-q", type=float, help="Min mean Phred Q (FASTQ)")
    ap.add_argument("--min-pos-q", type=int, help="Min per-position Phred Q (FASTQ)")
    args = ap.parse_args()

    stats = count_inserts(args.input, args.left.upper(), args.right.upper(),
                          args.min_mean_q, args.min_pos_q)
    write_counts(stats.counter, args.out)
    write_diagnostics(stats)
    if args.metadata:
        with open(args.metadata, "w", encoding="utf-8") as fh:
            fh.write(f"# input        {args.input}\n")
            fh.write(f"# left         {args.left.upper()}\n")
            fh.write(f"# right        {args.right.upper()}\n")
            fh.write(f"# min_mean_q   {args.min_mean_q}\n")
            fh.write(f"# min_pos_q    {args.min_pos_q}\n")
            fh.write(f"# out          {args.out}\n")
            write_diagnostics(stats, fh=fh)


if __name__ == "__main__":
    main()
