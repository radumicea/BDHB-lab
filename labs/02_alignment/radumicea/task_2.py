#!/usr/bin/env python3
from pathlib import Path
import argparse
from Bio import SeqIO
from Bio.Align import PairwiseAligner


def load_two_sequences(fasta_path: Path, i1: int, i2: int):
    recs = list(SeqIO.parse(str(fasta_path), "fasta"))
    if len(recs) < 2:
        raise SystemExit("[error] The file must contain at least 2 sequences.")
    if not (0 <= i1 < len(recs) and 0 <= i2 < len(recs)):
        raise SystemExit(f"[error] Invalid indices (0..{len(recs)-1}).")
    return str(recs[i1].seq), str(recs[i2].seq), recs[i1].id, recs[i2].id


def run_alignment(mode, s1, s2):
    aligner = PairwiseAligner()
    aligner.mode = mode

    aligner.match_score = 1
    aligner.mismatch_score = 0
    aligner.open_gap_score = -1
    aligner.extend_gap_score = -0.5

    alns = aligner.align(s1, s2)
    return alns[0]


def summarize_alignment(label, aln, s1, s2):
    blocks1, blocks2 = aln.aligned
    start1, end1 = blocks1[0][0], blocks1[-1][1]
    start2, end2 = blocks2[0][0], blocks2[-1][1]

    aligned_len = sum(e - s for (s, e) in blocks1)

    g1, g2 = str(aln[0]), str(aln[1])
    gaps1 = g1.count("-")
    gaps2 = g2.count("-")

    aln_len = len(g1)

    print(f"{label} alignment")
    print(f"  Score: {aln.score}")
    print(f"  Alignment columns (with gaps): {aln_len}")
    print(f"  Seq1 aligned region: {start1}-{end1} (len {end1 - start1}) of {len(s1)}")
    print(f"  Seq2 aligned region: {start2}-{end2} (len {end2 - start2}) of {len(s2)}")
    print(f"  Total aligned (matches+subs, ungapped): {aligned_len}")
    print(f"  Gaps in seq1 (row 0): {gaps1}")
    print(f"  Gaps in seq2 (row 1): {gaps2}")


def find_internal_gap_center(g1, g2):
    """
    Find an internal position where at least one row has a gap ('-').
    Avoids very ends. If none exists, return center.
    """
    n = min(len(g1), len(g2))
    left = 5
    right = max(left, n - 5)

    for i in range(left, right):
        if g1[i] == "-" or g2[i] == "-":
            return i
    return n // 2


def print_alignment_fragment(label, aln, width=40):
    """
    Print a small fragment (~width columns) of the alignment.
    We try to center around an internal gap; if none, show a middle chunk.
    """
    g1, g2 = str(aln[0]), str(aln[1])
    n = len(g1)

    center = find_internal_gap_center(g1, g2)
    start = max(0, center - width // 2)
    end = min(n, start + width)

    print(f"\n{label} alignment fragment (columns {start}:{end}):")
    print(aln[:, start:end])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--i1", type=int, default=0)
    ap.add_argument("--i2", type=int, default=1)
    args = ap.parse_args()

    s1, s2, id1, id2 = load_two_sequences(Path(args.fasta), args.i1, args.i2)

    global_aln = run_alignment("global", s1, s2)
    local_aln = run_alignment("local", s1, s2)

    print(f"Sequences: {id1} vs {id2}\n")

    summarize_alignment("Global", global_aln, s1, s2)
    summarize_alignment("Local",  local_aln,  s1, s2)

    print_alignment_fragment("Global", global_aln, width=60)
    print_alignment_fragment("Local",  local_aln,  width=60)

    print("\nComparison:")
    print(f"  Global score: {global_aln.score}")
    print(f"  Local  score: {local_aln.score}")


if __name__ == "__main__":
    main()
