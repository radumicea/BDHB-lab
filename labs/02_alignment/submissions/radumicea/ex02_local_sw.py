#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Exercise: Local alignment (Smith–Waterman)

Purpose:
  - Load two sequences from the file downloaded in Lab 1 (data/work/<handle>/lab01/).
  - Implement the basic steps of the SW algorithm to obtain a local alignment.

TODO:
  - Initialize the local scoring matrix.
  - Compute cell scores (match, mismatch, gap, with max(0, ...)).

Example run:
  python labs/02_alignment/ex02_local_sw.py --fasta data/work/<handle>/lab01/my_tp53.fa --i1 0 --i2 1
"""

from pathlib import Path
import argparse
from Bio import SeqIO


# ===================== Matrix initialization =========================================

def init_score_matrix_local(m: int, n: int):
    """
    TODO: Initialize the (m+1) x (n+1) matrix with all values = 0.
    Hint: use list comprehension or simple loops.
    """
    return [[0] * (n + 1) for _ in range(m + 1)]


def score_cell_local(score, i: int, j: int, a: str, b: str, match: int, mismatch: int, gap: int):
    """
    TODO: Compute the score for a single cell (i, j).
    Steps:
      - diagonal = score[i-1][j-1] + (match if a == b else mismatch)
      - up       = score[i-1][j] + gap
      - left     = score[i][j-1] + gap
    Result = max(0, diagonal, up, left).
    """
    diagonal = score[i-1][j-1] + (match if a == b else mismatch)
    up       = score[i-1][j] + gap
    left     = score[i][j-1] + gap

    return max(0, diagonal, up, left)


def smith_waterman(seq1: str, seq2: str, match=3, mismatch=-3, gap=-2):
    # Too big. Killed by OS
    seq1 = seq1[:1000]
    seq2 = seq2[:1000]

    # Simplified implementation of the Smith–Waterman algorithm.
    m, n = len(seq1), len(seq2)

    # Matrix initialization
    score = init_score_matrix_local(m, n)

    max_score = 0
    max_pos = (0, 0)

    # Fill the matrix cell by cell
    # + store the highest score and its position
    for i in range(1, m + 1):
        ai = seq1[i - 1]
        for j in range(1, n + 1):
            bj = seq2[j - 1]
            score[i][j] = score_cell_local(score, i, j, ai, bj, match, mismatch, gap)

            if score[i][j] > max_score:
                max_score = score[i][j]
                max_pos = (i, j)

    # ================== Backtracking ==================
    # start from the cell with the maximum score and trace back
    align1, align2 = "", ""
    i, j = max_pos
    while i > 0 and j > 0 and score[i][j] > 0:
        # recompute neighbor scores to decide the direction
        diag = score[i - 1][j - 1] + (match if seq1[i - 1] == seq2[j - 1] else mismatch)
        up   = score[i - 1][j] + gap
        left = score[i][j - 1] + gap

        if score[i][j] == diag:
            align1 = seq1[i - 1] + align1
            align2 = seq2[j - 1] + align2
            i -= 1; j -= 1
        elif score[i][j] == up:
            align1 = seq1[i - 1] + align1
            align2 = "-" + align2
            i -= 1
        else:  # left
            align1 = "-" + align1
            align2 = seq2[j - 1] + align2
            j -= 1

    return align1, align2, max_score


def load_two_sequences(fasta_path: Path, i1: int, i2: int):
    """
    Load the FASTA file and select two sequences by index.
    """
    recs = list(SeqIO.parse(str(fasta_path), "fasta"))
    if len(recs) < 2:
        raise SystemExit("[error] The file must contain at least 2 sequences.")
    if not (0 <= i1 < len(recs) and 0 <= i2 < len(recs)):
        raise SystemExit(f"[error] Invalid indices (0..{len(recs)-1}).")
    return str(recs[i1].seq), str(recs[i2].seq), recs[i1].id, recs[i2].id


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True, help="Path to your FASTA file in data/work/<handle>/lab01/")
    ap.add_argument("--i1", type=int, default=0, help="Index of the first sequence (default 0)")
    ap.add_argument("--i2", type=int, default=1, help="Index of the second sequence (default 1)")
    args = ap.parse_args()

    fasta_path = Path(args.fasta)
    if not fasta_path.exists():
        raise SystemExit(f"[error] File not found: {fasta_path}")

    s1, s2, id1, id2 = load_two_sequences(fasta_path, args.i1, args.i2)
    a1, a2, sc = smith_waterman(s1, s2)

    print("=== Local Alignment (Smith–Waterman) ===")
    print(f"{id1}  vs  {id2}")
    print(a1)
    print(a2)
    print("Score:", sc)


if __name__ == "__main__":
    main()
