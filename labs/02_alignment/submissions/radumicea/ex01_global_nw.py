#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Ex: Global alignment (Needleman–Wunsch)
Purpose:
  - Load two sequences from the file downloaded in Lab 1 (data/work/<handle>/lab01/).
  - Implement the basic steps of the NW algorithm to obtain a global alignment.

TODO:
Goals:
  - Initialize the scoring matrix for global alignment.
  - Compute cell scores (match, mismatch, gap).

Example run:
  python labs/02_alignment/ex01_global_nw.py --fasta data/work/<handle>/lab01/my_tp53.fa --i1 0 --i2 1
"""

from pathlib import Path
import argparse
from Bio import SeqIO


# ===================== TODO: Scoring matrix initialization =========================================

def init_score_matrix_global(m: int, n: int, gap: int):
    """
    Initialize the (m+1) x (n+1) matrix for global alignment.
    Steps:
      - Create a list of lists filled with zeros (size (m+1)x(n+1)).
      - First column: [i * gap] for i=0..m
      - First row: [j * gap] for j=0..n
    Return the matrix.
    """
    
    score = [[0] * (n + 1) for _ in range(m + 1)]

    for i in range(0, m + 1):
        score[i][0] = i * gap
    for j in range(0, n + 1):
        score[0][j] = j * gap

    return score


def score_cell_global(score, i: int, j: int, a: str, b: str, match: int, mismatch: int, gap: int):
    """
    TODO: Compute the score for one cell (i, j).
    Steps:
      - diagonal = score[i-1][j-1] + (match if a == b else mismatch)
      - up       = score[i-1][j] + gap
      - left     = score[i][j-1] + gap
    Return max(diagonal, up, left).
    """

    diagonal = score[i-1][j-1] + (match if a == b else mismatch)
    up       = score[i-1][j] + gap
    left     = score[i][j-1] + gap

    return max(diagonal, up, left)


def needleman_wunsch(seq1: str, seq2: str, match=1, mismatch=-1, gap=-2):
    # Too big. Killed by OS
    seq1 = seq1[:1000]
    seq2 = seq2[:1000]

    # Simplified implementation of the Needleman–Wunsch algorithm.
    m, n = len(seq1), len(seq2)

    # Initialize score matrix (you will complete this function)
    score = init_score_matrix_global(m, n, gap)

    # Fill the matrix cell by cell
    # Note: loops start from 1 to sequence length
    for i in range(1, m + 1):
        ai = seq1[i - 1]
        for j in range(1, n + 1):
            bj = seq2[j - 1]
            # call the scoring function
            score[i][j] = score_cell_global(score, i, j, ai, bj, match, mismatch, gap)

    # ================== Backtracking ==================
    # start from bottom-right corner (score[m][n])
    align1, align2 = "", ""
    i, j = m, n
    while i > 0 and j > 0:
        current = score[i][j]
        # recompute neighbor scores to decide where we came from
        diag = score[i - 1][j - 1] + (match if seq1[i - 1] == seq2[j - 1] else mismatch)
        up   = score[i - 1][j] + gap
        left = score[i][j - 1] + gap

        if current == diag:  # came from diagonal
            align1 = seq1[i - 1] + align1
            align2 = seq2[j - 1] + align2
            i -= 1; j -= 1
        elif current == up:  # came from above
            align1 = seq1[i - 1] + align1
            align2 = "-" + align2
            i -= 1
        else:                # came from left
            align1 = "-" + align1
            align2 = seq2[j - 1] + align2
            j -= 1

    # if we reached an edge, fill with gaps
    while i > 0:
        align1 = seq1[i - 1] + align1
        align2 = "-" + align2
        i -= 1
    while j > 0:
        align1 = "-" + align1
        align2 = seq2[j - 1] + align2
        j -= 1

    return align1, align2, score[m][n]


def load_two_sequences(fasta_path: Path, i1: int, i2: int):
    """
    Load sequences from a FASTA file and select two by index.
    """
    recs = list(SeqIO.parse(str(fasta_path), "fasta"))
    if len(recs) < 2:
        raise SystemExit("[error] The file must contain at least 2 sequences.")
    if not (0 <= i1 < len(recs) and 0 <= i2 < len(recs)):
        raise SystemExit(f"[error] Invalid indices (0..{len(recs)-1}).")
    return str(recs[i1].seq), str(recs[i2].seq), recs[i1].id, recs[i2].id


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True, help="Path to your FASTA file from data/work/<handle>/lab01/")
    ap.add_argument("--i1", type=int, default=0, help="Index of the first sequence (default 0)")
    ap.add_argument("--i2", type=int, default=1, help="Index of the second sequence (default 1)")
    args = ap.parse_args()

    fasta_path = Path(args.fasta)
    if not fasta_path.exists():
        raise SystemExit(f"[error] File not found: {fasta_path}")

    s1, s2, id1, id2 = load_two_sequences(fasta_path, args.i1, args.i2)
    a1, a2, sc = needleman_wunsch(s1, s2)

    print("=== Global Alignment (Needleman–Wunsch) ===")
    print(f"{id1}  vs  {id2}")
    print(a1)
    print(a2)
    print("Score:", sc)


if __name__ == "__main__":
    main()
