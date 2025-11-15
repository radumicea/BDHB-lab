#!/usr/bin/env python
import argparse
from itertools import combinations
from Bio import SeqIO


def hamming_distance(s1: str, s2: str) -> int:
    if len(s1) != len(s2):
        raise ValueError("Hamming distance requires equal-length sequences")
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))


def p_distance(s1: str, s2: str) -> float:
    L = min(len(s1), len(s2))
    s1 = s1[:L]
    s2 = s2[:L]
    diffs = hamming_distance(s1, s2)
    return diffs / float(L) if L > 0 else 0.0


def compute_distance_matrix(recs, use_hamming=False):
    n = len(recs)
    seqs = [str(r.seq) for r in recs]
    ids = [rec.id for rec in recs]
    mat = [[0.0] * n for _ in range(n)]

    for (i, j) in combinations(range(len(seqs)), 2):
        s1, s2 = seqs[i], seqs[j]
        if use_hamming:
            d = hamming_distance(s1, s2)
        else:
            d = p_distance(s1, s2)
        mat[i][j] = d
        mat[j][i] = d

    return ids, mat


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True, help="path to the FASTA file")
    args = ap.parse_args()

    recs = list(SeqIO.parse(args.fasta, "fasta"))

    ids, mat = compute_distance_matrix(recs)
    
    print("\t" + "\t".join(ids))
    for i, id_i in enumerate(ids):
        row_vals = [id_i]
        for j in range(len(ids)):
            if j <= i:
                row_vals.append("0")
            else:
                row_vals.append(f"{mat[i][j]:.4f}")
        print("\t".join(row_vals))


if __name__ == "__main__":
    main()
