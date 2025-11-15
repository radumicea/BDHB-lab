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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--i1", type=int, default=0)
    ap.add_argument("--i2", type=int, default=1)
    args = ap.parse_args()

    s1, s2, id1, id2 = load_two_sequences(Path(args.fasta), args.i1, args.i2)

    global_aln = run_alignment("global", s1, s2)
    local_aln  = run_alignment("local",  s1, s2)

    print("Global alignment")
    print(id1, "vs", id2)
    # print(global_aln)
    print("Score:", global_aln.score)

    print("\nLocal alignment")
    print(id1, "vs", id2)
    # print(local_aln)
    print("Score:", local_aln.score)


if __name__ == "__main__":
    main()
