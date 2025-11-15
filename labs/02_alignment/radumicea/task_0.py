# Slice the TP53 regions to be compared

from pathlib import Path
from Bio import SeqIO
import argparse

def extract_tp53(in_fasta: Path, out_fasta: Path):
    recs = {rec.id: rec for rec in SeqIO.parse(in_fasta, "fasta")}

    out_records = []

    # NG_017013.2: plus strand, 5001..24149
    ng = recs["NG_017013.2"]
    tp53_ng_seq = ng.seq[5001-1:24149]  # 1-based -> 0-based
    ng_rec = ng[:]
    ng_rec.seq = tp53_ng_seq
    ng_rec.id = "NG_017013.2_TP53_5001_24149"
    ng_rec.description = "TP53 region NG_017013.2:5001-24149 (forward)"
    out_records.append(ng_rec)

    # NC_000017.11: minus strand, 7668421..7687490
    nc17 = recs["NC_000017.11"]
    tp53_nc17_seq = nc17.seq[7668421-1:7687490].reverse_complement()
    nc17_rec = nc17[:]
    nc17_rec.seq = tp53_nc17_seq
    nc17_rec.id = "NC_000017.11_TP53_7668421_7687490_rev"
    nc17_rec.description = "TP53 region NC_000017.11:7668421-7687490 (reverse complement)"
    out_records.append(nc17_rec)

    # NC_060941.1: minus strand, 7572544..7591594
    nc60 = recs["NC_060941.1"]
    tp53_nc60_seq = nc60.seq[7572544-1:7591594].reverse_complement()
    nc60_rec = nc60[:]
    nc60_rec.seq = tp53_nc60_seq
    nc60_rec.id = "NC_060941.1_TP53_7572544_7591594_rev"
    nc60_rec.description = "TP53 region NC_060941.1:7572544-7591594 (reverse complement)"
    out_records.append(nc60_rec)

    SeqIO.write(out_records, out_fasta, "fasta")
    print(f"[ok] Wrote {len(out_records)} TP53-region sequences to: {out_fasta}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in_fasta", required=True)
    ap.add_argument("--out_fasta", required=True)
    args = ap.parse_args()

    extract_tp53(args.in_fasta, args.out_fasta)


if __name__ == "__main__":
    main()