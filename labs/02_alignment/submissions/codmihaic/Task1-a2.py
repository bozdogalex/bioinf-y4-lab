#!/usr/bin/env python
import argparse
from itertools import combinations
from Bio import SeqIO

def hamming_equal(a, b):
    return sum(x != y for x, y in zip(a, b))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--ids", type=str, default="", help="CSV de IDs din FASTA: NG_017013.2,NC_000017.11,NC_060941.1")
    ap.add_argument("--maxlen", type=int, default=2000, help="Taie secventele la maxlen pentru calcul, daca sunt prea mari")
    args = ap.parse_args()

    wanted = {x.strip() for x in args.ids.split(",") if x.strip()} if args.ids else None

    recs = [r for r in SeqIO.parse(args.fasta, "fasta") if (not wanted or r.id in wanted)]
    ids = [r.id for r in recs]
    seqs = [str(r.seq)[:args.maxlen] for r in recs]  # limitare pentru viteza si memorie cat de cat eficienta

    print("pair,hamming,p_distance,len_used")
    for (i, j) in combinations(range(len(seqs)), 2):
        a, b = seqs[i], seqs[j]
        L = min(len(a), len(b))
        a2, b2 = a[:L], b[:L]
        d_h = hamming_equal(a2, b2)
        d_p = d_h / float(L) if L > 0 else 0.0
        print(f"{ids[i]}-{ids[j]},{d_h},{d_p:.6f},{L}")

if __name__ == "__main__":
    main()
