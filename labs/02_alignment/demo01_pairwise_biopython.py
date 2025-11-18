#!/usr/bin/env python
"""
Demo: aliniere globală și locală cu Biopython (pairwise).
- Refolosim datele din data/sample/; extragem subsecvențe scurte pentru debugging.
Rulare:
  python labs/02_alignment/demo01_pairwise_biopython.py --fasta data/sample/tp53_dna_multi.fasta
"""
import argparse
from Bio import SeqIO, pairwise2

def take_two_short_subseqs(fasta_path, k=7):
    recs = [r for r in SeqIO.parse(fasta_path, "fasta")]
    if len(recs) < 2:
        raise ValueError("Need at least 2 sequences in the FASTA.")
    a = str(recs[0].seq)[:k]
    b = str(recs[1].seq)[:k]
    return a, b

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--k", type=int, default=7, help="substring length")
    args = ap.parse_args()

    A, B = take_two_short_subseqs(args.fasta, k=args.k)

    # Scor simplu: +1 match, -1 mismatch, -1 gap
    global_alignments = pairwise2.align.globalms(A, B, 1, -1, -1, -1)
    local_alignments = pairwise2.align.localms(A, B, 1, -1, -1, -1)

    print("[INPUT]")
    print("A:", A)
    print("B:", B)

    print("\n[GLOBAL] top alignment:")
    print(pairwise2.format_alignment(*global_alignments[0]))

    print("[LOCAL] top alignment:")
    print(pairwise2.format_alignment(*local_alignments[0]))

if __name__ == "__main__":
    main()
