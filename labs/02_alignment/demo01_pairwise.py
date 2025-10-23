#!/usr/bin/env python
"""
Demo: aliniere globală și locală cu Biopython (PairwiseAligner).
- Înlocuiește pairwise2 (depreciat) cu noul modul Bio.Align.
- Folosește aceleași scoruri: +1 match, -1 mismatch, -1 gap.
Rulare:
  python labs/02_alignment/demo01_pairwise_biopython_new.py --fasta data/sample/tp53_dna_multi.fasta
"""

import argparse
from Bio import SeqIO
from Bio.Align import PairwiseAligner

def take_two_short_subseqs(fasta_path, k=7):
    """Extrage primele două secvențe dintr-un FASTA și le taie la lungimea k."""
    recs = [r for r in SeqIO.parse(fasta_path, "fasta")]
    if len(recs) < 2:
        raise ValueError("Need at least 2 sequences in the FASTA.")
    a = str(recs[0].seq)[:k]
    b = str(recs[1].seq)[:k]
    return a, b

def run_alignment(A, B, mode="global"):
    """Rulează alinierea între A și B folosind PairwiseAligner."""
    aligner = PairwiseAligner()
    aligner.mode = mode  # 'global' sau 'local'
    aligner.match_score = 1
    aligner.mismatch_score = -1
    aligner.open_gap_score = -1
    aligner.extend_gap_score = -1
    alignments = aligner.align(A, B)
    return alignments

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True, help="path to FASTA file")
    ap.add_argument("--k", type=int, default=7, help="substring length")
    args = ap.parse_args()

    A, B = take_two_short_subseqs(args.fasta, k=args.k)

    print("[INPUT]")
    print("A:", A)
    print("B:", B)

    print("\n[GLOBAL] top alignment:")
    global_alignments = run_alignment(A, B, mode="global")
    print(global_alignments[0])

    print("[LOCAL] top alignment:")
    local_alignments = run_alignment(A, B, mode="local")
    print(local_alignments[0])

if __name__ == "__main__":
    main()
