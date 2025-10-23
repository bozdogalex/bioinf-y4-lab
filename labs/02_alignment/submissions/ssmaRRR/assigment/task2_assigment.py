#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Exercițiu: Aliniere globală și locală TP53 (Task 2)
Exemplu rulare:
  python labs/02_alignment/submissions/<handle>/assignment/task2_assigment.py --fasta data/sample/tp53_dna_multi.fasta
"""

from Bio import SeqIO, pairwise2
from Bio.pairwise2 import format_alignment
from pathlib import Path
import argparse

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True, help="Calea către fișierul FASTA")
    ap.add_argument("--i1", type=int, default=0, help="Index prima secvență (implicit 0)")
    ap.add_argument("--i2", type=int, default=1, help="Index a doua secvență (implicit 1)")
    args = ap.parse_args()

    fasta_path = Path(args.fasta)
    if not fasta_path.exists():
        raise SystemExit(f"[eroare] Fișierul nu a fost găsit: {fasta_path}")

    recs = list(SeqIO.parse(str(fasta_path), "fasta"))
    s1, s2 = str(recs[args.i1].seq), str(recs[args.i2].seq)
    id1, id2 = recs[args.i1].id, recs[args.i2].id

    print(f"=== Aliniere pairwise ===")
    print(f"{id1}  vs  {id2}\n")

    # ---------------- GLOBAL ALIGNMENT ----------------
    print(">>> Global alignment (Needleman–Wunsch)\n")
    global_alignments = pairwise2.align.globalms(s1, s2, 2, -1, -2, -0.5)
    top_global = global_alignments[0]
    print(format_alignment(*top_global, full_sequences=False))
    print(f"Global score: {top_global.score}\n")

    # ---------------- LOCAL ALIGNMENT ----------------
    print(">>> Local alignment (Smith–Waterman)\n")
    local_alignments = pairwise2.align.localms(s1, s2, 2, -1, -2, -0.5)
    top_local = local_alignments[0]
    print(format_alignment(*top_local, full_sequences=False))
    print(f"Local score: {top_local.score}\n")

    # ---------------- FRAGMENTS ----------------
    print(">>> Fragmente selectate pentru comparație\n")

    # Fragment global (cu gap-uri)
    global_fragment_seq1 = top_global.seqA[:80]
    global_fragment_seq2 = top_global.seqB[:80]

    # Fragment local (fără gap-uri)
    local_fragment_seq1 = top_local.seqA.replace("-", "")[:80]
    local_fragment_seq2 = top_local.seqB.replace("-", "")[:80]

    print("Global fragment (cu gap-uri):")
    print(f"{id1}: {global_fragment_seq1}")
    print(f"{id2}: {global_fragment_seq2}\n")

    print("Local fragment (fără gap-uri):")
    print(f"{id1}: {local_fragment_seq1}")
    print(f"{id2}: {local_fragment_seq2}\n")

    print("=== Comparație ===")
    if top_local.score > top_global.score:
        print("Local alignment are un scor mai mare pe regiunea conservată.")
    else:
        print("Global alignment acoperă întreaga secvență, dar are mai multe gap-uri.")


if __name__ == "__main__":
    main()