#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Exercițiu: Aliniere globală, locală și semiglobală (TP53)
Exemplu rulare:
python labs/02_alignment/submissions/<handle>/assignment/task_bonus_assigment.py --fasta data/sample/tp53_dna_multi.fasta  
"""

from Bio import SeqIO, pairwise2
from Bio.pairwise2 import format_alignment
from pathlib import Path
import argparse

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True, help="Calea către fișierul FASTA (ex: data/sample/tp53_dna_multi.fasta)")
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

    frag_g1 = top_global.seqA[:100]
    frag_g2 = top_global.seqB[:100]
    print("🔹 Fragment Global:")
    print(f"{id1}: {frag_g1}")
    print(f"{id2}: {frag_g2}\n")

    # ---------------- LOCAL ALIGNMENT ----------------
    print(">>> Local alignment (Smith–Waterman)\n")
    local_alignments = pairwise2.align.localms(s1, s2, 2, -1, -2, -0.5)
    top_local = local_alignments[0]
    print(format_alignment(*top_local, full_sequences=False))
    print(f"Local score: {top_local.score}\n")

    frag_l1 = top_local.seqA[:100]
    frag_l2 = top_local.seqB[:100]
    print("🔹 Fragment Local:")
    print(f"{id1}: {frag_l1}")
    print(f"{id2}: {frag_l2}\n")

    # ---------------- SEMIGLOBAL ALIGNMENT ----------------
    print(">>> Semiglobal alignment (no end-gap penalties)\n")
    semi_alignments = pairwise2.align.globalms(
        s1, s2,
        2,   # match
        -1,  # mismatch
        -2,  # gap open
        -0.5,# gap extend
        penalize_end_gaps=(False, False)  # cheie: fără penalizări la capete
    )
    top_semi = semi_alignments[0]
    print(format_alignment(*top_semi, full_sequences=False))
    print(f"Semiglobal score: {top_semi.score}\n")

    frag_s1 = top_semi.seqA[:100]
    frag_s2 = top_semi.seqB[:100]
    print("Fragment Semiglobal:")
    print(f"{id1}: {frag_s1}")
    print(f"{id2}: {frag_s2}\n")

    # ---------------- COMPARAȚIE ----------------
    print("=== Comparație ===")
    if top_local.score > top_global.score:
        print("Local alignment are un scor mai mare pe regiunea conservată.")
    else:
        print("Global alignment acoperă întreaga secvență, dar are mai multe gap-uri.")

    if top_semi.score > top_global.score:
        print("Semiglobal alignment evită penalizarea capetelor și e util pentru secvențe parțiale.")
    else:
        print("Global alignment rămâne mai potrivit pentru secvențe complete omoloage.")


if __name__ == "__main__":
    main()