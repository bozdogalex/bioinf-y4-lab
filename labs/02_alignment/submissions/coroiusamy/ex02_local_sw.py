#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Exercițiu: Aliniere locală (Smith–Waterman)

Scop:
  - Încărcați două secvențe din fișierul descărcat în Lab 1 (data/work/<handle>/lab01/).
  - Implementați pașii de bază ai algoritmului SW pentru a obține alinierea locală.

TODO:
  - Inițializarea matricei locale (Rezolvat)
  - Scoring celulelor (match, mismatch, gap, cu max(0,...)). (Rezolvat)

Exemplu Rulare:
  python labs/02_alignment/ex03_local_sw.py --fasta data/work/<handle>/lab01/my_tp53.fa --i1 0 --i2 1
"""

from pathlib import Path
import argparse
from Bio import SeqIO


# ===================== Matrix initiation =========================================

def init_score_matrix_local(m: int, n: int):
    """
    TODO 1: Inițializați matricea (m+1) x (n+1) cu toate valorile = 0.
    """
    # REZOLVAT: Creează o listă de (m+1) liste, fiecare conținând (n+1) zerouri.
    return [[0 for _ in range(n + 1)] for _ in range(m + 1)]


def score_cell_local(score, i: int, j: int, a: str, b: str, match: int, mismatch: int, gap: int):
    """
    TODO 2: Calculați scorul unei celule (i, j).
    Pași:
      - diagonal = score[i-1][j-1] + (match dacă a == b altfel mismatch)
      - sus      = score[i-1][j] + gap
      - stânga   = score[i][j-1] + gap
    Rezultat = max(0, diagonal, sus, stânga).
    """
    # REZOLVAT: Calculăm scorul de potrivire/nepotrivire
    match_mismatch_score = match if a == b else mismatch
    
    # Calculăm scorurile posibile
    diagonal = score[i - 1][j - 1] + match_mismatch_score
    sus = score[i - 1][j] + gap
    stanga = score[i][j - 1] + gap
    
    # Algoritmul Smith-Waterman resetează la 0 dacă toate scorurile sunt negative
    return max(0, diagonal, sus, stanga)



def smith_waterman(seq1: str, seq2: str, match=3, mismatch=-3, gap=-2):
    # Implementare simplificată Smith–Waterman.
    m, n = len(seq1), len(seq2)

    # Inițializare matrice
    score = init_score_matrix_local(m, n)

    max_score = 0
    max_pos = (0, 0)

    # Umplem matricea celulă cu celulă
    # + memorăm cea mai mare valoare și poziția asociată
    for i in range(1, m + 1):
        ai = seq1[i - 1]
        for j in range(1, n + 1):
            bj = seq2[j - 1]
            
            # TODO 2 a fost înlocuit cu apelul funcției implementate
            score[i][j] = score_cell_local(score, i, j, ai, bj, match, mismatch, gap)

            if score[i][j] > max_score:
                max_score = score[i][j]
                max_pos = (i, j)

    # ================== Backtracking ==================
    # pornim din celula cu scor maxim și mergem înapoi
    align1, align2 = "", ""
    i, j = max_pos
    while i > 0 and j > 0 and score[i][j] > 0:
        # recalculăm scorurile vecinilor pentru a decide direcția
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
        else:  # stânga
            align1 = "-" + align1
            align2 = seq2[j - 1] + align2
            j -= 1

    return align1, align2, max_score


def load_two_sequences(fasta_path: Path, i1: int, i2: int):
    """
    Încărcăm FASTA și alegem două secvențe după index.
    """
    recs = list(SeqIO.parse(str(fasta_path), "fasta"))
    if len(recs) < 2:
        raise SystemExit("[eroare] Fișierul trebuie să conțină cel puțin 2 secvențe.")
    if not (0 <= i1 < len(recs) and 0 <= i2 < len(recs)):
        raise SystemExit(f"[eroare] Indici invalizi (0..{len(recs)-1}).")
    return str(recs[i1].seq), str(recs[i2].seq), recs[i1].id, recs[i2].id


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True, help="Cale către FASTA-ul propriu din data/work/<handle>/lab01/")
    ap.add_argument("--i1", type=int, default=0, help="Index prima secvență (implicit 0)")
    ap.add_argument("--i2", type=int, default=1, help="Index a doua secvență (implicit 1)")
    
    # parametrul --limit
    ap.add_argument("--limit", type=int, default=0, 
                    help="Limitează secvențele la primele N baze (0 = fără limită, implicit 0)")

    args = ap.parse_args()

    fasta_path = Path(args.fasta)
    if not fasta_path.exists():
        raise SystemExit(f"[eroare] Nu găsesc fișierul: {fasta_path}")

    s1, s2, id1, id2 = load_two_sequences(fasta_path, args.i1, args.i2)
    
    # Logica pentru a trunchia secvențele dacă --limit este folosit
    if args.limit > 0:
        s1 = s1[:args.limit]
        s2 = s2[:args.limit]
        print(f"[Info] Secvențele au fost limitate la primele {args.limit} baze.")
    
    print("======================================================================")
    print("Aliniere locală Smith-Waterman")
    print("======================================================================")
    print(f"Secvență 1: {id1} (lungime: {len(s1)} bp)")
    print(f"Secvență 2: {id2} (lungime: {len(s2)} bp)")
    print("Se calculează alinierea...")

    a1, a2, sc = smith_waterman(s1, s2)

    print("\n=== Rezultat Aliniere ===")
    print(f"{id1} vs {id2}")
    print(a1)
    print(a2)
    print("Score:", sc)


if __name__ == "__main__":
    main()