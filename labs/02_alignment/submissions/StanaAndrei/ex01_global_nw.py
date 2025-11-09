#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Exercițiu: Aliniere globală (Needleman–Wunsch)

Scop:
  - Încărcați două secvențe din fișierul descărcat în Lab 1 (data/work/<handle>/lab01/).
  - Implementați pașii de bază ai algoritmului NW pentru a obține alinierea globală.

TODO:
  - Inițializarea matricei de scoruri pentru aliniere globală.
  - Calculul scorurilor celulelor (match, mismatch, gap).

Exemplu rulare:
  python labs/02_alignment/ex02_global_nw.py --fasta data/work/<handle>/lab01/my_tp53.fa --i1 0 --i2 1
"""

from pathlib import Path
import argparse
from Bio import SeqIO


# ===================== TODO: Scoring matrix init =========================================

def init_score_matrix_global(m: int, n: int, gap: int):
    """
    Inițializați matricea (m+1) x (n+1) pentru aliniere globală.
    Pași:
      - creați o listă de liste plină cu 0 (dimensiune (m+1)x(n+1)).
      - prima coloană: [i * gap] pentru i=0..m
      - prima linie: [j * gap] pentru j=0..n
    Returnati matricea.
    """
    # creare matrice cu 0
    score = [[0 for _ in range(n + 1)] for _ in range(m + 1)]

    # prima coloană
    for i in range(1, m + 1):
        score[i][0] = i * gap

    # prima linie
    for j in range(1, n + 1):
        score[0][j] = j * gap

    return score


def score_cell_global(score, i: int, j: int, a: str, b: str, match: int, mismatch: int, gap: int):
    """
    Calculați scorul unei celule (i, j).
    Pași:
      - diagonal = score[i-1][j-1] + (match dacă a == b altfel mismatch)
      - sus      = score[i-1][j] + gap
      - stânga   = score[i][j-1] + gap
    Returnati max(diagonal, sus, stânga).
    """
    # diagonal
    diag_score = score[i - 1][j - 1] + (match if a == b else mismatch)
    # sus (up)
    up_score = score[i - 1][j] + gap
    # stânga (left)
    left_score = score[i][j - 1] + gap

    return max(diag_score, up_score, left_score)


def needleman_wunsch(seq1: str, seq2: str, match=1, mismatch=-1, gap=-2):
    # Implementare simplificată Needleman–Wunsch.
    m, n = len(seq1), len(seq2)

    # Inițializare matrice (voi completați funcția)
    score = init_score_matrix_global(m, n, gap)

    # Umplem matricea celulă cu celulă
    # Observație: buclele merg de la 1 la lungimea secvenței
    for i in range(1, m + 1):
        ai = seq1[i - 1]
        for j in range(1, n + 1):
            bj = seq2[j - 1]
            # apelăm funcția pentru scor
            score[i][j] = score_cell_global(score, i, j, ai, bj, match, mismatch, gap)

    # ================== Backtracking ==================
    # pornim din colțul dreapta-jos (scor[m][n])
    align1, align2 = "", ""
    i, j = m, n
    while i > 0 and j > 0:
        current = score[i][j]
        # recalculăm scorurile vecinilor pentru a decide de unde am venit
        diag = score[i - 1][j - 1] + (match if seq1[i - 1] == seq2[j - 1] else mismatch)
        up   = score[i - 1][j] + gap
        left = score[i][j - 1] + gap

        if current == diag:  # am venit pe diagonală
            align1 = seq1[i - 1] + align1
            align2 = seq2[j - 1] + align2
            i -= 1; j -= 1
        elif current == up:  # am venit de sus
            align1 = seq1[i - 1] + align1
            align2 = "-" + align2
            i -= 1
        else:                # am venit de la stânga
            align1 = "-" + align1
            align2 = seq2[j - 1] + align2
            j -= 1

    # dacă am ajuns la margine, completăm cu gap-uri
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
    args = ap.parse_args()

    fasta_path = Path(args.fasta)
    if not fasta_path.exists():
        raise SystemExit(f"[eroare] Nu găsesc fișierul: {fasta_path}")

    s1, s2, id1, id2 = load_two_sequences(fasta_path, args.i1, args.i2)
    a1, a2, sc = needleman_wunsch(s1, s2)

    print("=== Aliniere globală (NW) ===")
    print(f"{id1}  vs  {id2}")
    print(a1)
    print(a2)
    print("Score:", sc)


if __name__ == "__main__":
    main()
