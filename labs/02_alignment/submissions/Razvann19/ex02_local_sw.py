gghhhh#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Exercițiu: Aliniere locală (Smith–Waterman)

Scop:
  - Încărcați două secvențe din fișierul descărcat în Lab 1 (data/work/<handle>/lab01/).
  - Implementați pașii de bază ai algoritmului SW pentru a obține alinierea locală.

TODO:
  - Inițializarea matricei locale
  - Scoring celulelor (match, mismatch, gap, cu max(0,...)).

Exemplu Rulare:
  python labs/02_alignment/ex03_local_sw.py --fasta data/work/<handle>/lab01/my_tp53.fa --i1 0 --i2 1
"""

from pathlib import Path
import argparse
from Bio import SeqIO


# ===================== Matrix initiation =========================================

def init_score_matrix_local(m: int, n: int):
    # matrice (m+1) x (n+1) inițializată cu 0
    return [[0] * (n + 1) for _ in range(m + 1)]



def score_cell_local(score, i: int, j: int, a: str, b: str, match: int, mismatch: int, gap: int):
    diag = score[i - 1][j - 1] + (match if a == b else mismatch)
    up   = score[i - 1][j] + gap
    left = score[i][j - 1] + gap
    return max(0, diag, up, left)



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
    args = ap.parse_args()

    fasta_path = Path(args.fasta)
    if not fasta_path.exists():
        raise SystemExit(f"[eroare] Nu găsesc fișierul: {fasta_path}")

    s1, s2, id1, id2 = load_two_sequences(fasta_path, args.i1, args.i2)
    a1, a2, sc = smith_waterman(s1, s2)

    print("=== Aliniere locală (SW) ===")
    print(f"{id1}  vs  {id2}")
    print(a1)
    print(a2)
    print("Score:", sc)


if __name__ == "__main__":
    main()

"""

Am facut o copie la my_tp53.fa mai sccurta doar cu primele 500 de elemnte pentru ca nu am reusit sa rulez pe tot fisierul
Imi da "Terminated" dupa cateva secunde

=== Aliniere locală (SW) ===
NG_017013.2  vs  NC_060941.1
CTCGAACTCCTGACCTCAGGTGATCC--A-CCT-GCCTCAGCCT--CCCAAAGTGCTGGGATTA--C-AGGA-GTCAGCC--ACCGC--ACCC-AGCCCCAA--CTAATTTTTGTATTTGA--GTAGAGACAGGGTTTTACCATGTTGGCCAGGCT-GGTCTAAAACTCTTC-ACCTCAGGTGATCC-ACCCATCTCAGCCT--CCC-AAAGTGTTGGGATTACAGGCGTGAGCC--ACCGT--GCCTGGCCCTGGATTTCACTCTTGCCC--ACCCATAAACCATTCACTCTTCTGTTTTAA----AACTCTCTTGGCCGGGCGCAGTGGCTCATGCCTGTAATCCC-AGCACT
CT--AAC-CCTAACC-CA--TAACCCTAACCCTAACCT-ACCCTAACCCTAA-CCCT---A--ACCCTA--ACCTAACCCTAACC-CTAACCCTAACCCTAACCCTAA---CCCTA----ACCCT--A-AC----CCTAACC---CTAACC---CTAACCCT--AAC-C--CAACC-C---T-AACCTAACC--CTAACCCTAACCCTAAA-CCCT---A--AC---CCTAACCCTAACCCTAACCCTAACCCT--A----ACCCTAACCCTAACCC-T-AACC-CTAAC-C--C-----TAACCCTAAC-C-C-TAACC---C-TA-ACCCT-A-ACC-CTAA-CCCTAACCCT
Score: 125

"""