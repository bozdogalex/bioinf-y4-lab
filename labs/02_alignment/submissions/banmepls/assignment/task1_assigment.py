#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Exercițiu: Distanțe perechi (Hamming / p-distance)

Scop:
  - Calculați distanța Hamming (pentru secvențe de aceeași lungime)
    sau p-distance (proporția pozițiilor diferite).
  - Produceți matricea de distanțe pentru ≥3 secvențe din fișierul FASTA.

Exemplu rulare:
  python labs/02_alignment/submissions/<handle>/assignment/task1_assigment.py --fasta data/sample/tp53_dna_multi.fasta
"""

from pathlib import Path
import argparse
from Bio import SeqIO
import pandas as pd
import math


# ====================== Funcții utile =====================================

def hamming_distance(seq1, seq2):
    """Calculează distanța Hamming între două secvențe de aceeași lungime."""
    if len(seq1) != len(seq2):
        raise ValueError("Hamming distance requires equal-length sequences.")
    return sum(a != b for a, b in zip(seq1, seq2))


def p_distance(seq1, seq2):
    """Calculează proporția pozițiilor diferite (trunchiată la lungimea minimă)."""
    m = min(len(seq1), len(seq2))
    if m == 0:
        return math.nan
    diffs = sum(a != b for a, b in zip(seq1[:m], seq2[:m]))
    return diffs / m


# ====================== Calcul matrice =====================================

def compute_distance_matrix(records):
    ids = list(records.keys())
    n = len(ids)
    matrix = [[None] * n for _ in range(n)]
    pair_info = []

    for i in range(n):
        for j in range(i + 1, n):
            s1, s2 = records[ids[i]], records[ids[j]]
            try:
                if len(s1) == len(s2):
                    dist = hamming_distance(s1, s2)
                    method = "Hamming"
                else:
                    dist = p_distance(s1, s2)
                    method = "p-distance (truncated)"
            except Exception as e:
                dist, method = math.nan, f"Error: {e}"

            matrix[i][j] = dist
            pair_info.append((ids[i], ids[j], dist, method))

    df = pd.DataFrame(matrix, index=ids, columns=ids)
    return df, pair_info


# ====================== Funcție principală =====================================

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True, help="Cale către FASTA-ul propriu din data/work/<handle>/lab01/")
    args = ap.parse_args()

    fasta_path = Path(args.fasta)
    if not fasta_path.exists():
        raise SystemExit(f"[eroare] Nu găsesc fișierul: {fasta_path}")

    records = {rec.id: str(rec.seq) for rec in SeqIO.parse(str(fasta_path), "fasta")}
    if len(records) < 3:
        raise SystemExit("[eroare] Fișierul trebuie să conțină cel puțin 3 secvențe.")

    df, pairs = compute_distance_matrix(records)
    df.to_csv("distance_matrix.csv", float_format="%.5f")

    # găsim perechea cu distanța minimă
    valid = [(a, b, d, m) for a, b, d, m in pairs if d == d]  # eliminăm NaN
    if not valid:
        print("[info] Nicio pereche validă de comparat.")
        return

    closest = min(valid, key=lambda x: x[2])
    print("=== Matrice de distanțe (triunghi superior) ===")
    print(df.fillna("-"))
    print("\nCea mai apropiată pereche:")
    print(f"{closest[0]} vs {closest[1]} → distanță = {closest[2]:.5f} ({closest[3]})")


if __name__ == "__main__":
    main()
