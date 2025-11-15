#!/usr/bin/env python3
"""
Task 1: Hamming Distance Matrix
Calculează distanțele Hamming pentru toate perechile de secvențe de aceeași lungime.
"""

from Bio import SeqIO
import numpy as np
import os


def hamming_distance(seq1, seq2):
    """
    Calculează distanța Hamming între două secvențe de aceeași lungime.
    
    Args:
        seq1: Prima secvență (string sau Seq object)
        seq2: A doua secvență (string sau Seq object)
    
    Returns:
        Distanța Hamming (numărul de poziții diferite)
    """
    if len(seq1) != len(seq2):
        raise ValueError(f"Secvențele trebuie să aibă aceeași lungime: {len(seq1)} vs {len(seq2)}")
    
    return sum(c1 != c2 for c1, c2 in zip(seq1, seq2))


def calculate_distance_matrix(sequences):
    """
    Calculează matricea de distanțe Hamming pentru un set de secvențe.
    
    Args:
        sequences: Lista de tuple (id, seq)
    
    Returns:
        Matricea de distanțe ca numpy array
    """
    n = len(sequences)
    matrix = np.zeros((n, n))
    
    for i in range(n):
        for j in range(i + 1, n):
            dist = hamming_distance(sequences[i][1], sequences[j][1])
            matrix[i][j] = dist
            matrix[j][i] = dist
    
    return matrix


def print_distance_matrix(sequences, matrix):
    """
    Afișează matricea de distanțe într-un format lizibil.
    """
    n = len(sequences)
    ids = [seq[0] for seq in sequences]
    
    # Header
    print("\nMatricea de distanțe Hamming:")
    print("=" * 80)
    print(f"{'':15}", end="")
    for seq_id in ids:
        print(f"{seq_id:>15}", end="")
    print()
    
    # Rows
    for i in range(n):
        print(f"{ids[i]:15}", end="")
        for j in range(n):
            if i == j:
                print(f"{'0':>15}", end="")
            elif j > i:
                print(f"{int(matrix[i][j]):>15}", end="")
            else:
                print(f"{'':>15}", end="")
        print()
    print("=" * 80)


def find_closest_pair(sequences, matrix):
    """
    Găsește perechea de secvențe cu cea mai mică distanță.
    """
    n = len(sequences)
    min_dist = float('inf')
    min_pair = None
    
    for i in range(n):
        for j in range(i + 1, n):
            if matrix[i][j] < min_dist:
                min_dist = matrix[i][j]
                min_pair = (i, j)
    
    return min_pair, min_dist


def main():
    # Folosim secvențele proteice TP53 din fișierul multi-FASTA
    # Acestea au lungimi similare și sunt ideale pentru Hamming distance
    fasta_file = "tp53_protein_multi.fasta"
    
    if not os.path.exists(fasta_file):
        print(f"EROARE: Fișierul {fasta_file} nu există!")
        print("Asigură-te că fișierul este în același director cu scriptul.")
        return
    
    # Citește secvențele
    records = list(SeqIO.parse(fasta_file, "fasta"))
    
    # Selectează primele 3 secvențe pentru analiză
    # Vom folosi doar primele 3 secvențe proteice (Human, Mouse, XPO1 partial)
    selected_records = records[:3]
    
    print(f"Fișier procesat: {fasta_file}")
    print(f"Număr total de secvențe: {len(records)}")
    print(f"Secvențe selectate pentru analiză: {len(selected_records)}")
    print()
    
    # Afișează informații despre secvențe
    for i, rec in enumerate(selected_records, 1):
        print(f"{i}. {rec.id}: {len(rec.seq)} caractere")
    
    # Verifică dacă toate secvențele au aceeași lungime
    lengths = [len(rec.seq) for rec in selected_records]
    if len(set(lengths)) > 1:
        print("\n⚠️  ATENȚIE: Secvențele au lungimi diferite!")
        print("Pentru Hamming distance, vom trunchia la lungimea minimă.")
        min_len = min(lengths)
        print(f"Lungime minimă: {min_len}")
        sequences = [(rec.id, str(rec.seq[:min_len])) for rec in selected_records]
    else:
        sequences = [(rec.id, str(rec.seq)) for rec in selected_records]
    
    # Calculează matricea de distanțe
    matrix = calculate_distance_matrix(sequences)
    
    # Afișează matricea
    print_distance_matrix(sequences, matrix)
    
    # Găsește perechea cea mai apropiată
    closest_pair, min_dist = find_closest_pair(sequences, matrix)
    
    print(f"\n📊 REZULTATE:")
    print(f"Perechea cea mai apropiată:")
    print(f"  - {sequences[closest_pair[0]][0]}")
    print(f"  - {sequences[closest_pair[1]][0]}")
    print(f"  - Distanța Hamming: {int(min_dist)}")
    print(f"  - Similaritate: {100 * (1 - min_dist / len(sequences[0][1])):.2f}%")
    
    print("\n💡 INTERPRETARE:")
    print("Distanța Hamming reprezintă numărul de poziții la care cele două secvențe diferă.")
    print("Cu cât distanța este mai mică, cu atât secvențele sunt mai asemănătoare.")
    print("Perechea cu distanța minimă sugerează o relație evolutivă mai apropiată")
    print("sau o funcție biologică mai conservată.")


if __name__ == "__main__":
    main()
