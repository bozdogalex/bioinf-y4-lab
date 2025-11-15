#!/usr/bin/env python3
"""
Task 2: Pairwise Alignments (Global vs Local)
Compară aliniamentele globale și locale pentru două secvențe.
"""

from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio import SeqIO
import os


def compare_alignments(seq1, seq2, id1, id2):
    """
    Compară aliniamentele globale și locale pentru două secvențe.
    """
    print(f"\n{'='*80}")
    print(f"ANALIZĂ PAIRWISE: {id1} vs {id2}")
    print(f"{'='*80}")
    print(f"Lungime seq1: {len(seq1)}")
    print(f"Lungime seq2: {len(seq2)}")
    
    # ========== ALINIERE GLOBALĂ ==========
    print(f"\n{'─'*80}")
    print("1. ALINIERE GLOBALĂ (Needleman-Wunsch)")
    print(f"{'─'*80}")
    
    # Aliniere globală simplă (match=1, mismatch=0, gap=0)
    global_alignments = pairwise2.align.globalxx(seq1, seq2, one_alignment_only=True)
    
    if global_alignments:
        global_align = global_alignments[0]
        print(format_alignment(*global_align, full_sequences=True))
        
        # Calculează statistici
        aligned1, aligned2, score, begin, end = global_align
        matches = sum(c1 == c2 for c1, c2 in zip(aligned1, aligned2) if c1 != '-' and c2 != '-')
        gaps1 = aligned1.count('-')
        gaps2 = aligned2.count('-')
        
        print(f"\n📊 Statistici aliniere globală:")
        print(f"  - Scor: {score}")
        print(f"  - Lungime aliniere: {len(aligned1)}")
        print(f"  - Match-uri: {matches}")
        print(f"  - Gap-uri în seq1: {gaps1}")
        print(f"  - Gap-uri în seq2: {gaps2}")
        print(f"  - Identitate: {100 * matches / len(aligned1):.2f}%")
    
    # ========== ALINIERE LOCALĂ ==========
    print(f"\n{'─'*80}")
    print("2. ALINIERE LOCALĂ (Smith-Waterman)")
    print(f"{'─'*80}")
    
    # Aliniere locală simplă
    local_alignments = pairwise2.align.localxx(seq1, seq2, one_alignment_only=True)
    
    if local_alignments:
        local_align = local_alignments[0]
        print(format_alignment(*local_align, full_sequences=False))
        
        # Calculează statistici
        aligned1, aligned2, score, begin, end = local_align
        matches = sum(c1 == c2 for c1, c2 in zip(aligned1, aligned2) if c1 != '-' and c2 != '-')
        gaps1 = aligned1.count('-')
        gaps2 = aligned2.count('-')
        
        print(f"\n📊 Statistici aliniere locală:")
        print(f"  - Scor: {score}")
        print(f"  - Lungime aliniere: {len(aligned1)}")
        print(f"  - Poziție start: {begin}")
        print(f"  - Poziție end: {end}")
        print(f"  - Match-uri: {matches}")
        print(f"  - Gap-uri în seq1: {gaps1}")
        print(f"  - Gap-uri în seq2: {gaps2}")
        print(f"  - Identitate: {100 * matches / len(aligned1):.2f}%")
    
    # ========== COMPARAȚIE ==========
    print(f"\n{'─'*80}")
    print("3. COMPARAȚIE GLOBAL vs LOCAL")
    print(f"{'─'*80}")
    
    if global_alignments and local_alignments:
        global_align = global_alignments[0]
        local_align = local_alignments[0]
        
        print(f"Diferențe cheie:")
        print(f"  • Global forțează alinierea întregii secvențe")
        print(f"  • Local găsește doar cea mai bună regiune de similaritate")
        print(f"  • Scor global: {global_align[2]:.1f} vs Scor local: {local_align[2]:.1f}")
        print(f"  • Lungime global: {len(global_align[0])} vs local: {len(local_align[0])}")
        
        # Găsește un fragment unde există diferențe majore
        print(f"\n🔍 Exemplu de diferență:")
        print(f"  Alinierea globală poate introduce gap-uri la capete pentru a")
        print(f"  alinia întreaga secvență, în timp ce alinierea locală se")
        print(f"  concentrează doar pe regiunea cu cea mai mare similaritate.")


def advanced_alignment(seq1, seq2, id1, id2):
    """
    Aliniere cu parametri mai sofisticați (scoring matrix).
    """
    print(f"\n{'='*80}")
    print(f"ALINIERE AVANSATĂ CU SCORING")
    print(f"{'='*80}")
    
    # Parametri: match=2, mismatch=-1, gap_open=-2, gap_extend=-0.5
    print("\nParametri: match=+2, mismatch=-1, gap_open=-2, gap_extend=-0.5")
    
    # Global cu scoring
    global_scored = pairwise2.align.globalms(
        seq1, seq2, 
        match=2, mismatch=-1, 
        open=-2, extend=-0.5,
        one_alignment_only=True
    )
    
    if global_scored:
        print(f"\n🔹 Global scored alignment:")
        print(format_alignment(*global_scored[0], full_sequences=False))
        print(f"Scor: {global_scored[0][2]:.2f}")
    
    # Local cu scoring
    local_scored = pairwise2.align.localms(
        seq1, seq2,
        match=2, mismatch=-1,
        open=-2, extend=-0.5,
        one_alignment_only=True
    )
    
    if local_scored:
        print(f"\n🔹 Local scored alignment:")
        print(format_alignment(*local_scored[0], full_sequences=False))
        print(f"Scor: {local_scored[0][2]:.2f}")


def main():
    # Folosim secvențele proteice TP53
    fasta_file = "tp53_protein_multi.fasta"
    
    if not os.path.exists(fasta_file):
        print(f"EROARE: Fișierul {fasta_file} nu există!")
        print("Asigură-te că fișierul este în același director cu scriptul.")
        return
    
    # Citește secvențele
    records = list(SeqIO.parse(fasta_file, "fasta"))
    
    if len(records) < 2:
        print("EROARE: Sunt necesare cel puțin 2 secvențe!")
        return
    
    # Selectează primele două secvențe (Human și Mouse TP53)
    seq1 = str(records[0].seq)
    seq2 = str(records[1].seq)
    id1 = records[0].id
    id2 = records[1].id
    
    print(f"Fișier procesat: {fasta_file}")
    print(f"Secvența 1: {id1}")
    print(f"Secvența 2: {id2}")
    
    # Compară aliniamentele
    compare_alignments(seq1, seq2, id1, id2)
    
    # Aliniere avansată
    advanced_alignment(seq1, seq2, id1, id2)
    
    print(f"\n{'='*80}")
    print("💡 CONCLUZIE")
    print(f"{'='*80}")
    print("Alinierea GLOBALĂ este utilă când:")
    print("  • Secvențele au lungimi similare")
    print("  • Vrem să comparăm întreaga structură")
    print("  • Căutăm relații evolutive între gene întregi")
    print()
    print("Alinierea LOCALĂ este utilă când:")
    print("  • Căutăm domenii conservate")
    print("  • Secvențele au lungimi foarte diferite")
    print("  • Ne interesează doar regiunile similare (ex: motive funcționale)")


if __name__ == "__main__":
    main()
