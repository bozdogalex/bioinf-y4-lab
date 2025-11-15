#!/usr/bin/env python3
"""
BONUS: Semiglobal Alignment
Implementare și demonstrație a alinierii semiglobale.
"""

from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio import SeqIO
import os


def semiglobal_alignment_manual(seq1, seq2, match=1, mismatch=-1, gap=-1):
    """
    Implementare simplificată a alinierii semiglobale.
    
    Semiglobal = nu penalizăm gap-urile la capete.
    Util când:
    - O secvență este subsecvență a alteia
    - Căutăm un fragment într-o secvență mai mare
    - Comparăm gene cu lungimi foarte diferite
    """
    m, n = len(seq1), len(seq2)
    
    # Inițializare matrice de scoring
    # Prima linie și coloană nu sunt penalizate (semiglobal!)
    score = [[0 for _ in range(n + 1)] for _ in range(m + 1)]
    
    # Nu penalizăm gap-urile la început (semiglobal)
    # Spre deosebire de global care ar avea: score[i][0] = i * gap
    
    # Fill matrix
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            if seq1[i-1] == seq2[j-1]:
                diag = score[i-1][j-1] + match
            else:
                diag = score[i-1][j-1] + mismatch
            
            up = score[i-1][j] + gap
            left = score[i][j-1] + gap
            
            score[i][j] = max(diag, up, left, 0)  # 0 pentru a permite restart
    
    # Găsește maximul în ultima linie/coloană (semiglobal)
    max_score = 0
    max_pos = (m, n)
    
    # Check ultima linie
    for j in range(n + 1):
        if score[m][j] > max_score:
            max_score = score[m][j]
            max_pos = (m, j)
    
    # Check ultima coloană
    for i in range(m + 1):
        if score[i][n] > max_score:
            max_score = score[i][n]
            max_pos = (i, n)
    
    return max_score, max_pos


def demonstrate_semiglobal():
    """
    Demonstrație a diferenței dintre global, local și semiglobal.
    """
    # Exemple ilustrative
    seq1 = "ACGTACGTACGT"  # secvență mai scurtă
    seq2 = "TTTTACGTACGTACGTAAAA"  # secvență mai lungă cu seq1 în mijloc
    
    print(f"{'='*80}")
    print("DEMONSTRAȚIE: GLOBAL vs LOCAL vs SEMIGLOBAL")
    print(f"{'='*80}\n")
    
    print(f"Secvența 1: {seq1} (lungime: {len(seq1)})")
    print(f"Secvența 2: {seq2} (lungime: {len(seq2)})")
    print()
    print("Observație: seq1 este conținută în seq2, dar cu flanking regions")
    print()
    
    # ========== GLOBAL ==========
    print(f"{'─'*80}")
    print("1. ALINIERE GLOBALĂ")
    print(f"{'─'*80}")
    print("Forțează alinierea întregii secvențe, penalizând gap-urile")
    print()
    
    global_align = pairwise2.align.globalms(
        seq1, seq2,
        match=1, mismatch=-1,
        open=-2, extend=-0.5,
        one_alignment_only=True
    )
    
    if global_align:
        print(format_alignment(*global_align[0]))
        print(f"Scor: {global_align[0][2]:.2f}")
        print("⚠️  Penalizează gap-urile la capete!")
    
    # ========== LOCAL ==========
    print(f"\n{'─'*80}")
    print("2. ALINIERE LOCALĂ")
    print(f"{'─'*80}")
    print("Găsește doar cea mai bună regiune de similaritate")
    print()
    
    local_align = pairwise2.align.localms(
        seq1, seq2,
        match=1, mismatch=-1,
        open=-2, extend=-0.5,
        one_alignment_only=True
    )
    
    if local_align:
        print(format_alignment(*local_align[0]))
        print(f"Scor: {local_align[0][2]:.2f}")
        print("✓ Găsește match-ul perfect din mijloc")
    
    # ========== SEMIGLOBAL ==========
    print(f"\n{'─'*80}")
    print("3. ALINIERE SEMIGLOBALĂ (aproximare)")
    print(f"{'─'*80}")
    print("Nu penalizează gap-urile la capete")
    print()
    
    # Biopython nu are direct semiglobal, dar putem simula
    # folosind parametri care nu penalizează gap-urile terminale
    max_score, max_pos = semiglobal_alignment_manual(seq1, seq2)
    print(f"Scor maxim: {max_score}")
    print(f"Poziție finală: {max_pos}")
    print("✓ Permite gap-uri nepanalizate la capete")


def real_world_example():
    """
    Exemplu real cu secvențe din fișier.
    """
    print(f"\n{'='*80}")
    print("EXEMPLU REAL: TP53 SEQUENCES")
    print(f"{'='*80}\n")
    
    fasta_file = "tp53_protein_multi.fasta"
    
    if not os.path.exists(fasta_file):
        print(f"⚠️  Fișierul {fasta_file} nu există")
        print("Asigură-te că fișierul este în același director cu scriptul.")
        return
    
    records = list(SeqIO.parse(fasta_file, "fasta"))
    
    if len(records) < 2:
        return
    
    # Folosim primele 100 caractere pentru demonstrație
    seq1 = str(records[0].seq[:100])
    seq2 = str(records[1].seq[:100])
    
    print(f"Comparăm primele 100 aa din:")
    print(f"  - {records[0].id}")
    print(f"  - {records[1].id}")
    print()
    
    # Global
    global_align = pairwise2.align.globalms(
        seq1, seq2,
        match=2, mismatch=-1,
        open=-2, extend=-0.5,
        one_alignment_only=True
    )
    
    # Local
    local_align = pairwise2.align.localms(
        seq1, seq2,
        match=2, mismatch=-1,
        open=-2, extend=-0.5,
        one_alignment_only=True
    )
    
    print(f"Scoruri comparative:")
    if global_align:
        print(f"  • Global:      {global_align[0][2]:>8.2f}")
    if local_align:
        print(f"  • Local:       {local_align[0][2]:>8.2f}")
    
    max_score, _ = semiglobal_alignment_manual(seq1, seq2)
    print(f"  • Semiglobal:  {max_score:>8.2f}")


def print_use_cases():
    """
    Afișează cazurile de utilizare pentru fiecare tip de aliniere.
    """
    print(f"\n{'='*80}")
    print("📚 CÂND SĂ FOLOSEȘTI FIECARE TIP DE ALINIERE")
    print(f"{'='*80}\n")
    
    print("🔹 GLOBAL (Needleman-Wunsch):")
    print("   Când: secvențe de lungimi similare, complet înrudite")
    print("   Exemple:")
    print("     • Comparare gene ortologe între specii apropiate")
    print("     • Analiza variantelor allelice ale aceluiași gene")
    print("     • Secvențe proteice din aceeași familie")
    print()
    
    print("🔹 LOCAL (Smith-Waterman):")
    print("   Când: căutăm regiuni similare în secvențe neînrudite")
    print("   Exemple:")
    print("     • Identificare domenii conservate")
    print("     • Căutare motive funcționale")
    print("     • Comparare proteine cu arhitecturi diferite")
    print("     • Database search (BLAST-like)")
    print()
    
    print("🔹 SEMIGLOBAL:")
    print("   Când: o secvență este subsecvență a alteia")
    print("   Exemple:")
    print("     • Aliniere read-uri la genom")
    print("     • Căutare prime/adapter în secvențe NGS")
    print("     • Comparare gene parțiale vs complete")
    print("     • mRNA vs genom (fără introni)")
    print("     • Proteine cu domenii adăugate/șterse")
    print()
    
    print("💡 RECOMANDARE PRACTICĂ:")
    print("   • Începe cu GLOBAL dacă lungimile sunt similare")
    print("   • Folosește LOCAL pentru secvențe foarte diferite")
    print("   • Alege SEMIGLOBAL când lungimile diferă mult și")
    print("     suspectezi că o secvență conține cealaltă")


def main():
    print(f"{'='*80}")
    print("BONUS TASK: SEMIGLOBAL ALIGNMENT")
    print(f"{'='*80}\n")
    
    # Demonstrație
    demonstrate_semiglobal()
    
    # Exemplu real
    real_world_example()
    
    # Cazuri de utilizare
    print_use_cases()
    
    print(f"\n{'='*80}")
    print("💡 PENTRU NOTES.PDF (BONUS):")
    print(f"{'='*80}\n")
    
    print("Preferă SEMIGLOBAL vs GLOBAL/LOCAL când:")
    print()
    print("1. Lungimile diferă semnificativ:")
    print("   • O secvență are 100 aa, cealaltă 500 aa")
    print("   • Nu vrei să penalizezi diferența de lungime")
    print()
    print("2. Suspectezi conținere:")
    print("   • O secvență este un fragment din cealaltă")
    print("   • Ex: exon vs gene complet, domeniu vs proteină")
    print()
    print("3. Aliniezi la o referință:")
    print("   • Read-uri scurte la genom lung")
    print("   • Secvențe parțiale la database complete")
    print()
    print("4. Vrei să ignori capetele:")
    print("   • Calitatea scăzută la capete (NGS)")
    print("   • Diferențe cunoscute în regiuni terminale")


if __name__ == "__main__":
    main()
