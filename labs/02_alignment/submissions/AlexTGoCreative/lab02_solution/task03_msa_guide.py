#!/usr/bin/env python3
"""
Task 3: Multiple Sequence Alignment (MSA) - Guide and Analysis
Ghid pentru utilizarea Clustal Omega și analiza rezultatelor MSA.
"""

from Bio import SeqIO
import os


def prepare_fasta_for_msa(input_file, output_file, num_sequences=3):
    """
    Pregătește un fișier FASTA cu un subset de secvențe pentru MSA online.
    """
    records = list(SeqIO.parse(input_file, "fasta"))
    
    # Selectează primele N secvențe
    selected = records[:num_sequences]
    
    # Salvează în format FASTA
    SeqIO.write(selected, output_file, "fasta")
    
    print(f"✓ Creat {output_file} cu {len(selected)} secvențe pentru MSA")
    print(f"\nSecvențe incluse:")
    for i, rec in enumerate(selected, 1):
        print(f"  {i}. {rec.id} ({len(rec.seq)} caractere)")
    
    return selected


def print_msa_instructions():
    """
    Afișează instrucțiuni pentru utilizarea Clustal Omega online.
    """
    print(f"\n{'='*80}")
    print("GHID PENTRU MSA CU CLUSTAL OMEGA (EBI)")
    print(f"{'='*80}\n")
    
    print("📝 PAȘI PENTRU RULAREA MSA ONLINE:")
    print()
    
    print("1. Accesează Clustal Omega pe EBI:")
    print("   🔗 https://www.ebi.ac.uk/Tools/msa/clustalo/")
    print()
    
    print("2. Pregătire secvențe:")
    print("   • Deschide fișierul 'sequences_for_msa.fasta' generat")
    print("   • Copiază conținutul")
    print()
    
    print("3. Configurare tool:")
    print("   • Lipește secvențele în câmpul de text")
    print("   • Sau încarcă fișierul direct (Upload file)")
    print("   • Lasă setările default (sunt optime pentru majoritatea cazurilor)")
    print()
    
    print("4. Rulare:")
    print("   • Click pe butonul 'Submit'")
    print("   • Așteaptă procesarea (de obicei < 1 minut)")
    print()
    
    print("5. Descărcare rezultate:")
    print("   • Click pe 'Download Alignment File' (format ClustalW sau FASTA)")
    print("   • Salvează ca 'msa_result.aln' în directorul curent")
    print()
    
    print("6. Vizualizare:")
    print("   • Rezultatul afișat online arată:")
    print("     * - pozițiile identice în toate secvențele")
    print("     : - substituții conservative")
    print("     . - substituții semi-conservative")
    print()
    
    print(f"{'─'*80}")
    print("🔍 CE SĂ CAUȚI ÎN REZULTATUL MSA:")
    print(f"{'─'*80}\n")
    
    print("1. REGIUNI CONSERVATE (multe '*'):")
    print("   • Indică zone funcționale importante")
    print("   • Site-uri active, domenii structurale")
    print("   • Presiune evolutivă pentru conservare")
    print()
    
    print("2. REGIUNI VARIABILE (puține match-uri):")
    print("   • Zone cu presiune evolutivă mai mică")
    print("   • Pot indica adaptări specifice speciei")
    print()
    
    print("3. PATTERN-URI DE GAP-URI:")
    print("   • Gap-uri la capete: diferențe în lungime")
    print("   • Gap-uri în mijloc: inserții/deleții evolutive")
    print()
    
    print(f"{'─'*80}")
    print("📊 AVANTAJE MSA vs PAIRWISE:")
    print(f"{'─'*80}\n")
    
    print("✓ MSA oferă:")
    print("  • Vedere de ansamblu asupra mai multor secvențe simultan")
    print("  • Identificare mai bună a regiunilor conservate")
    print("  • Context evolutiv mai bogat")
    print("  • Detecție de motive funcționale comune")
    print()
    
    print("✓ Pairwise este mai bun pentru:")
    print("  • Comparații directe între două secvențe")
    print("  • Analiza detaliată a diferențelor")
    print("  • Când vrei control precis asupra parametrilor")
    print()


def analyze_msa_result(msa_file):
    """
    Analizează un fișier MSA rezultat (dacă există).
    """
    if not os.path.exists(msa_file):
        print(f"\n⚠️  Fișierul {msa_file} nu există încă.")
        print("Rulează MSA online și salvează rezultatul aici pentru analiză automată.")
        return
    
    print(f"\n{'='*80}")
    print(f"ANALIZA REZULTATULUI MSA: {msa_file}")
    print(f"{'='*80}\n")
    
    # Citește alinierea
    try:
        from Bio import AlignIO
        alignment = AlignIO.read(msa_file, "clustal")
        
        print(f"✓ Aliniere încărcată cu succes")
        print(f"  - Număr secvențe: {len(alignment)}")
        print(f"  - Lungime aliniere: {alignment.get_alignment_length()}")
        print()
        
        # Calculează statistici de conservare
        length = alignment.get_alignment_length()
        conserved_positions = 0
        
        for i in range(length):
            column = alignment[:, i]
            if len(set(column)) == 1 and '-' not in column:
                conserved_positions += 1
        
        print(f"📊 Statistici conservare:")
        print(f"  - Poziții complet conservate: {conserved_positions}")
        print(f"  - Procent conservare: {100 * conserved_positions / length:.2f}%")
        
    except Exception as e:
        print(f"⚠️  Eroare la citirea MSA: {e}")
        print("Asigură-te că fișierul este în format Clustal (.aln)")


def main():
    print(f"{'='*80}")
    print("TASK 3: MULTIPLE SEQUENCE ALIGNMENT (MSA)")
    print(f"{'='*80}\n")
    
    # Pregătește fișierul pentru MSA
    input_file = "tp53_protein_multi.fasta"
    output_file = "sequences_for_msa.fasta"
    
    if os.path.exists(input_file):
        selected = prepare_fasta_for_msa(input_file, output_file, num_sequences=3)
        
        print(f"\n📄 Conținut {output_file}:")
        print(f"{'─'*80}")
        with open(output_file, 'r') as f:
            content = f.read()
            # Afișează primele 500 caractere
            if len(content) > 500:
                print(content[:500] + "\n... (truncat)")
            else:
                print(content)
        print(f"{'─'*80}")
    
    # Afișează instrucțiuni
    print_msa_instructions()
    
    # Încearcă să analizeze rezultatul (dacă există)
    analyze_msa_result("msa_result.aln")

if __name__ == "__main__":
    main()
