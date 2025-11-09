"""
Script utilitar pentru Laboratorul 4.

Acest script descarcă 3 secvențe FASTA (GenBank) specificate de pe serverul
NCBI folosind Bio.Entrez și le combină într-un singur fișier 
multi-FASTA, necesar ca input pentru exercițiul din Laboratorul 4
(ex05_phylo_tree.py).

Acesta înlocuiește scriptul de conversie FASTQ, folosind secvențe
complete, mai potrivite pentru o analiză filogenetică.

ACTUALIZARE: Folosim 3 gene 16S rRNA, care sunt scurte și comparabile,
perfecte pentru un exercițiu de aliniere și filogenetică.
"""

from pathlib import Path
from Bio import Entrez, SeqIO
import sys

# --- Configurare ---
# TODO: asigură-te că acesta este handle-ul tău corect
handle = "RAZVAN24RR"

# TODO: Setează un email valid. NCBI cere acest lucru.
Entrez.email = "student@unitate.ro" 

# ID-uri GenBank (FASTA) pentru 3 gene 16S rRNA din tulpini diferite
# de Staphylococcus aureus. Acestea sunt scurte și comparabile.
ids_to_fetch = [
    "NR_036928.1",  # S. aureus NCTC 8325 16S ribosomal RNA
    "NR_036932.1",  # S. aureus Mu50 16S ribosomal RNA
    "NR_025946.1"   # S. aureus USA300 16S ribosomal RNA
]

# --- Sfârșit Configurare ---

# Definirea căii de output (unde scriem)
out_lab_folder = Path(f"data/work/{handle}/lab04")
out_fasta = out_lab_folder / "your_sequences.fasta"

# Ne asigurăm că folderul de output (lab04) există
try:
    out_lab_folder.mkdir(parents=True, exist_ok=True)
except Exception as e:
    print(f"Nu am putut crea folderul de output: {e}")
    sys.exit(1)

print("--- Script de descărcare și combinare Multi-FASTA ---")
print(f"Se descarcă {len(ids_to_fetch)} secvențe (16S rRNA) de pe NCBI...")
print(f"ID-uri: {', '.join(ids_to_fetch)}")

try:
    # Contactăm serverul NCBI și cerem secvențele
    # db="nucleotide" -> baza de date GenBank
    # id=... -> lista de ID-uri
    # rettype="fasta" -> vrem formatul FASTA
    # retmode="text" -> vrem ca text
    fetch_handle = Entrez.efetch(
        db="nucleotide", 
        id=ids_to_fetch, 
        rettype="fasta", 
        retmode="text"
    )
    
    # Citim toate datele primite
    fasta_data = fetch_handle.read()
    fetch_handle.close()

    # Scriem datele direct în fișierul nostru de output
    with open(out_fasta, "w") as output_handle:
        output_handle.write(fasta_data)

    print("\n[OK] Descărcare și combinare finalizate.")
    print(f"Fișierul (nealiniat) a fost salvat în: {out_fasta.resolve()}")

except Exception as e:
    print(f"\nA apărut o eroare la contactarea NCBI sau la scrierea fișierului: {e}")
    print("Verifică conexiunea la internet și dacă email-ul este setat.")