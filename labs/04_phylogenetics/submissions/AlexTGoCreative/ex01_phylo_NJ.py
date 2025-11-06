"""
Exercițiul 5 — Construirea unui arbore Neighbor-Joining

Instrucțiuni (de urmat în laborator):
1. Refolosiți secvențele din laboratoarele anterioare (FASTA din Lab 2 sau FASTQ→FASTA din Lab 3).
2. Dacă aveți doar fișiere FASTA cu o singură secvență, combinați cel puțin 3 într-un fișier multi-FASTA:
3. Salvați fișierul multi-FASTA în: data/work/<handle>/lab04/your_sequences.fasta
4. Completați pașii de mai jos:
   - încărcați multi-FASTA-ul,
   - calculați matricea de distanțe,
   - construiți arborele NJ,
   - salvați rezultatul în format Newick (.nwk).
"""

from pathlib import Path
from Bio import SeqIO, AlignIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Align.Applications import ClustalOmegaCommandline
import shutil

if __name__ == "__main__":
    # TODO 1: Încărcați fișierul multi-FASTA propriu
    # Folosim secvențele TP53 protein multi-FASTA din data/sample
    source_fasta = Path("../../../../data/sample/tp53_protein_multi.fasta")
    work_dir = Path("../../../../data/work/AlexTGoCreative/lab04")
    work_dir.mkdir(parents=True, exist_ok=True)
    
    fasta = work_dir / "your_sequences.fasta"
    aligned_fasta = work_dir / "your_sequences_aligned.fasta"
    
    # Copiem fișierul multi-FASTA în directorul nostru de lucru
    shutil.copy(source_fasta, fasta)
    print(f"[OK] Am copiat {source_fasta} în {fasta}")
    
    # Verificăm câte secvențe avem
    records = list(SeqIO.parse(fasta, "fasta"))
    print(f"[OK] Fișierul conține {len(records)} secvențe TP53")
    
    # Pentru arbore filogenetic, avem nevoie de secvențe de aceeași lungime
    # Vom face padding pentru a alinia secvențele la aceeași lungime
    max_length = max(len(rec.seq) for rec in records)
    print(f"[INFO] Lungimea maximă a secvențelor: {max_length}")
    
    # Creăm un alignment simplu prin padding
    from Bio.Align import MultipleSeqAlignment
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq
    
    aligned_records = []
    for rec in records:
        # Adăugăm gap-uri (-) la sfârșitul secvențelor mai scurte
        padded_seq = str(rec.seq) + '-' * (max_length - len(rec.seq))
        aligned_rec = SeqRecord(Seq(padded_seq), id=rec.id, description=rec.description)
        aligned_records.append(aligned_rec)
    
    alignment = MultipleSeqAlignment(aligned_records)
    
    # Salvăm alinierea
    AlignIO.write(alignment, aligned_fasta, "fasta")
    print(f"[OK] Am creat alinierea cu padding în {aligned_fasta}")
    print(f"[OK] Aliniere cu {len(alignment)} secvențe de lungime {alignment.get_alignment_length()}")

    # TODO 2: Calculați matricea de distanțe
    calculator = DistanceCalculator("identity")
    dm = calculator.get_distance(alignment)
    print("[OK] Am calculat matricea de distanțe")
    print(f"Matrice de distanțe:\n{dm}")

    # TODO 3: Construiți arborele NJ
    constructor = DistanceTreeConstructor()
    tree = constructor.nj(dm)
    print("[OK] Am construit arborele Neighbor-Joining")

    # TODO 4: Salvați arborele în format Newick
    output = Path("tree_AlexTGoCreative.nwk")
    output.parent.mkdir(parents=True, exist_ok=True)
    Phylo.write(tree, output, "newick")
    print(f"[OK] Am salvat arborele în {output.resolve()}")

    # TODO 5: Vizualizați arborele
    print("\n[Vizualizare arbore ASCII:]")
    Phylo.draw_ascii(tree)
