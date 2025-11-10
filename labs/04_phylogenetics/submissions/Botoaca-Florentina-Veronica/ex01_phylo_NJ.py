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
from Bio import AlignIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
import matplotlib.pyplot as plt
import os

if __name__ == "__main__":
    # TODO 1: Încărcați fișierul multi-FASTA propriu
    # Înlocuiește <handle> cu numele tău sau un identificator unic
    handle = "Botoaca-Florentina-Veronica"  # Înlocuiește cu numele tău real
    fasta_path = Path(f"/workspaces/bioinf-y4-lab/labs/01_intro&databases/data/work/Botoaca-Florentina-Veronica/lab01/vera_multi.fasta")
    
    # Creează directorul dacă nu există
    os.makedirs(fasta_path.parent, exist_ok=True)
    
    # Încarcă alinierea
    alignment = AlignIO.read(fasta_path, "fasta")
    print(f"A fost încărcat un aliniament cu {len(alignment)} secvențe")
    print(f"Lungimea aliniamentului: {alignment.get_alignment_length()}")

    # TODO 2: Calculați matricea de distanțe
    calculator = DistanceCalculator('identity')  # Poți folosi și 'blosum62', etc.
    distance_matrix = calculator.get_distance(alignment)
    print("\nMatricea de distanțe:")
    print(distance_matrix)

    # TODO 3: Construiți arborele NJ
    constructor = DistanceTreeConstructor()
    nj_tree = constructor.nj(distance_matrix)
    
    # TODO 4: Salvați arborele în format Newick
    output_nwk = Path(f"/workspaces/bioinf-y4-lab/labs/04_phylogenetics/submissions/Botoaca-Florentina-Veronica/your_tree.nwk")
    os.makedirs(output_nwk.parent, exist_ok=True)
    
    Phylo.write(nj_tree, output_nwk, "newick")
    print(f"\nArborele a fost salvat în: {output_nwk}")

    # TODO 5: Vizualizați arborele
    print("\nVizualizare arbore:")
    Phylo.draw_ascii(nj_tree)
    
    # Vizualizare grafică (opțional)
    try:
        fig = plt.figure(figsize=(10, 8))
        ax = fig.add_subplot(111)
        Phylo.draw(nj_tree, axes=ax, do_show=False)
        plt.title("Arbore Filogenetic - Neighbor Joining")
        plt.savefig(f"/workspaces/bioinf-y4-lab/labs/04_phylogenetics/submissions/Botoaca-Florentina-Veronica/tree_visualization.png")
        plt.show()
        print(f"\nVizualizarea grafică a fost salvată ca PNG")
    except Exception as e:
        print(f"Vizualizarea grafică a eșuat: {e}")
        print("Folosind doar vizualizarea ASCII")