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

if __name__ == "__main__":
    fasta_file = Path("data/work/blejdeaxoxo/lab04/your_seq.fasta")
    
    if not fasta_file.exists():
        raise FileNotFoundError(f"Fișierul {fasta_file} nu a fost găsit.")

    alignment = AlignIO.read(fasta_file, "fasta")
    
    calculator = DistanceCalculator('identity')  
    distance_matrix = calculator.get_distance(alignment)
    print("Matricea de distanțe:")
    print(distance_matrix)

    constructor = DistanceTreeConstructor()
    nj_tree = constructor.nj(distance_matrix)

    output_file = Path("labs/04_phylogenetics/submissions/blejdeaxoxo/nj_tree.nwk")
    Phylo.write(nj_tree, output_file, "newick")
    print(f"Arborele NJ a fost salvat în {output_file}")

    Phylo.draw_ascii(nj_tree)