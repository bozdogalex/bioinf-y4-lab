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
import os

if __name__ == "__main__":
    
    handle = "coroiusamy"

    # TODO 1: Încărcați fișierul multi-FASTA propriu
    fasta_path = Path("/workspaces/bioinf-y4-lab/data/work/coroiusamy/lab04/my_sequences.fasta")

    if not fasta_path.exists():
        print(f"EROARE: Fișierul FASTA nu a fost găsit la: {fasta_path}")
        print("Verifică dacă calea este corectă și dacă ai înlocuit <handle>.")
    else:
        print(f"Am încărcat secvențele din: {fasta_path}")
        # Citim alinierea
        aln = AlignIO.read(fasta_path, "fasta")

        # TODO 2: Calculați matricea de distanțe
        # Folosim 'identity' care e similar cu p-distance pentru secvențe aliniate
        calculator = DistanceCalculator('identity')
        dm = calculator.get_distance(aln)
        print("\nMatricea de distanțe calculată:")
        print(dm)

        # TODO 3: Construiți arborele NJ
        constructor = DistanceTreeConstructor(calculator)
        tree = constructor.nj(dm) 
        
        print("\nArborele Neighbor-Joining a fost construit.")

        # TODO 4: Salvați arborele în format Newick
        # Creeam directorul de submissions dacă nu există
        output_dir = Path(f"labs/04_phylogenetics/submissions/{handle}")
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Definim calea fișierului .nwk
        output_path = output_dir / f"tree_{handle}.nwk"
        
        Phylo.write(tree, output_path, "newick")
        print(f"\nArborele a fost salvat în format Newick la:")
        print(output_path)

        # TODO 5 Vizualizați arborele
        print("\nVizualizare arbore (format ASCII):")
        Phylo.draw_ascii(tree)