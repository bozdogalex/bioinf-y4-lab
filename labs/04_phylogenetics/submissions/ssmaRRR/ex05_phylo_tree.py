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
    # 1. Încarcă fișierul multi-FASTA (înlocuiește <handle> cu userul tău real)
    fasta = Path("data/work/ssmaRRR/lab04/your_sequences.fasta")

    if not fasta.exists():
        raise FileNotFoundError(f"Nu găsesc fișierul: {fasta}")

    # 2. Citește alinierea (toate secvențele trebuie să aibă aceeași lungime)
    alignment = AlignIO.read(fasta, "fasta")

    # 3. Calculează matricea de distanțe (model 'identity' -> distanță = 1 - identitate)
    calculator = DistanceCalculator("identity")
    distance_matrix = calculator.get_distance(alignment)

    # 4. Construieste arborele Neighbor-Joining
    constructor = DistanceTreeConstructor()
    nj_tree = constructor.nj(distance_matrix)

    # 5. Salvează arborele în format Newick (.nwk)
    newick_path = fasta.with_suffix(".nwk")
    Phylo.write(nj_tree, newick_path, "newick")
    print(f"Arborele NJ a fost salvat în: {newick_path}")

    # 6. (Optional laborator) Vizualizare text în consolă
    Phylo.draw_ascii(nj_tree)
