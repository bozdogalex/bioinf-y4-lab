"""
Exercițiul 5 — Construirea unui arbore Neighbor-Joining

Instrucțiuni (de urmat în laborator):
1. Refolosiți secvențele din laboratoarele anterioare (FASTA din Lab 2 sau FASTQ→FASTA din Lab 3).
2. Dacă aveți doar fișiere FASTA cu o singură secvență, combinați cel puțin 3 într-un fișier multi-FASTA:
3. Salvați fișierul multi-FASTA în: data/sample/tp53_dna_multi.fasta
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
    # Încărcați fișierul multi-FASTA propriu
    fasta = Path("labs/04_phylogenetics/submissions/banmepls/sequence2.aln-fasta")

    # Exemplu (decomentați după ce înlocuiți <handle>):
    alignment = AlignIO.read(fasta, "fasta")

    # Calculați matricea de distanțe
    calculator = DistanceCalculator('identity')
    dm = calculator.get_distance(alignment)

    # Construiți arborele NJ
    constructor = DistanceTreeConstructor()
    tree = constructor.nj(dm)

    # Salvați arborele în format Newick
    Phylo.write(tree, "tree.nwk", "newick")

    # Vizualizați arborele
    Phylo.draw_ascii(tree)
