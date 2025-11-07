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
    # TODO 1: Încărcați fișierul multi-FASTA propriu
    fasta = Path("data/work/<handle>/lab04/your_sequences.fasta")

    # Exemplu (decomentați după ce înlocuiți <handle>):

    # TODO 2: Calculați matricea de distanțe

    # TODO 3: Construiți arborele NJ

    # TODO 4: Salvați arborele în format Newick

    # TODO 5 Vizualizați arborele
