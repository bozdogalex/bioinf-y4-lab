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
import argparse
import sys
import shutil

from Bio import Phylo, SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor

if __name__ == "__main__":
   # TODO 1: Încărcați fișierul multi-FASTA propriu
    src = Path("data/work/marabonea16/lab01/my_tp53.fa")
    dest = Path("data/work/marabonea16/lab04/your_sequences.fasta")
    dest.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(str(src), str(dest))
   
    records = list(SeqIO.parse(str(dest), "fasta"))
    if len(records) < 3:
        sys.exit(3)

    lengths = [len(r.seq) for r in records]
    if len(set(lengths)) == 1:
        alignment = MultipleSeqAlignment(records)
    else:
        min_len = min(lengths)
        trimmed = []
        for r in records:
            new_r = r[:min_len]
            trimmed.append(new_r)
        alignment = MultipleSeqAlignment(trimmed)

    # TODO 2: Calculați matricea de distanțe
    calculator = DistanceCalculator("identity")
    dm = calculator.get_distance(alignment)

    # TODO 3: Construiți arborele NJ
    constructor = DistanceTreeConstructor()
    tree = constructor.nj(dm)

    # TODO 4: Salvați arborele în format Newick
    output = Path("labs/04_phylogenetics/submissions/marabonea16/tree_marabonea16.nwk")
    output.parent.mkdir(parents=True, exist_ok=True)
    Phylo.write(tree, output, "newick")

    # TODO 5 Vizualizați arborele
    Phylo.draw_ascii(tree)
