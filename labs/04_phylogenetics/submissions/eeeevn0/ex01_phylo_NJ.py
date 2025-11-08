"""
Exercițiul 5 — Construirea unui arbore Neighbor-Joining

Instrucțiuni (de urmat în laborator):
1. Refolosiți secvențele din laboratoarele anterioare (FASTA din Lab 2 sau FASTQ→FASTA din Lab 3).
2. Dacă aveți doar fișiere FASTA cu o singură secvență, combinați cel puțin 3 într-un fișier multi-FASTA.
3. Salvați fișierul multi-FASTA în: data/work/<handle>/lab04/your_sequences.fasta
4. Completați pașii de mai jos:
   - încărcați multi-FASTA-ul,
   - calculați matricea de distanțe,
   - construiți arborele NJ,
   - salvați rezultatul în format Newick (.nwk).
"""

from pathlib import Path
from Bio import Phylo, SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
import sys


if __name__ == "__main__":
    # TODO 1: Încărcați fișierul multi-FASTA propriu
    fasta = Path("data/work/eeeevn0/lab04/my_sequences.fasta")

    if not fasta.exists():
        print(f"Fișierul {fasta} nu există.")
        sys.exit(1)

    records = list(SeqIO.parse(str(fasta), "fasta"))
    if len(records) < 3:
        print("Ai nevoie de cel puțin 3 secvențe pentru a construi un arbore NJ.")
        sys.exit(1)

    names = [r.id for r in records]
    print(f"Am încărcat {len(records)} secvențe din {fasta}")
    print("ID-uri secvențe:", ", ".join(names))

    lengths = [len(r.seq) for r in records]
    min_len = min(lengths)

    if len(set(lengths)) != 1:
        print(f"\nSecvențele au lungimi diferite: {lengths}")
        print(f"Trunchiem toate secvențele la lungimea minimă: {min_len} baze")
        for r in records:
            r.seq = r.seq[:min_len]
    else:
        print("\nToate secvențele au deja aceeași lungime.")

    # TODO 2: Calculați matricea de distanțe 
    alignment = MultipleSeqAlignment(records)
    calculator = DistanceCalculator("identity")
    dm = calculator.get_distance(alignment)

    print("\nMatricea de distanțe (identity):")
    print(dm)

    # TODO 3: Construiți arborele NJ
    constructor = DistanceTreeConstructor()
    tree = constructor.nj(dm)

    # TODO 4: Salvați arborele în format Newick 
    output_dir = Path("labs/04_phylogenetics/submissions/eeeevn0")
    output_dir.mkdir(parents=True, exist_ok=True)
    out_nwk = output_dir / "tree_eeeevn0.nwk"

    Phylo.write(tree, str(out_nwk), "newick")
    print(f"\nArbore NJ salvat în: {out_nwk}")

    # TODO 5: Vizualizați arborele
    print("\nArbore (ASCII):")
    Phylo.draw_ascii(tree)

    print("\nArbore (Newick):")
    Phylo.write(tree, sys.stdout, "newick")
    print()
