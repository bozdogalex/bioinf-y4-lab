"""
Exercițiul 5 — Construirea unui arbore Neighbor-Joining

Pași:
- încărcați multi-FASTA-ul,
- calculați matricea de distanțe,
- construiți arborele NJ,
- salvați rezultatul în format Newick (.nwk).
"""

from pathlib import Path
from Bio import AlignIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
import matplotlib.pyplot as plt

HANDLE = "olezhka0809"

if __name__ == "__main__":
    # TODO 1: Încărcați fișierul multi-FASTA propriu
    # Varianta simplă: fișierul este în același director cu scriptul
    fasta = Path("/workspaces/bioinf-y4-lab/data/work/olezhka0809/lab04/your_sequences_aligned.fasta")

    if not fasta.exists():
        raise FileNotFoundError(f"Nu am găsit {fasta.resolve()}")

    # Dacă fișierul este aliniament multiplu (MSA) în format FASTA:
    alignment = AlignIO.read(str(fasta), "fasta")
    print(f"[INFO] Am încărcat {len(alignment)} secvențe din {fasta.name}")

    # TODO 2: Calculați matricea de distanțe
    # Pentru proteine folosim de obicei un model de tip BLOSUM62
    calculator = DistanceCalculator("blosum62")
    distance_matrix = calculator.get_distance(alignment)
    print("[INFO] Matrice de distanțe calculată:")
    print(distance_matrix)

    # TODO 3: Construiți arborele NJ
    constructor = DistanceTreeConstructor(calculator, "nj")
    nj_tree = constructor.build_tree(alignment)
    print("[INFO] Arbore NJ construit.")

    # TODO 4: Salvați arborele în format Newick
    out_nwk = Path(f"tree_{HANDLE}.nwk")
    Phylo.write(nj_tree, str(out_nwk), "newick")
    print(f"[OK] Arbore NJ salvat în: {out_nwk.resolve()}")

    # TODO 5 (opțional): Vizualizați arborele (ASCII în terminal)
    print("\nArbore NJ (ASCII):")
    Phylo.draw_ascii(nj_tree)

    

