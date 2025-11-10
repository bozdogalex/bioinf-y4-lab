
from pathlib import Path
from Bio import AlignIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
import os

if __name__ == "__main__":
    # TODO 1: Încărcați fișierul multi-FASTA propriu
   
    handle = "Ana-Maria-Bojan"  
    fasta_path = Path(f"/workspaces/bioinf-y4-lab/data/work/Ana-Maria-Bojan/lab03/mysequence.fasta")
    
    # Creează directorul dacă nu există
    os.makedirs(fasta_path.parent, exist_ok=True)
    
    alignment=AlignIO.read(fasta_path, "fasta")

    # TODO 2: Calculați matricea de distanțe
    print("\nCalculăm matricea de distanțe...")
    calculator = DistanceCalculator('identity')  # 'identity', 'blastn', 'translation'
    dm = calculator.get_distance(alignment)
    
    print("Matricea de distanțe:")
    print(dm)

    # TODO 3: Construiți arborele NJ
    print("\nConstruim arborele Neighbor-Joining...")
    constructor = DistanceTreeConstructor()
    tree = constructor.nj(dm)
    
    print("Arbore construit!")
    print(f"Număr de clade/ramuri: {len(tree.get_terminals())}")

    # TODO 4: Salvați arborele în format Newick
    output_nwk = Path(f"labs/04_phylogenetics/submissions/{handle}/tree_{handle}.nwk")
    Phylo.write(tree, output_nwk, "newick")
    print(f"\nArbore salvat în: {output_nwk}")

    # TODO 5: Vizualizați arborele
    print("\nVizualizare arbore (text):")
    Phylo.draw_ascii(tree)
    
    # Afișează arborele într-un format mai detaliat
    print("\nStructura arborelui:")
    for clade in tree.find_clades():
        if clade.name:
            print(f"  - {clade.name}: branch_length = {clade.branch_length}")
