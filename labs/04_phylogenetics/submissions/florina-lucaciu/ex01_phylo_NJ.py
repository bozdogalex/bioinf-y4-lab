from pathlib import Path
from Bio import AlignIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor

if __name__ == "__main__":
    # TODO 1: Încărcați fișierul multi-FASTA propriu
    fasta = Path("data/work/florina-lucaciu/lab04/your_sequences.fasta")


    alignment = AlignIO.read(fasta, "fastq")
    print(f"Aliniament incarcat din: {fasta}")
    print(f"Nr. secvente: {len(alignment)}")

    # TODO 2: Calculați matricea de distanțe
    calculator = DistanceCalculator("identity")
    dm = calculator.get_distance(alignment)
    print("\nMatricea de distante calculata:")
    print(dm)

    # TODO 3: Construiți arborele NJ
    constructor = DistanceTreeConstructor()
    tree = constructor.nj(dm)
    print("\nArborele Neighbor-Joining construit.")

    # TODO 4: Salvați arborele în format Newick
    output = Path("labs/04_phylogenetics/submissions/florina-lucaciu/tree_florina-lucaciu.nwk")
    # Se creeaza directorul daca nu exista
    output.parent.mkdir(parents=True, exist_ok=True) 
    
    Phylo.write(tree, output, "newick")
    print(f"Arborele a fost salvat în format Newick la: {output}")

    # TODO 5 Vizualizați arborele
    print("\nVizualizare arbore (format ASCII):")
    Phylo.draw_ascii(tree)