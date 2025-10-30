from pathlib import Path
from Bio import AlignIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor

if __name__ == "__main__":
    handle = "MarioCld"

    fasta = Path(f"data/work/MarioCld/lab04/sequence.fasta")
    alignment = AlignIO.read(fasta, "fasta")

    calculator = DistanceCalculator("identity")
    dm = calculator.get_distance(alignment)
    print("\nMatricea de distanțe:")
    print(dm)

    constructor = DistanceTreeConstructor()
    tree = constructor.nj(dm)

    output = Path(f"labs/04_phylogenetics/submissions/MarioCld/tree_MarioCld.nwk")
    output.parent.mkdir(parents=True, exist_ok=True)
    Phylo.write(tree, output, "newick")
    print(f"\nArborele a fost salvat în: {output}")

    Phylo.draw_ascii(tree)
