"""
Lab 04 - Neighbor-Joining tree construction for handle Papaia462819.

Steps:
1. Load multi-FASTA alignment from data/work/<handle>/lab04/your_sequences.fasta
2. Compute pairwise distances (identity distance) with Biopython
3. Build Neighbor-Joining tree
4. Save tree to labs/04_phylogenetics/submissions/<handle>/tree_<handle>.nwk
"""

from pathlib import Path

from Bio import AlignIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor


HANDLE = "Papaia462819"


def main() -> None:
    fasta_path = Path("data/work") / HANDLE / "lab04" / "your_sequences.fasta"
    output_dir = Path("labs/04_phylogenetics/submissions") / HANDLE
    tree_path = output_dir / f"tree_{HANDLE}.nwk"

    if not fasta_path.exists():
        raise FileNotFoundError(
            f"Expected multi-FASTA file at {fasta_path}. "
            "Create it with at least three aligned sequences."
        )

    alignment = AlignIO.read(fasta_path, "fasta")

    # Compute the pairwise distance matrix with the simple identity (p-distance) metric.
    calculator = DistanceCalculator("identity")
    distance_matrix = calculator.get_distance(alignment)

    constructor = DistanceTreeConstructor()
    tree = constructor.nj(distance_matrix)

    output_dir.mkdir(parents=True, exist_ok=True)
    Phylo.write(tree, tree_path, "newick")

    print(f"Loaded {len(alignment)} sequences from {fasta_path}")
    print("Distance matrix:")
    print(distance_matrix)
    print(f"Neighbor-Joining tree saved to {tree_path}")

    # Print an ASCII representation to quickly inspect topology in the terminal.
    Phylo.draw_ascii(tree)


if __name__ == "__main__":
    main()
