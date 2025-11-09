
from pathlib import Path
from Bio import SeqIO, AlignIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import sys


def load_sequences(fasta_path: Path):
    """Load sequences from multi-FASTA file."""
    if not fasta_path.exists():
        raise FileNotFoundError(f"Fișierul nu există: {fasta_path}")
    
    records = list(SeqIO.parse(fasta_path, "fasta"))
    print(f"Loaded {len(records)} sequences:")
    for rec in records:
        print(f"  - {rec.id}: {len(rec.seq)} bp")
    
    return records


def calculate_distance_matrix(records):
    """
    Calculate p-distance matrix from sequences.
    p-distance = proportion of different positions
    """
    import numpy as np
    
    n = len(records)
    sequences = [str(rec.seq) for rec in records]
    
    # Find minimum length (to avoid index errors)
    min_len = min(len(seq) for seq in sequences)
    print(f"\nMinimum sequence length: {min_len} bp")
    print("Using first {min_len} bp for distance calculation...")
    
    # Truncate all sequences to same length
    truncated = [seq[:min_len] for seq in sequences]
    
    # Calculate pairwise p-distances
    matrix = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            # Count differences
            differences = sum(a != b for a, b in zip(truncated[i], truncated[j]))
            # Calculate proportion
            p_dist = differences / min_len
            matrix[i, j] = matrix[j, i] = p_dist
    
    print("\nDistance Matrix (p-distance):")
    print("Sequences:", [rec.id for rec in records])
    print(matrix)
    
    return matrix, truncated


def create_alignment_from_sequences(records, truncated_seqs):
    """
    Create a MultipleSeqAlignment object for Biopython Phylo.
    This pads sequences to the same length if needed.
    """
    # All sequences already same length (truncated)
    max_len = len(truncated_seqs[0])
    
    aligned_records = []
    for rec, seq in zip(records, truncated_seqs):
        # Pad if needed (shouldn't be necessary after truncation)
        padded = seq.ljust(max_len, '-')
        aligned_rec = SeqRecord(
            Seq(padded),
            id=rec.id,
            description=rec.description
        )
        aligned_records.append(aligned_rec)
    
    alignment = MultipleSeqAlignment(aligned_records)
    return alignment


def build_nj_tree(alignment):
    """
    Build Neighbor-Joining tree using Biopython.
    """
    print("\nCalculating phylogenetic distances...")
    calculator = DistanceCalculator('identity')
    dm = calculator.get_distance(alignment)
    
    print("\nConstructing Neighbor-Joining tree...")
    constructor = DistanceTreeConstructor(calculator, 'nj')
    tree = constructor.build_tree(alignment)
    
    return tree


def save_tree(tree, output_path: Path):
    """Save tree in Newick format."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    Phylo.write(tree, output_path, 'newick')
    print(f"\n✓ Tree saved to: {output_path}")


def main():
    # Handle for GitHub username
    handle = "IrisDanila"
    
    # Input: TP53 mRNA sequences from Lab 1
    fasta_path = Path(f"data/work/{handle}/lab04/your_sequences.fasta")
    
    # Output: Newick tree file
    output_tree = Path(f"labs/04_phylogenetics/submissions/{handle}/tree_{handle}.nwk")
    
    print("="*70)
    print("LAB 4 - Phylogenetic Tree Construction (Neighbor-Joining)")
    print("="*70)
    
    # Step 1: Load sequences
    print("\n[Step 1] Loading sequences...")
    records = load_sequences(fasta_path)
    
    if len(records) < 3:
        print("\n⚠ WARNING: Need at least 3 sequences for meaningful tree!")
        print("You have only", len(records), "sequences.")
        sys.exit(1)
    
    # Step 2: Calculate distance matrix
    print("\n[Step 2] Calculating distance matrix...")
    matrix, truncated = calculate_distance_matrix(records)
    
    # Step 3: Create alignment for Biopython
    print("\n[Step 3] Creating alignment...")
    alignment = create_alignment_from_sequences(records, truncated)
    
    # Step 4: Build NJ tree
    print("\n[Step 4] Building Neighbor-Joining tree...")
    tree = build_nj_tree(alignment)
    
    # Step 5: Save tree
    print("\n[Step 5] Saving tree...")
    save_tree(tree, output_tree)
    
    print("\n" + "="*70)
    print("TREE CONSTRUCTION COMPLETE!")
    print("="*70)
    print(f"\nTree file: {output_tree}")
    print("="*70)


if __name__ == "__main__":
    main()
