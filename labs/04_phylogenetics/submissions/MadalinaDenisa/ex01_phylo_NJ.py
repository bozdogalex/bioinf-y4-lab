from pathlib import Path
from Bio import AlignIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor

# --- Config ---
handle = "MadalinaDenisa"
fasta = Path(f"data/work/{handle}/lab04/your_sequences_aligned.fasta")
output = Path(f"labs/04_phylogenetics/submissions/{handle}/tree_{handle}.nwk")


alignment = AlignIO.read(fasta, "fasta")
print(f"Loaded alignment with {len(alignment)} sequences")

# --- Calculează matricea de distanțe ---
calculator = DistanceCalculator("identity")
dm = calculator.get_distance(alignment)

# --- Construiește arborele NJ ---
constructor = DistanceTreeConstructor()
tree = constructor.nj(dm)

# --- Salvează arborele în format Newick ---
output.parent.mkdir(parents=True, exist_ok=True)
Phylo.write(tree, output, "newick")
print(f"Arborele a fost salvat în {output}")

# --- Vizualizează arborele în ASCII ---
Phylo.draw_ascii(tree)

