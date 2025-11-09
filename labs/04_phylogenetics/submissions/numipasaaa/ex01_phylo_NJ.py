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
from Bio import Phylo, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor

if __name__ == "__main__":
    handle = "numipasaaa" 

    # TODO 1: Încărcați fișierul multi-FASTA propriu
    fasta = Path(f"data/work/{handle}/lab04/your_sequences.fasta")
    sequences = list(SeqIO.parse(fasta, "fasta"))

    max_len = max(len(record.seq) for record in sequences)
    aligned_sequences = []

    for record in sequences:
         padding = str(record.seq) + '-' * (max_len - len(record.seq))
         aligned_record = SeqRecord(seq=Seq(padding), id=record.id, description=record.description)
         aligned_sequences.append(aligned_record)

    alignment = MultipleSeqAlignment(aligned_sequences)

    # TODO 2: Calculați matricea de distanțe
    dist_calculator = DistanceCalculator('identity')
    distance_matrix = dist_calculator.get_distance(alignment)

    # TODO 3: Construiți arborele NJ
    constructor = DistanceTreeConstructor()
    nj_tree = constructor.nj(distance_matrix)

    # TODO 4: Salvați arborele în format Newick
    Phylo.write(nj_tree, f"labs/04_phylogenetics/submissions/{handle}/tree_{handle}.nwk", "newick")

    # TODO 5 Vizualizați arborele
    Phylo.draw_ascii(nj_tree)