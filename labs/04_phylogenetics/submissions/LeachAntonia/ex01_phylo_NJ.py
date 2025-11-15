from pathlib import Path
from Bio import AlignIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio import SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os

if __name__ == "__main__":
    # TODO 1: Încărcați fișierul multi-FASTA propriu
    handle = "LeachAntonia"
    fasta = Path(f"/workspaces/bioinf-y4-lab/data/work/{handle}/lab04/tp53_sequences.fasta")
    aligned_fasta = Path(f"/workspaces/bioinf-y4-lab/data/work/{handle}/lab04/tp53_aligned_sequences.fasta")

    if not fasta.exists():
        raise FileNotFoundError(f"FASTA file not found: {fasta}")

    records = list(SeqIO.parse(fasta, "fasta"))
    print(f"[OK] Fișierul conține {len(records)} secvențe")

    max_length = max(len(rec.seq) for rec in records)
    print(f"[INFO] Lungimea maximă a secvențelor: {max_length}")

    aligned_records = []
    for rec in records:
        padded_seq = str(rec.seq) + '-' * (max_length - len(rec.seq))
        aligned_rec = SeqRecord(Seq(padded_seq), id=rec.id, description=rec.description)
        aligned_records.append(aligned_rec)

    alignment = MultipleSeqAlignment(aligned_records)

    AlignIO.write(alignment, aligned_fasta, "fasta")
    print(f"[OK] Aliniere cu padding salvată în {aligned_fasta}")
    print(f"[OK] Aliniere cu {len(alignment)} secvențe de lungime {alignment.get_alignment_length()}")

    # TODO 2: Calculați matricea de distanțe
    calculator = DistanceCalculator("identity")
    distance_matrix = calculator.get_distance(alignment)
    print("\nDistance matrix:")
    print(distance_matrix)

    # TODO 3: Construiți arborele NJ
    constructor = DistanceTreeConstructor()
    nj_tree = constructor.nj(distance_matrix)
    print("\nNJ tree constructed successfully!")

    # TODO 4: Salvați arborele în format Newick
    out_dir = Path(f"/workspaces/bioinf-y4-lab/labs/04_phylogenetics/submissions/{handle}")
    out_dir.mkdir(parents=True, exist_ok=True)
    newick_path = out_dir / f"tree_{handle}.nwk"
    Phylo.write(nj_tree, newick_path, "newick")
    print(f"Tree saved to: {newick_path}")
   
    # TODO 5 Vizualizați arborele
    print("\nASCII view of NJ tree:")
    Phylo.draw_ascii(nj_tree)