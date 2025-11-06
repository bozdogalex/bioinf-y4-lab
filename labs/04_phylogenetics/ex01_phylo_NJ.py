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
from Bio import AlignIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor, DistanceMatrix
from Bio import SeqIO, pairwise2
from Bio.Align import MultipleSeqAlignment
import sys

if __name__ == "__main__":
    # TODO 1: Încărcați fișierul multi-FASTA propriu
    fasta = Path("data/work/eeeevn0/lab04/your_sequences.fasta")
    if not fasta.exists():
        print(f"Fișierul {fasta} nu există. Înlocuiți <handle> și asigurați-vă că fișierul există.")
        sys.exit(1)

    # Încărcăm secvențele ca SeqRecord
    records = list(SeqIO.parse(str(fasta), "fasta"))
    if len(records) < 3:
        print("Aveți nevoie de cel puțin 3 secvențe pentru a construi un arbore NJ.")
        sys.exit(1)



    # TODO 2: Calculați matricea de distanțe
    names = [r.id for r in records]

    # Dacă secvențele sunt deja aliniate (aceeași lungime) folosim DistanceCalculator
    same_length = all(len(r.seq) == len(records[0].seq) for r in records)
    if same_length:
        alignment = MultipleSeqAlignment(records)
        calculator = DistanceCalculator('identity')
        dm = calculator.get_distance(alignment)
    else:
        # Dacă nu sunt aliniate, calculează distanțe pereche cu pairwise globalxx (identitate)
        n = len(records)
        full = [[0.0]*n for _ in range(n)]
        for i in range(n):
            seq_i = str(records[i].seq)
            for j in range(i+1, n):
                seq_j = str(records[j].seq)
                aln = pairwise2.align.globalxx(seq_i, seq_j, one_alignment_only=True)[0]
                aln_seqA, aln_seqB, score = aln[0], aln[1], aln[2]
                aln_len = len(aln_seqA)
                # distanța = 1 - (număr potriviri / lungimea aliniamentului)
                dist = 1.0 - (score / aln_len) if aln_len > 0 else 1.0
                full[i][j] = full[j][i] = dist
        # Construim matricea în formatul așteptat (lower-triangular)
        lower = []
        for i in range(n):
            row = [full[i][j] for j in range(i)]
            lower.append(row)
        dm = DistanceMatrix(names, lower)


    # TODO 3: Construiți arborele NJ
    constructor = DistanceTreeConstructor()
    tree = constructor.nj(dm)

    # TODO 4: Salvați arborele în format Newick
    out_nwk = fasta.with_suffix(".nwk")
    Phylo.write(tree, str(out_nwk), "newick")
    print(f"Arbore salvat în: {out_nwk}")

    # TODO 5 Vizualizați arborele
    print("\nArbore (ASCII):")
    Phylo.draw_ascii(tree)

    try:
        Phylo.draw(tree)  # va deschide o fereastră matplotlib dacă e disponibil
    except Exception:
        print("matplotlib nu este disponibil sau afișarea grafică a eșuat; am afișat arborele ASCII.")
