from Bio import SeqIO, Align
import itertools
from Bio.Seq import Seq
from pathlib import Path
from Bio import pairwise2

#task 1

def p_distance_aligned(aln1, aln2):
    assert len(aln1) == len(aln2)
    diff = 0
    valid = 0
    for a, b in zip(aln1, aln2):
        if a == "-" or b == "-":
            continue
        valid += 1
        if a != b:
            diff += 1
    return diff / valid if valid > 0 else 0.0


def read_fasta(filepath, max_seqs: int = 5):
    seqs = {}
    for i, record in enumerate(SeqIO.parse(str(filepath), "fasta")):
        if i >= max_seqs:
            break
        seqs[record.id] = str(record.seq)
    return seqs


def compute_distance_matrix(fasta_file, max_seqs: int = 5):
    sequences = read_fasta(fasta_file, max_seqs=max_seqs)
    names = list(sequences.keys())
    n = len(names)

    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 1
    aligner.mismatch_score = -1
    aligner.open_gap_score = -1
    aligner.extend_gap_score = -0.5

    matrix = [[0.0 for _ in range(n)] for _ in range(n)]

    for i, j in itertools.combinations(range(n), 2):
        s1 = sequences[names[i]]
        s2 = sequences[names[j]]

        alignment = aligner.align(s1, s2)[0]

        aln1 = str(alignment[0])
        aln2 = str(alignment[1])

        pd = p_distance_aligned(aln1, aln2)
        matrix[i][j] = pd
        matrix[j][i] = pd

    return names, matrix

def print_matrix(names, matrix):
    print("\t" + "\t".join(names))
    for name, row in zip(names, matrix):
        print(name + "\t" + "\t".join(f"{v:.4f}" for v in row))

#task 2

def load_two_sequences(fasta_path: Path, i1: int, i2: int):
    recs = list(SeqIO.parse(str(fasta_path), "fasta"))
    if len(recs) < 2:
        raise SystemExit("[eroare] Fișierul trebuie să conțină cel puțin 2 secvențe.")
    if not (0 <= i1 < len(recs) and 0 <= i2 < len(recs)):
        raise SystemExit(f"[eroare] Indici invalizi (0..{len(recs)-1}).")
    return str(recs[i1].seq), str(recs[i2].seq), recs[i1].id, recs[i2].id

if __name__ == "__main__":
    fasta_path = "/workspaces/bioinf-y4-lab/data/work/blejdeaxoxo/lab01/mito_trna.fa"
    #task 1
    names, matrix = compute_distance_matrix(fasta_path, max_seqs=5)
    print_matrix(names, matrix)

    #task 2
    seq1, seq2, id1, id2 = load_two_sequences(fasta_path, 1, 3)
    global_alignments = pairwise2.align.globalxx(seq1, seq2)

    print("=== GLOBAL ALIGNMENT ===")
    for aln in global_alignments[:1]:
        print(pairwise2.format_alignment(*aln))

    local_alignments = pairwise2.align.localxx(seq1, seq2)

    print("=== LOCAL ALIGNMENT ===")
    for aln in local_alignments[:1]:
        print(pairwise2.format_alignment(*aln))
