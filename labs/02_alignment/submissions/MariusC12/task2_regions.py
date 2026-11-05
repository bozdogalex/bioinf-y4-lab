from Bio import SeqIO
from Bio.Align import PairwiseAligner

human = SeqIO.read("data/work/MariusC12/lab01/seq1.fasta", "fasta")
mouse = SeqIO.read("data/work/MariusC12/lab01/seq2.fasta", "fasta")

human_cds = str(human.seq)[142:1324]
mouse_cds = str(mouse.seq)[157:1330]

aligner = PairwiseAligner()
aligner.mode = "global"
aligner.match_score = 1
aligner.mismatch_score = -1
aligner.open_gap_score = -2
aligner.extend_gap_score = -2

alignment = aligner.align(human_cds, mouse_cds)[0]

print("Scor aliniere:", alignment.score)
print("\nBlocuri aliniate de minimum 50 nucleotide:")

human_blocks, mouse_blocks = alignment.aligned

for h, m in zip(human_blocks, mouse_blocks):
    length = min(h[1] - h[0], m[1] - m[0])

    if length >= 50:
        print(
            f"Human CDS [{h[0]}:{h[0] + length}] "
            f"| Mouse CDS [{m[0]}:{m[0] + length}] "
            f"| Lungime: {length}"
        )



human_region = human_cds[521:621]
mouse_region = mouse_cds[512:612]

print("\nRegiuni selectate:")
print("Human:", human_region)
print("Mouse:", mouse_region)
print("Human length:", len(human_region))
print("Mouse length:", len(mouse_region))

with open("labs/02_alignment/submissions/MariusC12/task2_regions.fasta", "w") as f:
    f.write(">Human_NM_000546.6_coordinates_663_763\n")
    f.write(human_region + "\n")
    f.write(">Mouse_NM_011640.3_coordinates_669_769\n")
    f.write(mouse_region + "\n")