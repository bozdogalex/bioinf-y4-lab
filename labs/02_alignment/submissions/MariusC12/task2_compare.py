from Bio import SeqIO, pairwise2
from ex01_global_nw import needleman_wunsch
from ex02_local_sw import smith_waterman

fasta = "labs/02_alignment/submissions/MariusC12/task2_regions.fasta"
records = list(SeqIO.parse(fasta, "fasta"))

A = str(records[0].seq)
B = str(records[1].seq)

nw_a, nw_b, nw_score = needleman_wunsch(A, B)
sw_a, sw_b, sw_score = smith_waterman(A, B, match=1, mismatch=-1, gap=-2)

bio_global = pairwise2.align.globalms(A, B, 1, -1, -2, -2)[0]
bio_local = pairwise2.align.localms(A, B, 1, -1, -2, -2)[0]

print("NEEDLEMAN-WUNSCH")
print(nw_a)
print(nw_b)
print("Scor manual:", nw_score)
print("Scor Biopython:", bio_global.score)

print("\nSMITH-WATERMAN")
print(sw_a)
print(sw_b)
print("Scor manual:", sw_score)
print("Scor Biopython:", bio_local.score)
print("Lungime aliniere locala:", len(sw_a))