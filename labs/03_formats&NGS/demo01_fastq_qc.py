"""
Demo 01 — FASTQ Quality Control (QC)

Acest script demonstrează cum putem citi un fișier FASTQ și calcula statistici simple:
- număr total de citiri
- lungimea medie a citirilor
- proporția bazelor 'N'
- scorul Phred mediu aproximativ

Date de intrare: un fișier FASTQ mic din data/sample/ (ex. sample.fastq)
"""

from Bio import SeqIO
import gzip

fastq_file = "/workspaces/bioinf-y4-lab/data/sample/CPCT12345678R_AHHKYHDSXX_S13_L001_R1_001.fastq.gz"

num_reads = 0
total_length = 0
total_n = 0
total_phred = 0
total_bases = 0
with gzip.open(fastq_file, "rt") as handle:
    for record in SeqIO.parse(handle, "fastq"):
        num_reads += 1
        seq = str(record.seq)
        total_length += len(seq)
        total_n += seq.count("N")
        # scoruri de calitate
        phred_scores = record.letter_annotations["phred_quality"]
        total_phred += sum(phred_scores)
        total_bases += len(phred_scores)

len_mean = total_length / num_reads if num_reads > 0 else 0
n_rate = total_n / total_length if total_length > 0 else 0
phred_mean = total_phred / total_bases if total_bases > 0 else 0

print(f"Reads: {num_reads}")
print(f"Mean length: {len_mean:.2f}")
print(f"N rate: {n_rate:.4f}")
print(f"Mean Phred: {phred_mean:.2f}")