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

fastq_file = "data/sample/sample.fastq"  # înlocuiți cu calea fișierului vostru FASTQ

num_reads = 0
total_length = 0
total_n = 0
total_phred = 0
total_bases = 0

for record in SeqIO.parse(fastq_file, "fastq"):
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
