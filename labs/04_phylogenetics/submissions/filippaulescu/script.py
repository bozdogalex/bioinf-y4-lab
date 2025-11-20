from pathlib import Path
from Bio import SeqIO
import gzip

# 1. Fisierul de intrare din Lab 3 (FASTQ.GZ)
input_fastq_gz = Path("data/work/filippaulescu/lab03/your_reads.fastq.gz")  # <-- schimba daca e in alta parte

# 2. Fisierul de iesire pentru Lab 4 (multi-FASTA)
output_fasta = Path("data/work/filippaulescu/lab04/my_sequences.fasta")
output_fasta.parent.mkdir(parents=True, exist_ok=True)

records = []

# Citim din FASTQ.GZ si luam primele 3 secvente
with gzip.open(input_fastq_gz, "rt") as handle:
    for i, rec in enumerate(SeqIO.parse(handle, "fastq")):
        records.append(rec)
        if i == 2:  # 0,1,2 => 3 secvente
            break

if len(records) < 3:
    raise ValueError(f"Am gasit doar {len(records)} secvente in {input_fastq_gz}, ai nevoie de cel putin 3.")

SeqIO.write(records, output_fasta, "fasta")

print(f"Am salvat {len(records)} secvente in {output_fasta}")