from pathlib import Path
from Bio import SeqIO
import gzip


input_fastq_gz = Path("data/work/eeeevn0/lab03/your_reads.fastq.gz") 

output_fasta = Path("data/work/eeeevn0/lab04/my_sequences.fasta")
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
