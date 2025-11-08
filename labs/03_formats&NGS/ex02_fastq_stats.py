import os
import gzip
from Bio import SeqIO

handle = "MadalinaNegru"
fastq_file = f"data/work/MadalinaNegru/lab03/your_reads.fastq.gz"
outdir = f"submissions/MadalinaNegru"
os.makedirs(outdir, exist_ok=True)
out_report = os.path.join(outdir, f"qc_report_MadalinaNegru.txt")

# Numărul maxim de secvențe pentru QC rapid
max_reads = 1000

num_reads = 0
total_len = 0
min_len = None
max_len = None

# Deschidere fișier gzipped
with gzip.open(fastq_file, "rt") as f:
    for record in SeqIO.parse(f, "fastq"):
        seq_len = len(record.seq)
        total_len += seq_len
        num_reads += 1
        if min_len is None or seq_len < min_len:
            min_len = seq_len
        if max_len is None or seq_len > max_len:
            max_len = seq_len
        if num_reads >= max_reads:
            break

avg_len = total_len / num_reads if num_reads else 0

# Scrie raport QC
with open(out_report, "w") as f:
    f.write(f"Număr secvențe analizate: {num_reads}\n")
    f.write(f"Lungime medie secvențe: {avg_len:.2f}\n")
    f.write(f"Lungime minimă secvență: {min_len}\n")
    f.write(f"Lungime maximă secvență: {max_len}\n")

print(f"QC complet! Raport salvat în {out_report}")
