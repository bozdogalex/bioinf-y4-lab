
import gzip
from pathlib import Path
from Bio import SeqIO

handle = "MariusC12"

in_fastq = Path(f"data/work/{handle}/lab03/your_reads.fastq.gz")
out_report = Path(
    f"labs/03_formats&NGS/submissions/{handle}/qc_report_{handle}.txt"
)

num_reads = 0
total_length = 0
total_n = 0
total_phred = 0
total_bases = 0

with gzip.open(in_fastq, "rt") as fastq:
    for record in SeqIO.parse(fastq, "fastq"):
        sequence = str(record.seq)
        phred = record.letter_annotations["phred_quality"]

        num_reads += 1
        total_length += len(sequence)
        total_n += sequence.upper().count("N")
        total_phred += sum(phred)
        total_bases += len(sequence)

if num_reads == 0 or total_bases == 0:
    raise ValueError("FASTQ-ul nu contine citiri.")

len_mean = total_length / num_reads
n_rate = total_n / total_bases
phred_mean = total_phred / total_bases

with open(out_report, "w", encoding="utf-8") as out:
    out.write(f"Reads: {num_reads}\n")
    out.write(f"Mean length: {len_mean:.2f}\n")
    out.write(f"N rate: {n_rate:.4f}\n")
    out.write(f"Mean Phred: {phred_mean:.2f}\n")

print(f"[OK] QC report -> {out_report.resolve()}")
