import os
import gzip
from pathlib import Path
from Bio import SeqIO

# TODO: înlocuiți <handle> cu username-ul vostru GitHub
handle = "Ana-Maria-Bojan"

in_fastq_plain = Path(f"data/work/{handle}/lab03/your_reads.fastq")
in_fastq_gz = Path(f"data/work/{handle}/lab03/your_reads.fastq.gz")
out_report = Path(f"labs/03_formats&NGS/submissions/{handle}/qc_report_{handle}.txt")
out_report.parent.mkdir(parents=True, exist_ok=True)

# Selectați fișierul existent
if in_fastq_plain.exists():
    reader = SeqIO.parse(str(in_fastq_plain), "fastq")
elif in_fastq_gz.exists():
    # Biopython citește din file-like; folosim gzip.open(..., "rt")
    reader = SeqIO.parse(gzip.open(in_fastq_gz, "rt"), "fastq")
else:
    raise FileNotFoundError(
        f"Nu am găsit nici {in_fastq_plain} nici {in_fastq_gz}. "
        f"Rulați întâi ex03_fetch_fastq.py sau copiați un FASTQ propriu."
    )

num_reads = 0
total_length = 0
total_n = 0
total_phred = 0
total_bases = 0

# TODO: completați logica de agregare
for record in reader:
    num_reads += 1
    seq_str = str(record.seq)
    L = len(seq_str)
    total_length += L
    total_n += seq_str.upper().count("N")

    phred = record.letter_annotations.get("phred_quality", [])
    total_phred += sum(phred)
    total_bases += len(phred)
    

# TODO: calculați valorile finale (atenție la împărțiri la zero)
len_mean = (total_length / num_reads) if num_reads > 0 else 0.0
n_rate = (total_n / total_length) if total_length > 0 else 0.0
phred_mean = (total_phred / total_bases) if total_bases > 0 else 0.0

with open(out_report, "w", encoding="utf-8") as out:
    out.write(f"Reads: {num_reads}\n")
    out.write(f"Mean length: {len_mean:.2f}\n")
    out.write(f"N rate: {n_rate:.4f}\n")
    out.write(f"Mean Phred: {phred_mean:.2f}\n")

print(f"[OK] QC report -> {out_report.resolve()}")