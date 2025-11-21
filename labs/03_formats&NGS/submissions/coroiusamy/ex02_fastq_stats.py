"""
Exercițiu 04 — FASTQ QC pe date proprii

TODO:
- Citiți fișierul vostru FASTQ din data/work/<handle>/lab03/:
    your_reads.fastq  sau  your_reads.fastq.gz
- Calculați statistici:
    * număr total de citiri
    * lungimea medie a citirilor
    * proporția bazelor 'N'
    * scorul Phred mediu
- Salvați raportul în:
    labs/03_formats&NGS/submissions/<handle>/qc_report_<handle>.txt
"""

import os
import gzip
from pathlib import Path
from Bio import SeqIO

# TODO 1: înlocuiți <handle> cu username-ul vostru GitHub
handle = "coroiusamy"

in_fastq_plain = Path(f"data/work/{handle}/lab03/your_reads.fastq")
in_fastq_gz = Path(f"data/work/{handle}/lab03/your_reads.fastq.gz")
out_report = Path(f"labs/03_formats&NGS/submissions/{handle}/qc_report_{handle}.txt")
out_report.parent.mkdir(parents=True, exist_ok=True)

# Selectați fișierul existent
if in_fastq_plain.exists():
    print(f"[Info] Se citește fișierul: {in_fastq_plain}")
    reader_handle = open(in_fastq_plain, "rt")
    reader = SeqIO.parse(reader_handle, "fastq")
elif in_fastq_gz.exists():
    print(f"[Info] Se citește fișierul: {in_fastq_gz}")
    # Biopython citește din file-like; folosim gzip.open(..., "rt")
    reader_handle = gzip.open(in_fastq_gz, "rt")
    reader = SeqIO.parse(reader_handle, "fastq")
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

print("[Info] Se calculează statisticile QC...")

# TODO 2: completați logica de agregare
for record in reader:
    num_reads += 1
    
    seq_str = str(record.seq)
    seq_len = len(seq_str)
    
    total_length += seq_len
    total_n += seq_str.count("N")
    
    phred_scores = record.letter_annotations["phred_quality"]
    total_phred += sum(phred_scores)
    total_bases += len(phred_scores)

# inchidere fisier
reader_handle.close()

# TODO 3: calculați valorile finale (atenție la împărțiri la zero)
len_mean = total_length / num_reads if num_reads > 0 else 0.0
n_rate = total_n / total_length if total_length > 0 else 0.0
phred_mean = total_phred / total_bases if total_bases > 0 else 0.0

with open(out_report, "w", encoding="utf-8") as out:
    out.write(f"Reads: {num_reads}\n")
    out.write(f"Mean length: {len_mean:.2f}\n")
    out.write(f"N rate: {n_rate:.4f}\n")
    out.write(f"Mean Phred: {phred_mean:.2f}\n")

print(f"[OK] QC report -> {out_report.resolve()}")
