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

# TODO: înlocuiți <handle> cu username-ul vostru GitHub
handle = "rauldavid35"
your_reads="ERR071289"

in_fastq_gz = Path(f"data/work/{handle}/lab03/{your_reads}.fastq.gz")
out_report = Path(f"labs/03_formats&NGS/submissions/{handle}/qc_report_{handle}.txt")
out_report.parent.mkdir(parents=True, exist_ok=True)

# Selectați fișierul existent
if in_fastq_gz.exists():
    # Biopython citește din file-like; folosim gzip.open(..., "rt")
    reader = SeqIO.parse(gzip.open(in_fastq_gz, "rt"), "fastq")
else:
    raise FileNotFoundError(
        f"Nu am găsit {in_fastq_gz}. "
        f"Rulați întâi ex01_fetch_fastq.py sau copiați un FASTQ propriu."
    )

num_reads = 0
total_length = 0
total_n = 0
total_phred = 0
total_bases = 0

# TODO: completați logica de agregare
for record in reader:
    # HINT:
    # seq_str = str(record.seq)
    # phred = record.letter_annotations["phred_quality"]
    num_reads+=1
    s=str(record.seq)
    L=len(s)
    total_length+=L
    total_bases+=L
    s_upp=s.upper()
    total_n+=s_upp.count("N")

    if "phred_quality" in record.letter_annotations:
        phred=record.letter_annotations["phred_quality"]
        total_phred+=sum(phred)

# TODO: calculați valorile finale (atenție la împărțiri la zero)
if num_reads:
    len_mean = total_length / num_reads
else:
    len_mean = 0.0

if total_bases:
    n_rate = total_n / total_bases
    phred_mean = total_phred / total_bases
else:
    n_rate = 0.0
    phred_mean = 0.0

with open(out_report, "w", encoding="utf-8") as out:
    out.write(f"Reads: {num_reads}\n")
    out.write(f"Mean length: {len_mean:.2f}\n")
    out.write(f"N rate: {n_rate:.4f}\n")
    out.write(f"Mean Phred: {phred_mean:.2f}\n")

print(f"[OK] QC report -> {out_report.resolve()}")
