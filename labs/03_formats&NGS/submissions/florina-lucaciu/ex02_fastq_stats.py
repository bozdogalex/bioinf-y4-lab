#!/usr/bin/env python3

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

#in terminal: python labs/'03_formats&NGS'/submissions/florina-lucaciu/ex02_fastq_stats.py

import os
import gzip
from pathlib import Path
from Bio import SeqIO

# TODO: înlocuiți <handle> cu username-ul vostru GitHub
handle = "florina-lucaciu" 

in_fastq_plain = Path(f"data/work/{handle}/lab03/your_reads.fastq")
in_fastq_gz = Path(f"data/work/{handle}/lab03/your_reads.fastq.gz")
out_report = Path(f"labs/03_formats&NGS/submissions/{handle}/qc_report_{handle}.txt")
out_report.parent.mkdir(parents=True, exist_ok=True)

# selectati fișierul existent
if in_fastq_plain.exists():
    print(f"Se citește fișierul: {in_fastq_plain}")
    reader = SeqIO.parse(str(in_fastq_plain), "fastq")
elif in_fastq_gz.exists():
    print(f"Se citește fișierul: {in_fastq_gz}")
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

print("Se calculează statisticile...")

# TODO: completați logica de agregare
for record in reader:
    # umaram citirea
    num_reads += 1
    
    # obtinem secventa ca string si lungimea ei
    seq_str = str(record.seq)
    current_length = len(seq_str)
    
    # adaugam la lungimea totala (numarul total de baze)
    total_length += current_length
    
    # um bazele 'N'
    total_n += seq_str.upper().count('N')
    
    # obt scorurile Phred si le insumam
    phred_scores = record.letter_annotations["phred_quality"]
    total_phred += sum(phred_scores)

# TODO: calculați valorile finale (atenție la împărțiri la zero)
len_mean = 0.0
n_rate = 0.0
phred_mean = 0.0


if num_reads > 0:
    len_mean = total_length / num_reads

if total_length > 0:
    n_rate = total_n / total_length
    phred_mean = total_phred / total_length

with open(out_report, "w", encoding="utf-8") as out:
    out.write(f"--- Raport QC pentru {handle} ---\n")
    out.write(f"Reads: {num_reads}\n")
    out.write(f"Mean length: {len_mean:.2f}\n")
    out.write(f"N rate: {n_rate:.4f}\n")
    out.write(f"Mean Phred: {phred_mean:.2f}\n")

print(f"[OK] QC report -> {out_report.resolve()}")