#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Exercițiu 04 — FASTQ QC pe date proprii

Scop:
- Citiți fișierul vostru FASTQ din data/work/<handle>/lab03/
- Calculați:
    * număr total de citiri
    * lungimea medie a citirilor
    * proporția bazelor 'N'
    * scorul Phred mediu
- Salvați raportul în:
    labs/03_formats&NGS/submissions/<handle>/qc_report_<handle>.txt
"""

import gzip
from pathlib import Path
from Bio import SeqIO

handle = "MarioCld" 

import glob

candidates = glob.glob(f"data/work/{handle}/lab03/*.fastq*")
if not candidates:
    raise FileNotFoundError(f"Nu am găsit niciun FASTQ în data/work/{handle}/lab03/")
in_fastq = Path(candidates[0])

if in_fastq.suffix == ".gz":
    reader = SeqIO.parse(gzip.open(in_fastq, "rt"), "fastq")
else:
    reader = SeqIO.parse(str(in_fastq), "fastq")

out_report = Path(f"labs/03_formats&NGS/submissions/{handle}/qc_report_{handle}.txt")
out_report.parent.mkdir(parents=True, exist_ok=True)

num_reads = 0
total_length = 0
total_n = 0
total_phred = 0
total_bases = 0

for record in reader:
    num_reads += 1
    seq_str = str(record.seq)
    phred = record.letter_annotations["phred_quality"]

    total_length += len(seq_str)
    total_n += seq_str.count("N")
    total_phred += sum(phred)
    total_bases += len(phred)

len_mean = total_length / num_reads if num_reads > 0 else 0.0
n_rate = total_n / total_length if total_length > 0 else 0.0
phred_mean = total_phred / total_bases if total_bases > 0 else 0.0

with open(out_report, "w", encoding="utf-8") as out:
    out.write(f"Reads: {num_reads}\n")
    out.write(f"Mean length: {len_mean:.2f}\n")
    out.write(f"N rate: {n_rate:.4f}\n")
    out.write(f"Mean Phred: {phred_mean:.2f}\n")

print(f"[OK] QC report -> {out_report.resolve()}")
