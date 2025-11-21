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
handle = "indibotoc"
in_fastq_plain = Path(f"data/work/{handle}/lab03/your_reads.fastq")
in_fastq_gz    = Path(f"data/work/{handle}/lab03/your_reads.fastq.gz")
out_report     = Path(f"labs/03_formats&NGS/submissions/{handle}/qc_report_{handle}.txt")
out_report.parent.mkdir(parents=True, exist_ok=True)

def parse_fastq_stream(path: Path):
    if path.suffix == ".gz":
        with gzip.open(path, "rt") as fh:
            for rec in SeqIO.parse(fh, "fastq"):
                yield rec
    else:
        with open(path, "rt") as fh:
            for rec in SeqIO.parse(fh, "fastq"):
                yield rec

in_path = in_fastq_plain if in_fastq_plain.exists() else in_fastq_gz
if not in_path or not in_path.exists():
    raise FileNotFoundError("File not found!")

num_reads   = 0
total_len   = 0
total_n     = 0
total_phred = 0

# TODO: completați logica de agregare
for rec in parse_fastq_stream(in_path):
    seq = str(rec.seq)
    qs  = rec.letter_annotations.get("phred_quality", [])
    L   = len(seq)
    if L == 0:
        continue
    num_reads   += 1
    total_len   += L
    total_n     += seq.upper().count("N")
    total_phred += sum(qs)

len_mean   = (total_len / num_reads) if num_reads else 0.0
n_rate     = (total_n   / total_len) if total_len else 0.0
phred_mean = (total_phred / total_len) if total_len else 0.0

with open(out_report, "w", encoding="utf-8") as out:
    out.write(f"FASTQ QC report for: {in_path}\n")
    out.write(f"Reads: {num_reads}\n")
    out.write(f"Mean length: {len_mean:.2f}\n")
    out.write(f"N rate: {n_rate:.6f}\n")
    out.write(f"Mean Phred: {phred_mean:.2f}\n")

print(f"[OK] QC report -> {out_report.resolve()}")
