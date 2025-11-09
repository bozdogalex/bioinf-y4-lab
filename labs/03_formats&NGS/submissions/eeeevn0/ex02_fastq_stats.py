"""
Exercițiu 04 — FASTQ QC pe date proprii

- Citește FASTQ din data/work/<handle>/lab03/ (your_reads.fastq sau your_reads.fastq.gz)
- Calculează: număr citiri, lungime medie, proporție 'N', Phred mediu
- Salvează raport în labs/03_formats&NGS/submissions/<handle>/qc_report_<handle>.txt
"""

import gzip
from pathlib import Path
from Bio import SeqIO


handle = "eeeevn0"

in_fastq_plain = Path(f"data/work/{handle}/lab03/your_reads.fastq")
in_fastq_gz = Path(f"data/work/{handle}/lab03/your_reads.fastq.gz")
out_report = Path(f"labs/03_formats&NGS/submissions/{handle}/qc_report_{handle}.txt")
out_report.parent.mkdir(parents=True, exist_ok=True)


def open_fastq(path_plain: Path, path_gz: Path):
    """
    Deschide fișier FASTQ (plain sau gzipped) în mod text.
    """
    if path_plain.exists():
        print(f"[INFO] Deschid fișier text: {path_plain}")
        return open(path_plain, "rt")

    if path_gz.exists():
        print(f"[INFO] Deschid fișier comprimat: {path_gz}")
        with open(path_gz, "rb") as bfh:
            head = bfh.read(2)
        if head == b"\x1f\x8b":
            return gzip.open(path_gz, "rt")  # <- cheia
        else:
            return open(path_gz, "rt")

    return None


fh = open_fastq(in_fastq_plain, in_fastq_gz)
if fh is None:
    raise FileNotFoundError(
        f"Nu am găsit nici {in_fastq_plain} nici {in_fastq_gz}. "
        "Rulați întâi ex03_fetch_fastq.py sau copiați un FASTQ propriu în acest director."
    )

reader = SeqIO.parse(fh, "fastq")

num_reads = 0
total_length = 0
total_n = 0
total_phred = 0
total_bases = 0

try:
    for record in reader:
        num_reads += 1
        seq_str = str(record.seq)
        L = len(seq_str)
        total_length += L
        total_n += seq_str.upper().count("N")
        phred = record.letter_annotations.get("phred_quality", [])
        total_phred += sum(phred)
        total_bases += len(phred)
finally:
    fh.close()

len_mean = total_length / num_reads if num_reads > 0 else 0.0
n_rate = total_n / total_length if total_length > 0 else 0.0
phred_mean = total_phred / total_bases if total_bases > 0 else 0.0

with open(out_report, "w") as out:
    out.write(f"Reads: {num_reads}\n")
    out.write(f"Mean length: {len_mean:.2f}\n")
    out.write(f"N rate: {n_rate:.4f}\n")
    out.write(f"Mean Phred: {phred_mean:.2f}\n")

print(f"[OK] QC report -> {out_report.resolve()}")
