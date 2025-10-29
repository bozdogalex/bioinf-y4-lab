"""
Exercitiu 04 - FASTQ QC pe date proprii

Obiectiv:
- Cititi fisierul vostru FASTQ din data/work/<handle>/lab03/:
    your_reads.fastq  sau  your_reads.fastq.gz
- Calculati statistici:
    * numar total de citiri
    * lungimea medie a citirilor
    * proportia bazelor 'N'
    * scorul Phred mediu
- Salvati raportul in:
    labs/03_formats&NGS/submissions/<handle>/qc_report_<handle>.txt
"""

import gzip
from contextlib import contextmanager
from pathlib import Path

from Bio import SeqIO

# Inlocuiti <handle> cu username-ul vostru GitHub
handle = "Papaia462819"

in_fastq_plain = Path(f"data/work/{handle}/lab03/your_reads.fastq")
in_fastq_gz = Path(f"data/work/{handle}/lab03/your_reads.fastq.gz")
out_report = Path(f"labs/03_formats&NGS/submissions/{handle}/qc_report_{handle}.txt")
out_report.parent.mkdir(parents=True, exist_ok=True)


@contextmanager
def open_fastq():
    """
    Context manager care intoarce un iterator SeqIO.parse pentru FASTQ-ul gasit.
    Prioritizeaza varianta necomprimata, dar accepta si .gz.
    """
    if in_fastq_plain.exists():
        with in_fastq_plain.open("rt", encoding="utf-8") as fh:
            yield SeqIO.parse(fh, "fastq")
        return
    if in_fastq_gz.exists():
        with gzip.open(in_fastq_gz, "rt") as fh:
            yield SeqIO.parse(fh, "fastq")
        return
    raise FileNotFoundError(
        f"Nu am gasit nici {in_fastq_plain} nici {in_fastq_gz}. "
        "Rulati ex01_fetch_fastq.py sau copiati un FASTQ propriu in locatie."
    )


num_reads = 0
total_length = 0
total_n = 0
total_phred = 0
total_bases = 0

with open_fastq() as reader:
    for record in reader:
        seq_str = str(record.seq)
        phred_scores = record.letter_annotations.get("phred_quality", [])

        read_len = len(seq_str)
        num_reads += 1
        total_length += read_len
        total_n += seq_str.upper().count("N")

        total_phred += sum(phred_scores)
        total_bases += len(phred_scores)

# Calculam valorile finale (tinand cont de divizarile la zero).
len_mean = total_length / num_reads if num_reads else 0.0
n_rate = total_n / total_length if total_length else 0.0
phred_mean = total_phred / total_bases if total_bases else 0.0

with out_report.open("w", encoding="utf-8") as out:
    out.write(f"Reads: {num_reads}\n")
    out.write(f"Mean length: {len_mean:.2f}\n")
    out.write(f"N rate: {n_rate:.4f}\n")
    out.write(f"Mean Phred: {phred_mean:.2f}\n")

print(f"[OK] QC report -> {out_report.resolve()}")
