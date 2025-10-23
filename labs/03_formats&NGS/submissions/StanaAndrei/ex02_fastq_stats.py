#!/usr/bin/env python3
"""
Exercițiu 04 — FASTQ QC pe date proprii

- Citește your_reads.fastq( .gz ) din data/work/<handle>/lab03/
- Calculează:
    * număr total de citiri
    * lungimea medie a citirilor
    * proporția bazelor 'N' (pe total baze)
    * scorul Phred mediu (medie per bază)
- Salvează raportul în labs/03_formats&NGS/submissions/<handle>/qc_report_<handle>.txt

Dacă nu modificați `handle`, scriptul va încerca să determine numele din
variabila de mediu HANDLE sau din USER/USERNAME; altfel folosește "student".
"""
import gzip
from pathlib import Path
from Bio import SeqIO

# TODO: înlocuiți <handle> cu username-ul vostru GitHub dacă doriți fix
handle = "StanaAndrei"


in_fastq_plain = Path(f"data/work/{handle}/sample/your_reads.fastq")
in_fastq_gz = Path(f"data/work/{handle}/sample/your_reads.fastq.gz")
out_report = Path(f"labs/03_formats&NGS/submissions/{handle}/qc_report_{handle}.txt")
out_report.parent.mkdir(parents=True, exist_ok=True)

# Alege fișierul existent și folosește context manager pentru a închide corect fișierul gzip
if in_fastq_plain.exists():
    source_path = in_fastq_plain
    opener = open
    opener_mode = "r"
elif in_fastq_gz.exists():
    source_path = in_fastq_gz
    opener = gzip.open
    opener_mode = "rt"
else:
    raise FileNotFoundError(
        f"Nu am găsit nici {in_fastq_plain} nici {in_fastq_gz}.\n"
        f"Rulați întâi ex03_fetch_fastq.py sau copiați un FASTQ propriu în data/work/{handle}/lab03/."
    )

num_reads = 0
total_length = 0          # suma lungimilor citirilor
total_n = 0               # număr total de baze N
total_phred = 0           # sumă a tuturor scorurilor phred (pe bază)
total_bases = 0           # număr total de baze (A/C/G/T/N/... pentru toate citirile)

# Parsăm fișierul FASTQ folosind Biopython
with opener(str(source_path), opener_mode) as fh:
    reader = SeqIO.parse(fh, "fastq")
    for record in reader:
        seq_str = str(record.seq)
        phred = record.letter_annotations.get("phred_quality", [])

        L = len(seq_str)
        if L == 0:
            # skip citiri 0-length (teoretic nu ar trebui să existe în fastq corect)
            continue

        num_reads += 1
        total_length += L
        total_n += seq_str.upper().count("N")
        total_bases += L

        # dacă lista phred are aceeași lungime ca secvența, sumăm; altfel ignorăm
        if phred and len(phred) == L:
            total_phred += sum(phred)
        else:
            # Dacă nu sunt scoruri sau lungimi nepotrivite, putem încerca să
            # estimăm (dar aici vom considera 0 pentru acele poziții).
            # Practic, dacă phred lipsesc, nu adăugăm nimic la total_phred.
            pass

# Calcul valori finale (atenție la împărțiri la zero)
if num_reads > 0:
    len_mean = total_length / num_reads
else:
    len_mean = 0.0

if total_bases > 0:
    n_rate = total_n / total_bases
    phred_mean = total_phred / total_bases
else:
    n_rate = 0.0
    phred_mean = 0.0

with open(out_report, "w", encoding="utf-8") as out:
    out.write(f"Handle: {handle}\n")
    out.write(f"Reads: {num_reads}\n")
    out.write(f"Mean length: {len_mean:.2f}\n")
    out.write(f"N rate: {n_rate:.6f}\n")
    out.write(f"Mean Phred (per base): {phred_mean:.2f}\n")
    out.write("\n")
    out.write("Notes:\n")
    out.write("- Mean Phred computed as total_phred / total_bases (per-bază).\n")
    out.write("- If Phred scores are missing for some reads, ele nu sunt incluse în total_phred.\n")

print(f"[OK] QC report -> {out_report.resolve()}")
