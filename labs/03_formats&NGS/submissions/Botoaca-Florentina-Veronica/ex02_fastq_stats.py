"""
Exercițiu 04 — FASTQ QC pe date proprii
"""

import os
import gzip
from pathlib import Path
from Bio import SeqIO

handle = "Botoaca-Florentina-Veronica"

in_fastq_plain = Path(f"data/work/{handle}/lab03/your_reads.fastq")
in_fastq_gz = Path(f"data/work/{handle}/lab03/your_reads.fastq.gz")

alternative_path = Path("/workspaces/bioinf-y4-lab/labs/01_intro&databases/data/work/Botoaca-Florentina-Veronica/lab01/1_control_psbA3_2019_minq7.fastq")

out_report = Path(f"/workspaces/bioinf-y4-lab/labs/03_formats&NGS/submissions/Botoaca-Florentina-Veronica/qc_report_Botoaca-Florentina-Veronica.txt")
out_report.parent.mkdir(parents=True, exist_ok=True)

# VERIFICĂM toate căile posibile
possible_paths = [
    in_fastq_plain,
    in_fastq_gz,
    alternative_path,
    Path(f"data/work/{handle}/lab03/your_reads.fastq"),
    Path(f"data/work/{handle}/lab03/your_reads.fastq.gz") 
]

fastq_file = None
file_type = None

for path in possible_paths:
    if path.exists():
        fastq_file = path
        # Detectăm tipul fișierului verificând primele bytes
        with open(path, 'rb') as f:
            magic = f.read(2)
            if magic == b'\x1f\x8b':  # GZIP magic number
                file_type = 'gzip'
            else:
                file_type = 'plain'
        print(f"✅ Fișier găsit: {path} (tip: {file_type})")
        break

if fastq_file is None:
    print("❌ Nu s-a găsit niciun fișier FASTQ. Căi verificate:")
    for path in possible_paths:
        print(f"   - {path}")
    print("\n👉 Rulați mai întâi: python ex03_fetch_fastq.py SRR2602530")
    exit(1)

# Citim fișierul în funcție de tipul detectat
if file_type == 'gzip':
    reader = SeqIO.parse(gzip.open(fastq_file, "rt"), "fastq")
    print(f"📖 Citind fișier GZIP: {fastq_file}")
else:
    reader = SeqIO.parse(str(fastq_file), "fastq")
    print(f"📖 Citind fișier plain: {fastq_file}")

num_reads = 0
total_length = 0
total_n = 0
total_phred = 0
total_bases = 0

print("🔬 Procesez citirile FASTQ...")
for record in reader:
    num_reads += 1
    
    # Extrage secvența ca string
    seq_str = str(record.seq)
    total_length += len(seq_str)
    
    # Numără bazele 'N'
    total_n += seq_str.count("N")
    
    # Extrage scorurile Phred
    phred_scores = record.letter_annotations["phred_quality"]
    total_phred += sum(phred_scores)
    total_bases += len(phred_scores)
    
    # Afișează progres la fiecare 1000 de citiri
    if num_reads % 1000 == 0:
        print(f"📊 Procesat {num_reads} citiri...")

print(f"✅ Procesare completă: {num_reads} citiri totale")

# Calculați valorile finale (cu verificare împărțire la zero)
if num_reads > 0:
    len_mean = total_length / num_reads
else:
    len_mean = 0.0

if total_length > 0:
    n_rate = total_n / total_length
else:
    n_rate = 0.0

if total_bases > 0:
    phred_mean = total_phred / total_bases
else:
    phred_mean = 0.0

# Scrie raportul
with open(out_report, "w", encoding="utf-8") as out:
    out.write("=== RAPORT CALITATE FASTQ ===\n")
    out.write(f"Fișier analizat: {fastq_file}\n")
    out.write(f"Tip fișier: {file_type}\n")
    out.write("=" * 40 + "\n")
    out.write(f"Reads: {num_reads}\n")
    out.write(f"Mean length: {len_mean:.2f}\n")
    out.write(f"N rate: {n_rate:.4f}\n")
    out.write(f"Mean Phred: {phred_mean:.2f}\n")
    out.write("=" * 40 + "\n")
    out.write(f"Total bases: {total_length}\n")
    out.write(f"Total N bases: {total_n}\n")
    out.write(f"N percentage: {n_rate * 100:.2f}%\n")

# Afișează rezultatele
print("\n📊 REZULTATE CALITATE FASTQ:")
print(f"📖 Reads: {num_reads}")
print(f"📏 Mean length: {len_mean:.2f}")
print(f"🔍 N rate: {n_rate:.4f} ({n_rate * 100:.2f}%)")
print(f"✅ Mean Phred: {phred_mean:.2f}")
print(f"💾 Raport salvat: {out_report.resolve()}")