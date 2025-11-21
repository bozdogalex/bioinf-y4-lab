"""
Exercițiu 04 — FASTQ QC pe date proprii

Obiectiv:
- Citiți fișierul vostru FASTQ din data/work/<handle>/lab03/
- Calculați statistici de calitate
- Salvați raportul în labs/03_formats&NGS/submissions/<handle>/

Exemplu folosire:
    python ex04_fastq_qc.py
"""

import os
import gzip
from pathlib import Path
from Bio import SeqIO

handle = "olezhka0809" 

in_fastq_plain = Path(f"data/work/{handle}/lab03/your_file.fastq")
in_fastq_gz = Path(f"data/work/{handle}/lab03/SRR000001.fastq.gz")
out_report = Path(f"qc_report_{handle}.txt")
out_report.parent.mkdir(parents=True, exist_ok=True)

print(f"Căutare fișier FASTQ în: data/work/{handle}/lab03/")

# Selectați fișierul existent
reader = None
input_file = None

if in_fastq_plain.exists():
    print(f"Găsit: {in_fastq_plain}")
    reader = SeqIO.parse(str(in_fastq_plain), "fastq")
    input_file = str(in_fastq_plain)
elif in_fastq_gz.exists():
    print(f"Găsit: {in_fastq_gz}")
    # Biopython citește din file-like; folosim gzip.open(..., "rt")
    reader = SeqIO.parse(gzip.open(in_fastq_gz, "rt"), "fastq")
    input_file = str(in_fastq_gz)
else:
    # Verificăm ce fișiere există în director
    lab03_dir = Path(f"data/work/{handle}/lab03/")
    if lab03_dir.exists():
        files = list(lab03_dir.glob("*.fastq*"))
        if files:
            print(f"Fișiere găsite în {lab03_dir}:")
            for f in files:
                print(f"  - {f.name}")
            print(f"\nSugestie: Redenumește unul din fișiere în 'your_reads.fastq' sau 'your_reads.fastq.gz'")
        else:
            print(f"Directorul {lab03_dir} este gol.")
    else:
        print(f"Directorul {lab03_dir} nu există.")
    
    raise FileNotFoundError(
        f"Nu am găsit nici {in_fastq_plain} nici {in_fastq_gz}. "
        f"Rulați întâi ex03_download_fastq.py sau copiați un FASTQ propriu."
    )

print("Procesare FASTQ...")

# Variabile pentru statistici
num_reads = 0
total_length = 0
total_n = 0
total_phred = 0
total_bases = 0

# Statistici adiționale (opțional)
min_length = float('inf')
max_length = 0
min_phred = float('inf')
max_phred = 0

# Procesare FASTQ
for record in reader:
    num_reads += 1
    
    # Extrage secvența ca string
    seq_str = str(record.seq)
    seq_length = len(seq_str)
    
    # Actualizare statistici lungime
    total_length += seq_length
    min_length = min(min_length, seq_length)
    max_length = max(max_length, seq_length)
    
    # Numără bazele 'N'
    total_n += seq_str.upper().count("N")
    
    # Extrage scorurile Phred de calitate
    phred = record.letter_annotations["phred_quality"]
    
    # Actualizare statistici Phred
    phred_sum = sum(phred)
    total_phred += phred_sum
    total_bases += len(phred)
    
    # Min/max Phred
    if phred:
        min_phred = min(min_phred, min(phred))
        max_phred = max(max_phred, max(phred))
    
    # Progress indicator (la fiecare 10000 citiri)
    if num_reads % 10000 == 0:
        print(f"  Procesate: {num_reads} citiri...", end="\r")

print(f"  Procesate: {num_reads} citiri... Done!")

# Calculează valorile finale (atenție la împărțiri la zero)
len_mean = total_length / num_reads if num_reads > 0 else 0.0
n_rate = total_n / total_length if total_length > 0 else 0.0
phred_mean = total_phred / total_bases if total_bases > 0 else 0.0

# Handle edge cases pentru min/max
if min_length == float('inf'):
    min_length = 0
if min_phred == float('inf'):
    min_phred = 0

# Scrie raportul
print(f"\nScriere raport în: {out_report}")

with open(out_report, "w", encoding="utf-8") as out:
    out.write("=" * 50 + "\n")
    out.write("FASTQ Quality Control Report\n")
    out.write("=" * 50 + "\n\n")
    
    out.write(f"Input file: {input_file}\n")
    out.write(f"Author: {handle}\n\n")
    
    out.write("--- Basic Statistics ---\n")
    out.write(f"Reads: {num_reads}\n")
    out.write(f"Mean length: {len_mean:.2f}\n")
    out.write(f"Min length: {min_length}\n")
    out.write(f"Max length: {max_length}\n")
    out.write(f"Total bases: {total_length}\n\n")
    
    out.write("--- Quality Statistics ---\n")
    out.write(f"N rate: {n_rate:.4f}\n")
    out.write(f"N count: {total_n}\n")
    out.write(f"Mean Phred: {phred_mean:.2f}\n")
    out.write(f"Min Phred: {min_phred}\n")
    out.write(f"Max Phred: {max_phred}\n\n")
    
    # Interpretare calitate
    out.write("--- Quality Interpretation ---\n")
    if phred_mean >= 30:
        out.write("Quality: EXCELLENT (Phred ≥ 30, error rate < 0.1%)\n")
    elif phred_mean >= 20:
        out.write("Quality: GOOD (Phred ≥ 20, error rate < 1%)\n")
    elif phred_mean >= 10:
        out.write("Quality: POOR (Phred ≥ 10, error rate < 10%)\n")
    else:
        out.write("Quality: VERY POOR (Phred < 10, error rate > 10%)\n")
    
    if n_rate > 0.05:
        out.write(f"WARNING: High N rate ({n_rate:.2%}) - many ambiguous bases\n")
    
    out.write("\n" + "=" * 50 + "\n")

print(f"[OK] QC report -> {out_report.resolve()}")

# Afișează și în terminal
print("\n" + "=" * 50)
print("REZUMAT:")
print(f"  Reads: {num_reads}")
print(f"  Mean length: {len_mean:.2f}")
print(f"  N rate: {n_rate:.4f} ({n_rate*100:.2f}%)")
print(f"  Mean Phred: {phred_mean:.2f}")
print("=" * 50)