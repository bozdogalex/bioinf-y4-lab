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
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')  # Use non-GUI backend

# TODO: înlocuiți <handle> cu username-ul vostru GitHub
handle = "AlexTGoCreative"
your_reads = "SRR684066"

in_fastq_gz = Path(f"../../../../data/work/{handle}/lab03/{your_reads}.fastq.gz")
out_report = Path(f"qc_report_{handle}.txt")
out_plot = Path(f"qc_plot_{handle}.png")
out_report.parent.mkdir(parents=True, exist_ok=True)

# Selectați fișierul existent (doar .gz)
if in_fastq_gz.exists():
    # Biopython citește din file-like; folosim gzip.open(..., "rt")
    reader = SeqIO.parse(gzip.open(in_fastq_gz, "rt"), "fastq")
else:
    raise FileNotFoundError(
        f"Nu am găsit {in_fastq_gz}. "
        f"Rulați întâi ex01_fetch_fastq.py pentru a descărca fișierul FASTQ."
    )

num_reads = 0
total_length = 0
total_n = 0
total_phred = 0
total_bases = 0

# Lists for plotting distributions
read_lengths = []
mean_phred_scores = []

# TODO: completați logica de agregare
for record in reader:
    # HINT:
    # seq_str = str(record.seq)
    # phred = record.letter_annotations["phred_quality"]
    
    num_reads += 1
    seq_str = str(record.seq)
    seq_length = len(seq_str)
    total_length += seq_length
    total_bases += seq_length
    
    # Contorizăm bazele 'N'
    total_n += seq_str.count('N') + seq_str.count('n')
    
    # Calculăm suma scorurilor Phred
    if "phred_quality" in record.letter_annotations:
        phred = record.letter_annotations["phred_quality"]
        total_phred += sum(phred)
        
        # Store for plotting
        read_lengths.append(seq_length)
        mean_phred_scores.append(sum(phred) / len(phred) if len(phred) > 0 else 0)

# TODO: calculați valorile finale (atenție la împărțiri la zero)
len_mean = total_length / num_reads if num_reads > 0 else 0.0
n_rate = total_n / total_bases if total_bases > 0 else 0.0
phred_mean = total_phred / total_bases if total_bases > 0 else 0.0

with open(out_report, "w", encoding="utf-8") as out:
    out.write(f"Reads: {num_reads}\n")
    out.write(f"Mean length: {len_mean:.2f}\n")
    out.write(f"N rate: {n_rate:.4f}\n")
    out.write(f"Mean Phred: {phred_mean:.2f}\n")

print(f"[OK] QC report -> {out_report.resolve()}")

# Generate visualization (Bonus)
if read_lengths and mean_phred_scores:
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    # Plot 1: Read length distribution
    ax1.hist(read_lengths, bins=50, color='steelblue', edgecolor='black', alpha=0.7)
    ax1.set_xlabel('Read Length (bp)', fontsize=12)
    ax1.set_ylabel('Frequency', fontsize=12)
    ax1.set_title(f'Read Length Distribution\n(n={num_reads} reads)', fontsize=13, fontweight='bold')
    ax1.grid(axis='y', alpha=0.3)
    ax1.axvline(len_mean, color='red', linestyle='--', linewidth=2, label=f'Mean: {len_mean:.1f} bp')
    ax1.legend()
    
    # Plot 2: Mean Phred score distribution
    ax2.hist(mean_phred_scores, bins=50, color='forestgreen', edgecolor='black', alpha=0.7)
    ax2.set_xlabel('Mean Phred Score', fontsize=12)
    ax2.set_ylabel('Frequency', fontsize=12)
    ax2.set_title(f'Phred Quality Score Distribution\n(Mean: {phred_mean:.2f})', fontsize=13, fontweight='bold')
    ax2.grid(axis='y', alpha=0.3)
    ax2.axvline(phred_mean, color='red', linestyle='--', linewidth=2, label=f'Mean: {phred_mean:.2f}')
    ax2.legend()
    
    # Add quality threshold lines (Q20 and Q30)
    ax2.axvline(20, color='orange', linestyle=':', linewidth=1.5, alpha=0.7, label='Q20 (99% accuracy)')
    ax2.axvline(30, color='darkgreen', linestyle=':', linewidth=1.5, alpha=0.7, label='Q30 (99.9% accuracy)')
    ax2.legend()
    
    plt.tight_layout()
    plt.savefig(out_plot, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"[OK] QC plot -> {out_plot.resolve()}")
