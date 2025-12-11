import os
import gzip
from pathlib import Path
from Bio import SeqIO

# Username GitHub
handle = "IrisDanila"

in_fastq_plain = Path(f"data/work/{handle}/lab03/your_reads.fastq")
in_fastq_gz = Path(f"data/work/{handle}/lab03/your_reads.fastq.gz")
out_report = Path(f"labs/03_formats&NGS/submissions/{handle}/qc_report_{handle}.txt")
out_report.parent.mkdir(parents=True, exist_ok=True)

# Selectați fișierul existent
if in_fastq_plain.exists():
    print(f"Reading: {in_fastq_plain}")
    reader = SeqIO.parse(str(in_fastq_plain), "fastq")
elif in_fastq_gz.exists():
    print(f"Reading: {in_fastq_gz}")
    reader = SeqIO.parse(gzip.open(in_fastq_gz, "rt"), "fastq")
else:
    raise FileNotFoundError(
        f"Nu am găsit nici {in_fastq_plain} nici {in_fastq_gz}. "
        f"Rulați întâi ex01_fetch_fastq.py!"
    )

num_reads = 0
total_length = 0
total_n = 0
total_phred = 0
total_bases = 0

# Agregare statistici
print("Processing reads...")
for record in reader:
    num_reads += 1
    
    # Secvență
    seq_str = str(record.seq)
    seq_len = len(seq_str)
    
    total_length += seq_len
    total_bases += seq_len
    
    # Număr baze N
    n_count = seq_str.upper().count('N')
    total_n += n_count
    
    # Scoruri Phred
    phred = record.letter_annotations["phred_quality"]
    total_phred += sum(phred)
    
    # Progress
    if num_reads % 10000 == 0:
        print(f"Processed {num_reads:,} reads...", end='\r')

print(f"\nTotal reads processed: {num_reads:,}")

# Calculare valori finale
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

# Scriere raport
with open(out_report, "w", encoding="utf-8") as out:
    out.write("=" * 70 + "\n")
    out.write("FASTQ Quality Control Report\n")
    out.write("=" * 70 + "\n")
    out.write(f"Reads: {num_reads}\n")
    out.write(f"Mean length: {len_mean:.2f}\n")
    out.write(f"N rate: {n_rate:.4f}\n")
    out.write(f"Mean Phred: {phred_mean:.2f}\n")
    out.write("=" * 70 + "\n")
    out.write(f"\nQuality Assessment:\n")
    if phred_mean >= 30:
        out.write("  ✓ PASS - High quality (Phred >= 30)\n")
    elif phred_mean >= 20:
        out.write("  ⚠ WARNING - Moderate quality (20 <= Phred < 30)\n")
    else:
        out.write("  ✗ FAIL - Low quality (Phred < 20)\n")
    
    if n_rate < 0.01:
        out.write("  ✓ PASS - Low N-rate (< 1%)\n")
    else:
        out.write("  ⚠ WARNING - High N-rate (>= 1%)\n")

print(f"[OK] QC report -> {out_report.resolve()}")

# Print summary to console
print("\n" + "=" * 70)
print("SUMMARY:")
print("=" * 70)
print(f"Reads: {num_reads:,}")
print(f"Mean length: {len_mean:.2f} bp")
print(f"N rate: {n_rate:.4f} ({n_rate*100:.2f}%)")
print(f"Mean Phred: {phred_mean:.2f}")
print("=" * 70)