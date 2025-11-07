# generate_my_fasta.py

from pathlib import Path
from Bio import Entrez, SeqIO


Entrez.email = "codreanumihaic@gmail.com" 
genes = {
    "NC_060941.1": "NC_060941.1",
    "NC_000017.11": "NC_000017.11",  
    "NG_017013.2": "NG_017013.2"
}

out_path = Path("data/work/codmihaic/lab04/your_sequences.fasta")
out_path.parent.mkdir(parents=True, exist_ok=True)
records = []

print("[INFO] Descărcare secvențe de la NCBI...")
for name, acc in genes.items():
    print(f"  - {name} ({acc})")
    with Entrez.efetch(db="nucleotide", id=acc, rettype="fasta", retmode="text") as handle:
        record = SeqIO.read(handle, "fasta")
        record.id = name
        record.description = f"{name} ({acc})"
        records.append(record)

SeqIO.write(records, out_path, "fasta")
print(f"[OK] Fișier salvat: {out_path.resolve()}")
