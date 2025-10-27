from pathlib import Path
from Bio import Entrez, SeqIO

Entrez.email = "student@example.com"

# Define query
QUERY = "BRCA1[Gene] AND Homo sapiens[Organism] AND mRNA"

# Create output folder
out_dir = Path("data/work/AlexTGoCreative")
out_dir.mkdir(parents=True, exist_ok=True)

# Search in NCBI nucleotide DB
with Entrez.esearch(db="nucleotide", term=QUERY, retmax=1) as search:
    ids = Entrez.read(search)["IdList"]

print(f"Găsite {len(ids)} rezultate.")
if not ids:
    raise SystemExit("Niciun rezultat pentru BRCA1.")

acc = ids[0]
print("Fetching GenBank record for:", acc)

# Fetch GenBank record
with Entrez.efetch(db="nucleotide", id=acc, rettype="gbwithparts", retmode="text") as handle:
    gb_text = handle.read()

# Save file
out_file = out_dir / "brca1.gb"
out_file.write_text(gb_text, encoding="utf-8")

# Parse GenBank
gb_record = SeqIO.read(out_file, "genbank")

seq = gb_record.seq
gc = (seq.count("G") + seq.count("C")) / len(seq)

# Output summary
print("ID:", acc)
print("Titlu:", gb_record.description)
print("Length:", len(seq), "bp")
print("GC fraction:", round(gc, 3))
print("First 50 nt:", seq[:50])
