from pathlib import Path
from Bio import Entrez, SeqIO

# NCBI cere un email valid
Entrez.email = "student@example.com"
# opțional: Entrez.api_key = "NCBI_API_KEY"

Path("data").mkdir(exist_ok=True)

search = Entrez.esearch(
    db="nucleotide",
    term="BRCA1[Gene] AND Homo sapiens[Organism]",
    retmax=1,
)
ids = Entrez.read(search)["IdList"]
print(f"Găsite {len(ids)} rezultate.")

if not ids:
    raise SystemExit("Niciun rezultat pentru BRCA1.")

acc = ids[0]
gb = Entrez.efetch(db="nucleotide", id=acc, rettype="gbwithparts", retmode="text")
out = Path("data/brca1.gb")
out.write_text(gb.read(), encoding="utf-8")

gb_record = SeqIO.read(out, "genbank")
gc = (gb_record.seq.count("G") + gb_record.seq.count("C")) / max(1, len(gb_record.seq))

print("Accession:", gb_record.id)
print("Length:", len(gb_record.seq))
print("GC fraction:", round(gc, 3))
print("First 50 nt:", gb_record.seq[:50])
