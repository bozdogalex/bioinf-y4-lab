from pathlib import Path
from io import StringIO
from Bio import Entrez, SeqIO
import os

Entrez.email = "student@example.com"
Path("data").mkdir(exist_ok=True)

search = Entrez.esearch(
    db="nucleotide",
    term="BRCA1[Gene] AND Homo sapiens[Organism] AND mRNA",
    retmax=1,
)
ids = Entrez.read(search)["IdList"]
print(f"Găsite {len(ids)} rezultate.")

# search
with Entrez.esearch(db="nucleotide", term=QUERY, retmax=1) as h:
    ids = Entrez.read(h)["IdList"]
print(f"Găsite {len(ids)} rezultate.")
if not ids:
    raise SystemExit("Niciun rezultat pentru BRCA1.")

acc = ids[0]
print("Fetching GenBank record for:", acc)

with Entrez.efetch(db="nucleotide", id=acc, rettype="gbwithparts", retmode="text") as handle:
    gb_text = handle.read()

out = Path("data/work/AlexTGoCreative/brca1.gb")
out.write_text(gb_text, encoding="utf-8")

gb_record = SeqIO.read(out, "genbank")

if not gb_record.seq or len(gb_record.seq) == 0:
    raise ValueError("Recordul nu conține secvență ADN (probabil ai obținut un entry fără seq).")

gc = (gb_record.seq.count("G") + gb_record.seq.count("C")) / len(gb_record.seq)

gc = gc_content(seq)
print("ID:", acc)
print("Titlu:", gb_record.description)
print("Length:", len(seq), "bp")
print("GC fraction:", round(gc, 3))
print("First 50 nt:", seq[:50])
