from pathlib import Path
from io import StringIO
from Bio import Entrez, SeqIO
import os

# setează email pentru NCBI (sau: export NCBI_EMAIL="emailul_tău")
Entrez.email = os.getenv("NCBI_EMAIL", "constantin.baragan@student.upt.ro")

DATA_DIR = Path("data")
DATA_DIR.mkdir(exist_ok=True)

QUERY = 'BRCA1[Gene] AND "Homo sapiens"[Organism]'
OUT_GB = DATA_DIR / "brca1.gb"

def gc_content(seq: str) -> float:
    s = seq.upper().replace("N", "")
    return 0.0 if not s else (s.count("G") + s.count("C")) / len(s)

# search
with Entrez.esearch(db="nucleotide", term=QUERY, retmax=1) as h:
    ids = Entrez.read(h)["IdList"]
print(f"Găsite {len(ids)} rezultate.")
if not ids:
    raise SystemExit("Niciun rezultat pentru BRCA1.")

acc = ids[0]

# GenBank
with Entrez.efetch(db="nucleotide", id=acc, rettype="gb", retmode="text") as h:
    OUT_GB.write_text(h.read(), encoding="utf-8")
gb_record = SeqIO.read(OUT_GB, "genbank")

# FASTA & GC
with Entrez.efetch(db="nucleotide", id=acc, rettype="fasta", retmode="text") as hf:
    fasta_rec = SeqIO.read(StringIO(hf.read()), "fasta")
seq = str(fasta_rec.seq)

gc = gc_content(seq)
print("ID:", acc)
print("Titlu:", gb_record.description)
print("Length:", len(seq), "bp")
print("GC fraction:", round(gc, 3))
print("First 50 nt:", seq[:50])
