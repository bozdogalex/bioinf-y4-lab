#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Exercițiu (Lab 1): Descărcare FASTA + calcul GC

Scop:
  1) Descărcați un fișier FASTA de la NCBI (nucleotide sau proteină).
  2) Salvați fișierul local în data/work/<handle>/lab01/ (NU îl urcați pe git).
  3) Calculați fracția GC pentru fiecare înregistrare din fișier.

Instrucțiuni:
  - Rulați scriptul cu argumentele necesare (exemple):
      python ex01_multifasta_gc.py --email student@example.com \
        --query "TP53[Gene] AND Homo sapiens[Organism]" \
        --retmax 3 \
        --out data/work/<handle>/lab01/my_tp53.fa

      python ex01_multifasta_gc.py --email student@example.com \
        --accession NM_000546 \
        --out data/work/<handle>/lab01/nm000546.fa

  - Pași de completat:
    1) Configurați Entrez cu email (și api_key opțional).
    2) Dacă primiți accession → descărcați acel record cu efetch.
    3) Dacă primiți query → faceți esearch pentru IdList, apoi efetch pentru acele ID-uri.
    4) Scrieți rezultatele în fișierul dat prin --out.
    5) Citiți fișierul FASTA local și calculați GC pentru fiecare secvență.
    6) Afișați rezultatele pe ecran: <id>\tGC=<valoare cu 3 zecimale>.
"""

#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
from pathlib import Path
import sys

from Bio import SeqIO
from Bio import Entrez  # ✅ DEBLOCAT


def gc_fraction(seq: str) -> float:
    """Fracție GC pentru o secvență; robust la litere mici/mari și non-ATGC."""
    s = seq.upper()
    atgc = [c for c in s if c in ("A", "T", "G", "C")]
    if not atgc:
        return 0.0
    g = atgc.count("G")
    c = atgc.count("C")
    return (g + c) / float(len(atgc))


def download_fasta(email: str, out_path: Path, query: str = None,
                   accession: str = None, db: str = "nuccore",
                   retmax: int = 3, api_key: str = None) -> int:
    """
    ✅ IMPLEMENTAT: Descărcare din NCBI
    """
    # Configurează Entrez
    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key
    
    ids = []
    
    if accession:
        # Cazul 1: avem un accession specific
        handle = Entrez.esearch(db=db, term=accession, retmax=1)
        record = Entrez.read(handle)
        handle.close()
        ids = record["IdList"]
        
    elif query:
        # Cazul 2: căutare prin query
        handle = Entrez.esearch(db=db, term=query, retmax=retmax)
        record = Entrez.read(handle)
        handle.close()
        ids = record["IdList"]
    
    if not ids:
        print("⚠️  Nu s-au găsit înregistrări!")
        return 0
    
    # Descarcă secvențele în format FASTA
    handle = Entrez.efetch(db=db, id=ids, rettype="fasta", retmode="text")
    fasta_data = handle.read()
    handle.close()
    
    # Scrie în fișier
    with open(out_path, "w") as f:
        f.write(fasta_data)
    
    return len(ids)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--email", required=True, help="Email obligatoriu pentru NCBI Entrez")
    ap.add_argument("--api_key", help="NCBI API key (opțional)")
    ap.add_argument("--query", help="Ex: 'TP53[Gene] AND Homo sapiens[Organism]'")
    ap.add_argument("--accession", help="Ex: NM_000546")
    ap.add_argument("--db", default="nuccore", choices=["nuccore", "protein"])
    ap.add_argument("--retmax", type=int, default=3)
    ap.add_argument("--out", required=True, help="Fișier FASTA de ieșire")
    args = ap.parse_args()

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    # ✅ TODO COMPLETAT: Descarcă fișierul FASTA
    n = download_fasta(args.email, out_path, query=args.query,
                       accession=args.accession, db=args.db,
                       retmax=args.retmax, api_key=args.api_key)
    print(f"[ok] Am scris {n} înregistrări în: {out_path}")

    # ✅ TODO COMPLETAT: Citește și calculează GC
    records = SeqIO.parse(out_path, "fasta")
    
    for rec in records:
        gc_val = gc_fraction(str(rec.seq))
        print(f"{rec.id}\tGC={gc_val:.3f}")


if __name__ == "__main__":
    main()