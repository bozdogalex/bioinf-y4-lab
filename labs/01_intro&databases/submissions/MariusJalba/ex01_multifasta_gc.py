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

import argparse
from pathlib import Path
import sys

from Bio import SeqIO
from Bio import Entrez  # TODO: deblocați și folosiți pentru descărcare


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
    TODO: Implementați descărcarea din NCBI.
    Pași:
      - Configurați Entrez cu email (și api_key opțional).
      - Dacă avem accession: descărcați acel record.
      - Altfel, dacă avem query: faceți esearch -> lista de id-uri, apoi efetch.
      - Scrieți rezultatele în out_path.
      - Returnați numărul de înregistrări scrise.
    """
    Entrez.email = email
    if api_key:
      Entrez.api_key = api_key
    
    ids = []

    if accession:
      ids = [accession]
    elif query:
      print(f"Caut '{query}'")
      handle = Entrez.esearch(db=db, term=query, retmax=retmax)
      record = Entrez.read(handle)
      handle.close()
      ids = record["IdList"]
    
    if not ids:
      print("Nu s-au gasit rezultate")
      return 0
    
    print(f"Am gasit {len(ids)} inregistrari, descarc FASTA")
    handle = Entrez.efetch(db=db, id=ids, rettype="fasta", retmode="text")
    data = handle.read()
    handle.close()
    
    out_path.write_text(data)
    print(f"Au fost adaugate {len(ids)} inregistrari in {out_path}")
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

    # TODO: Apelați funcția download_fasta(...) și salvați rezultatele
    n = download_fasta(args.email, out_path, query=args.query,
                       accession=args.accession, db=args.db,
                       retmax=args.retmax, api_key=args.api_key)

    # TODO: Citiți fișierul FASTA cu SeqIO.parse
    records = list(SeqIO.parse(out_path,"fasta"))

    # TODO: Calculați GC pentru fiecare secvență și afișați rezultatele
    for rec in records:
        gc = gc_fraction(str(rec.seq))
        print(f"{rec.id}\tGC={gc:.3f}")


if __name__ == "__main__":
    main()
