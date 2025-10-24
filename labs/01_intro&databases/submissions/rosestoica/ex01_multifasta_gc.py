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
      python /workspaces/bioinf-y4-lab/labs/"01_intro&databases"/submissions/rosestoica/ex01_multifasta_gc.py --email rosemariestoica13@gmail.com \
        --query "TP53[Gene] AND Homo sapiens[Organism]" \
        --retmax 3 \
        --out data/my_tp53.fa

      python /workspaces/bioinf-y4-lab/labs/"01_intro&databases"/submissions/rosestoica/ex01_multifasta_gc.py --email rosemariestoica13@gmail.com \
        --accession NM_000546 \
        --out data/nm000546.fa

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
from Bio import Entrez  
from io import StringIO



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
    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key
    ids = []
    
    if accession:
        ids = [accession]
    elif query:
        with Entrez.esearch(db=db, term=query, retmax=retmax) as handle:
            search_results = Entrez.read(handle)
            ids = search_results.get("IdList", [])
        if not ids:
            print(" Nu s-au găsit înregistrări pentru query-ul dat.")
            return 0
    else:
        raise ValueError("Trebuie să specificați fie accession, fie query.")
    records = []
    for record_id in ids:
        with Entrez.efetch(db=db, id=record_id, rettype="fasta", retmode="text") as fetch_handle:
            fasta_text = fetch_handle.read()
            fasta_rec = SeqIO.read(StringIO(fasta_text), "fasta")
            records.append(fasta_rec)
    SeqIO.write(records, out_path, "fasta")
    
    return len(records)


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

    n = download_fasta(email=args.email, out_path=out_path, query=args.query,
                       accession=args.accession, db=args.db,
                       retmax=args.retmax, api_key=args.api_key)
    print(f"[ok] Am scris {n} înregistrări în: {out_path}")

    records = list(SeqIO.parse(out_path, "fasta"))

    for rec in records:
       gc = gc_fraction(str(rec.seq))
       print(f"{rec.id}\tGC={gc:.3f}")


if __name__ == "__main__":
    main()
