#!/usr/bin/env python3
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
from io import StringIO

from Bio import SeqIO
from Bio import Entrez # TODO: deblocați și folosiți pentru descărcare


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
    Descarcă FASTA din NCBI.
      - dacă avem accession: efetch direct
      - altfel, dacă avem query: esearch -> efetch pentru IdList
    Scrie rezultatele în out_path și întoarce numărul de înregistrări.
    """
    # Configurăm Entrez
    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key

    out_path.parent.mkdir(parents=True, exist_ok=True)

    # Descărcare prin accession direct
    if accession:
        with Entrez.efetch(db=db, id=accession, rettype="fasta", retmode="text") as handle:
            data = handle.read()
        if not data.strip():
            print("N-am primit conținut FASTA pentru accession.", file=sys.stderr)
            return 0
        out_path.write_text(data, encoding="utf-8")
        return sum(1 for _ in SeqIO.parse(StringIO(data), "fasta"))

    # Descărcare prin query
    if query:
        with Entrez.esearch(db=db, term=query, retmax=retmax) as search:
            ids = Entrez.read(search)["IdList"]
        if not ids:
            print("Niciun rezultat găsit pentru query.", file=sys.stderr)
            return 0
        with Entrez.efetch(db=db, id=",".join(ids), rettype="fasta", retmode="text") as fetch:
            data = fetch.read()
        if not data.strip():
            print("eFetch nu a întors FASTA pentru ID-urile găsite.", file=sys.stderr)
            return 0
        out_path.write_text(data, encoding="utf-8")
        return sum(1 for _ in SeqIO.parse(StringIO(data), "fasta"))

    print("Trebuie specificat --query sau --accession.", file=sys.stderr)
    return 0


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

    
    if (args.query is None) == (args.accession is None):
        ap.error("Folosește exact UNUL dintre --query sau --accession.")

    out_path = Path(args.out)
    n = download_fasta(args.email, out_path,
                       query=args.query, accession=args.accession,
                       db=args.db, retmax=args.retmax, api_key=args.api_key)
    if n == 0:
        print("Niciun record scris. Verifică parametrii.", file=sys.stderr)
        sys.exit(1)

    print(f"[ok] Am scris {n} înregistrări în: {out_path}")

    #citim și afișăm GC pentru fiecare secvență
    records = list(SeqIO.parse(out_path, "fasta"))
    vals = []
    for rec in records:
        gc = gc_fraction(str(rec.seq))
        vals.append(gc)
        print(f"{rec.id}\tGC={gc:.3f}\tLength={len(rec.seq)}")

    if vals:
        print(f"Average GC={sum(vals)/len(vals):.3f}")


if __name__ == "__main__":
    main()
