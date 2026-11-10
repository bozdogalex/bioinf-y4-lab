#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Exercițiu (Lab 1): Descărcare FASTA + calcul GC
"""

import argparse
from pathlib import Path
import sys

from Bio import SeqIO
from Bio import Entrez


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

    if accession:
        ids = [accession]
    elif query:
        with Entrez.esearch(db=db, term=query, retmax=retmax) as handle:
            search_results = Entrez.read(handle)
        ids = search_results.get("IdList", [])
        if not ids:
            raise ValueError(f"Niciun rezultat găsit pentru query: {query!r}")
    else:
        raise ValueError("Trebuie furnizat fie --accession, fie --query.")

    with Entrez.efetch(db=db, id=ids, rettype="fasta", retmode="text") as handle:
        fasta_text = handle.read()

    out_path.write_text(fasta_text, encoding="utf-8")

    n = sum(1 for line in fasta_text.splitlines() if line.startswith(">"))
    return n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--email", required=True, help="Email obligatoriu pentru NCBI Entrez")
    ap.add_argument("--api_key", help="NCBI API key (opțional)")
    ap.add_argument("--query", help="Ex: 'TP53[Gene] AND Homo sapiens[Organism]'")
    ap.add_argument("--accession", help="Ex: NM_000546")
    ap.add_argument("--db", default="nuccore", choices=["nuccore"], help="GC is defined for nucleotides, not proteins")
    ap.add_argument("--retmax", type=int, default=3)
    ap.add_argument("--out", required=True, help="Fișier FASTA de ieșire")
    args = ap.parse_args()

    if not args.query and not args.accession:
        ap.error("Trebuie specificat fie --query, fie --accession.")

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    n = download_fasta(args.email, out_path, query=args.query,
                       accession=args.accession, db=args.db,
                       retmax=args.retmax, api_key=args.api_key)
    print(f"[ok] Am scris {n} înregistrări în: {out_path}")

    records = list(SeqIO.parse(str(out_path), "fasta"))

    if not records:
        print("[warn] Nicio secvență găsită în fișierul descărcat.", file=sys.stderr)
        return

    for rec in records:
        valoare = gc_fraction(str(rec.seq))
        print(f"{rec.id}\tGC={valoare:.3f}")


if __name__ == "__main__":
    main()
