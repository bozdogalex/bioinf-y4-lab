#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Exercițiu (Lab 1): Descărcare FASTA + calcul GC
"""

import argparse
from pathlib import Path
from Bio import Entrez, SeqIO


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
    Descărcare din NCBI cu Entrez.
    """
    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key

    ids = []
    if accession:
        ids = [accession]
    elif query:
        print(f"[i] Caut: {query}")
        with Entrez.esearch(db=db, term=query, retmax=retmax) as handle:
            result = Entrez.read(handle)
            ids = result["IdList"]
        print(f"[i] Găsit {len(ids)} înregistrări.")
    else:
        raise ValueError("Trebuie specificat fie --accession, fie --query")

    if not ids:
        raise ValueError("Nu s-au găsit înregistrări pentru criteriile date.")

    print(f"[i] Descarc {len(ids)} înregistrări din {db}...")
    with Entrez.efetch(db=db, id=ids, rettype="fasta", retmode="text") as handle:
        data = handle.read()

    out_path.write_text(data, encoding="utf-8")

    # Număr înregistrări FASTA scrise
    num_records = sum(1 for line in data.splitlines() if line.startswith(">"))
    return num_records


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

    # Descărcare FASTA
    n = download_fasta(
        args.email,
        out_path,
        query=args.query,
        accession=args.accession,
        db=args.db,
        retmax=args.retmax,
        api_key=args.api_key
    )
    print(f"[ok] Am scris {n} înregistrări în: {out_path}")

    # Citire FASTA local
    records = list(SeqIO.parse(str(out_path), "fasta"))
    if not records:
        print("[!] Fișierul FASTA este gol sau invalid.")
        return

    # Calcul GC pentru fiecare
    for rec in records:
        gc = gc_fraction(str(rec.seq))
        print(f"{rec.id}\tGC={gc:.3f}")


if __name__ == "__main__":
    main()
