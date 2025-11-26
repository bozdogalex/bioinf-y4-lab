#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Exercițiu (Lab 1): Descărcare FASTA + calcul GC
"""

import argparse
from pathlib import Path
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
    """
    Descărcare FASTA din NCBI (prin Entrez).
    """
    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key

    # Dacă avem accession, descărcăm direct acel record
    if accession:
        with Entrez.efetch(db=db, id=accession, rettype="fasta", retmode="text") as handle:
            data = handle.read()
        with open(out_path, "w") as f:
            f.write(data)
        records = list(SeqIO.parse(out_path, "fasta"))
        return len(records)

    # Dacă avem query, facem mai întâi căutarea
    elif query:
        with Entrez.esearch(db=db, term=query, retmax=retmax) as handle:
            result = Entrez.read(handle)
        id_list = result.get("IdList", [])
        if not id_list:
            print("[eroare] Nu am găsit rezultate pentru query.")
            return 0
        with Entrez.efetch(db=db, id=",".join(id_list), rettype="fasta", retmode="text") as handle:
            data = handle.read()
        with open(out_path, "w") as f:
            f.write(data)
        records = list(SeqIO.parse(out_path, "fasta"))
        return len(records)

    else:
        raise ValueError("Trebuie specificat fie --accession, fie --query.")


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

    # Descărcăm fișierul FASTA
    n = download_fasta(args.email, out_path, query=args.query,
                       accession=args.accession, db=args.db,
                       retmax=args.retmax, api_key=args.api_key)

    print(f"[ok] Am scris {n} înregistrări în: {out_path}")

    # Calculăm GC pentru fiecare secvență
    for rec in SeqIO.parse(out_path, "fasta"):
        gc = gc_fraction(str(rec.seq))
        print(f"{rec.id}\tGC={gc:.3f}")


if __name__ == "__main__":
    main()
