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
"""

import argparse
from pathlib import Path
import sys
import time
from io import StringIO

from Bio import SeqIO
from Bio import Entrez  # deblocat și folosit pentru descărcare


def gc_fraction(seq: str) -> float:
    """Fracție GC pentru o secvență; robust la litere mici/mari și non-ATGC."""
    s = seq.upper()
    atgc = [c for c in s if c in ("A", "T", "G", "C")]
    if not atgc:
        return 0.0
    g = atgc.count("G")
    c = atgc.count("C")
    return (g + c) / float(len(atgc))


def _entrez_read_safe(handle, retries=3, backoff=1.5):
    """Wrapper mic cu retry pentru Entrez.read (uneori mai dă erori tranzitorii)."""
    last = None
    for i in range(retries):
        try:
            return Entrez.read(handle)
        except Exception as e:
            last = e
            if i == retries - 1:
                raise
            time.sleep(backoff ** i)
    raise last


def download_fasta(email: str, out_path: Path, query: str = None,
                   accession: str = None, db: str = "nuccore",
                   retmax: int = 3, api_key: str = None) -> int:
    """
    Implementați descărcarea din NCBI.
    Pași:
      - Configurați Entrez cu email (și api_key opțional).
      - Dacă avem accession: descărcați acel record.
      - Altfel, dacă avem query: faceți esearch -> lista de id-uri, apoi efetch.
      - Scrieți rezultatele în out_path.
      - Returnați numărul de înregistrări scrise.
    """
    if not email:
        raise ValueError("Email este obligatoriu pentru NCBI Entrez (--email).")

    Entrez.email = email
    Entrez.tool = "lab01_multifasta_gc"
    if api_key:
        Entrez.api_key = api_key

    if accession and query:
        raise ValueError("Specificați fie --accession, fie --query (nu ambele).")
    if not accession and not query:
        raise ValueError("Trebuie să furnizați --accession sau --query.")

    # 1) Obținere ID-uri
    if accession:
        ids = [accession]
    else:
        with Entrez.esearch(db=db, term=query, retmax=retmax) as h:
            result = _entrez_read_safe(h)
        ids = result.get("IdList", [])
        if not ids:
            raise RuntimeError(f"Nu s-au găsit rezultate pentru query: {query!r}")

    # 2) Descărcare FASTA
    with Entrez.efetch(db=db, id=",".join(ids), rettype="fasta", retmode="text") as h:
        fasta_text = h.read()

    if not fasta_text.strip():
        raise RuntimeError("Descărcarea a returnat conținut gol.")

    # 3) Scriere pe disc
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(fasta_text, encoding="utf-8")

    # 4) Număr înregistrări
    n_records = sum(1 for line in fasta_text.splitlines() if line.startswith(">"))
    if n_records == 0:
        # fallback: parsăm în memorie
        n_records = sum(1 for _ in SeqIO.parse(StringIO(fasta_text), "fasta"))

    return n_records


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

    # Descărcare și scriere FASTA
    try:
        n = download_fasta(args.email, out_path, query=args.query,
                           accession=args.accession, db=args.db,
                           retmax=args.retmax, api_key=args.api_key)
        print(f"[ok] Am scris {n} înregistrări în: {out_path}", file=sys.stderr)
    except Exception as e:
        print(f"[eroare] {e}", file=sys.stderr)
        sys.exit(1)

    # Citire FASTA și calcul GC
    try:
        records = list(SeqIO.parse(str(out_path), "fasta"))
        if args.db == "protein":
            print("[atenție] db='protein': GC se calculează doar pe A/T/G/C din șirul proteic; "
                  "de regulă iese 0.000.", file=sys.stderr)
        for rec in records:
            val = gc_fraction(str(rec.seq))
            print(f"{rec.id}\tGC={val:.3f}")
    except Exception as e:
        print(f"[eroare] Nu am putut citi/parse fișierul FASTA: {e}", file=sys.stderr)
        sys.exit(2)


if __name__ == "__main__":
    main()
