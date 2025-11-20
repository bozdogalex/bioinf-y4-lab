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
import typing

from Bio import SeqIO
from Bio import Entrez
from Bio import BiopythonWarning
import warnings

warnings.simplefilter("ignore", BiopythonWarning)


def gc_fraction(seq: typing.Union[str, object]) -> float:
    """Fracție GC pentru o secvență; robust la litere mici/mari și non-ATGC."""
    s = str(seq).upper()
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
    Descarcă secvențe FASTA de la NCBI și le salvează în out_path.
    Returnează numărul de înregistrări scrise.
    """
    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key

    ids = []
    try:
        if accession:
            # accept comma-separated list or single accession
            ids = [a.strip() for a in accession.split(",") if a.strip()]
            # If accessions are not numeric IDs, efetch by accession still works when
            # rettype="fasta" and retmode="text".
            with Entrez.efetch(db=db, id=",".join(ids), rettype="fasta", retmode="text") as handle:
                records = list(SeqIO.parse(handle, "fasta"))
        elif query:
            with Entrez.esearch(db=db, term=query, retmax=retmax) as handle:
                search_results = Entrez.read(handle)
                ids = search_results.get("IdList", [])
            if not ids:
                return 0
            with Entrez.efetch(db=db, id=",".join(ids), rettype="fasta", retmode="text") as handle:
                records = list(SeqIO.parse(handle, "fasta"))
        else:
            raise ValueError("Trebuie specificat --query sau --accession")
    except Exception as e:
        raise RuntimeError(f"Entrez download failed: {e}")

    # write records using SeqIO to ensure valid FASTA formatting
    written = 0
    try:
        with out_path.open("w", encoding="utf-8") as fh:
            written = SeqIO.write(records, fh, "fasta")
    except Exception as e:
        raise RuntimeError(f"Failed to write output file {out_path}: {e}")

    return written


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--email", required=True, help="Email obligatoriu pentru NCBI Entrez")
    ap.add_argument("--api_key", help="NCBI API key (opțional)")
    ap.add_argument("--query", help="Ex: 'TP53[Gene] AND Homo sapiens[Organism]'")
    ap.add_argument("--accession", help="Ex: NM_000546 (comma-separated allowed)")
    ap.add_argument("--db", default="nuccore", choices=["nuccore", "protein"])
    ap.add_argument("--retmax", type=int, default=3)
    ap.add_argument("--out", required=True, help="Fișier FASTA de ieșire")
    args = ap.parse_args()

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    try:
        n = download_fasta(
            email=args.email,
            out_path=out_path,
            query=args.query,
            accession=args.accession,
            db=args.db,
            retmax=args.retmax,
            api_key=args.api_key
        )
    except Exception as e:
        print(f"[error] {e}", file=sys.stderr)
        sys.exit(1)

    print(f"[ok] Am scris {n} înregistrări în: {out_path}")

    # Citire FASTA și calcul GC
    try:
        records = list(SeqIO.parse(str(out_path), "fasta"))
    except Exception as e:
        print(f"[error] could not parse FASTA: {e}", file=sys.stderr)
        sys.exit(1)

    if not records:
        print("[warn] no records found in output FASTA", file=sys.stderr)
        return

    for rec in records:
        gc = gc_fraction(rec.seq)
        print(f"{rec.id}\tGC={gc:.3f}")


if __name__ == "__main__":
    main()
# filepath: /workspaces/bioinf-y4-lab/labs/01_intro&databases/ex01_multifasta_gc.py
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
import typing

from Bio import SeqIO
from Bio import Entrez
from Bio import BiopythonWarning
import warnings

warnings.simplefilter("ignore", BiopythonWarning)


def gc_fraction(seq: typing.Union[str, object]) -> float:
    """Fracție GC pentru o secvență; robust la litere mici/mari și non-ATGC."""
    s = str(seq).upper()
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
    Descarcă secvențe FASTA de la NCBI și le salvează în out_path.
    Returnează numărul de înregistrări scrise.
    """
    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key

    ids = []
    try:
        if accession:
            # accept comma-separated list or single accession
            ids = [a.strip() for a in accession.split(",") if a.strip()]
            # If accessions are not numeric IDs, efetch by accession still works when
            # rettype="fasta" and retmode="text".
            with Entrez.efetch(db=db, id=",".join(ids), rettype="fasta", retmode="text") as handle:
                records = list(SeqIO.parse(handle, "fasta"))
        elif query:
            with Entrez.esearch(db=db, term=query, retmax=retmax) as handle:
                search_results = Entrez.read(handle)
                ids = search_results.get("IdList", [])
            if not ids:
                return 0
            with Entrez.efetch(db=db, id=",".join(ids), rettype="fasta", retmode="text") as handle:
                records = list(SeqIO.parse(handle, "fasta"))
        else:
            raise ValueError("Trebuie specificat --query sau --accession")
    except Exception as e:
        raise RuntimeError(f"Entrez download failed: {e}")

    # write records using SeqIO to ensure valid FASTA formatting
    written = 0
    try:
        with out_path.open("w", encoding="utf-8") as fh:
            written = SeqIO.write(records, fh, "fasta")
    except Exception as e:
        raise RuntimeError(f"Failed to write output file {out_path}: {e}")

    return written


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--email", required=True, help="Email obligatoriu pentru NCBI Entrez")
    ap.add_argument("--api_key", help="NCBI API key (opțional)")
    ap.add_argument("--query", help="Ex: 'TP53[Gene] AND Homo sapiens[Organism]'")
    ap.add_argument("--accession", help="Ex: NM_000546 (comma-separated allowed)")
    ap.add_argument("--db", default="nuccore", choices=["nuccore", "protein"])
    ap.add_argument("--retmax", type=int, default=3)
    ap.add_argument("--out", required=True, help="Fișier FASTA de ieșire")
    args = ap.parse_args()

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    try:
        n = download_fasta(
            email=args.email,
            out_path=out_path,
            query=args.query,
            accession=args.accession,
            db=args.db,
            retmax=args.retmax,
            api_key=args.api_key
        )
    except Exception as e:
        print(f"[error] {e}", file=sys.stderr)
        sys.exit(1)

    print(f"[ok] Am scris {n} înregistrări în: {out_path}")

    # Citire FASTA și calcul GC
    try:
        records = list(SeqIO.parse(str(out_path), "fasta"))
    except Exception as e:
        print(f"[error] could not parse FASTA: {e}", file=sys.stderr)
        sys.exit(1)

    if not records:
        print("[warn] no records found in output FASTA", file=sys.stderr)
        return

    for rec in records:
        gc = gc_fraction(rec.seq)
        print(f"{rec.id}\tGC={gc:.3f}")


if __name__ == "__main__":
    main()