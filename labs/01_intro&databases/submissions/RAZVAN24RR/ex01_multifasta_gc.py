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
    Descarcă secvențe FASTA din NCBI folosind Entrez.
    
    Args:
        email: Adresa de email obligatorie pentru NCBI
        out_path: Calea către fișierul FASTA de ieșire
        query: Query string pentru căutare (ex: "TP53[Gene] AND Homo sapiens[Organism]")
        accession: Număr de accesiune direct (ex: NM_000546)
        db: Baza de date NCBI (nuccore sau protein)
        retmax: Numărul maxim de rezultate
        api_key: Cheie API NCBI (opțional, crește rata de request-uri)
    
    Returns:
        Numărul de înregistrări scrise în fișier
    """
    # Configurare Entrez
    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key
    
    id_list = []
    
    # Cazul 1: Avem un număr de accesiune direct
    if accession:
        print(f"[info] Descărcare directă pentru accession: {accession}")
        id_list = [accession]
    
    # Cazul 2: Avem un query de căutare
    elif query:
        print(f"[info] Căutare după query: {query}")
        try:
            # Efectuăm căutarea cu esearch
            handle = Entrez.esearch(db=db, term=query, retmax=retmax)
            record = Entrez.read(handle)
            handle.close()
            
            id_list = record["IdList"]
            count = int(record["Count"])
            
            print(f"[info] Găsite {count} rezultate, descărcăm {len(id_list)}")
            
            if not id_list:
                print("[warning] Nu s-au găsit rezultate pentru acest query")
                return 0
                
        except Exception as e:
            print(f"[error] Eroare la căutare: {e}", file=sys.stderr)
            raise
    
    else:
        print("[error] Trebuie specificat fie --query fie --accession", file=sys.stderr)
        return 0
    
    # Descărcăm secvențele cu efetch
    try:
        print(f"[info] Descărcare {len(id_list)} secvențe din {db}...")
        
        # Convertim id_list la string pentru efetch
        if isinstance(id_list[0], int) or id_list[0].isdigit():
            # ID-uri numerice din esearch
            id_string = ",".join(str(i) for i in id_list)
        else:
            # Accession numbers
            id_string = ",".join(id_list)
        
        handle = Entrez.efetch(
            db=db,
            id=id_string,
            rettype="fasta",
            retmode="text"
        )
        
        # Scriem direct în fișier
        fasta_data = handle.read()
        handle.close()
        
        with open(out_path, "w") as f:
            f.write(fasta_data)
        
        # Numărăm câte înregistrări am scris
        num_records = fasta_data.count(">")
        
        return num_records
        
    except Exception as e:
        print(f"[error] Eroare la descărcare: {e}", file=sys.stderr)
        raise


def main():
    ap = argparse.ArgumentParser(
        description="Descarcă secvențe FASTA de la NCBI și calculează fracția GC"
    )
    ap.add_argument("--email", required=True, help="Email obligatoriu pentru NCBI Entrez")
    ap.add_argument("--api_key", help="NCBI API key (opțional)")
    ap.add_argument("--query", help="Ex: 'TP53[Gene] AND Homo sapiens[Organism]'")
    ap.add_argument("--accession", help="Ex: NM_000546")
    ap.add_argument("--db", default="nuccore", choices=["nuccore", "protein"])
    ap.add_argument("--retmax", type=int, default=3)
    ap.add_argument("--out", required=True, help="Fișier FASTA de ieșire")
    args = ap.parse_args()

    # Validare: trebuie specificat fie query fie accession
    if not args.query and not args.accession:
        ap.error("Trebuie specificat fie --query fie --accession")
    
    if args.query and args.accession:
        ap.error("Specificați doar unul dintre --query sau --accession, nu ambele")

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    # Descărcăm secvențele
    try:
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
    except Exception as e:
        print(f"[error] Descărcarea a eșuat: {e}", file=sys.stderr)
        sys.exit(1)
    
    # Verificăm dacă fișierul există și nu e gol
    if not out_path.exists() or out_path.stat().st_size == 0:
        print("[error] Fișierul de ieșire nu există sau este gol", file=sys.stderr)
        sys.exit(1)

    # Citim fișierul FASTA și calculăm GC
    print("\n[info] Calculare fracție GC pentru fiecare secvență:")
    print("-" * 60)
    
    try:
        records = SeqIO.parse(out_path, "fasta")
        count = 0
        
        for rec in records:
            seq_str = str(rec.seq)
            gc_val = gc_fraction(seq_str)
            print(f"{rec.id}\tGC={gc_val:.3f}")
            count += 1
        
        if count == 0:
            print("[warning] Nu s-au găsit înregistrări în fișierul FASTA")
        else:
            print("-" * 60)
            print(f"[ok] Total: {count} secvențe procesate")
            
    except Exception as e:
        print(f"[error] Eroare la citirea/procesarea fișierului FASTA: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
