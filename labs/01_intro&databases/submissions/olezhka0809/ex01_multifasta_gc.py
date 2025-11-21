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
        email: Email obligatoriu pentru NCBI
        out_path: Path pentru fișierul de output
        query: Query de căutare (ex: "TP53[Gene] AND Homo sapiens[Organism]")
        accession: Accession number specific (ex: NM_000546)
        db: Baza de date NCBI (nuccore sau protein)
        retmax: Număr maxim de rezultate
        api_key: API key NCBI (opțional, pentru rate limiting mai bun)
    
    Returns:
        Numărul de înregistrări descărcate
    """
    # Configurare Entrez
    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key
    
    try:
        # Caz 1: Descărcare direct prin accession number
        if accession:
            print(f"[info] Descărcare accession: {accession} din {db}")
            handle = Entrez.efetch(
                db=db,
                id=accession,
                rettype="fasta",
                retmode="text"
            )
            fasta_data = handle.read()
            handle.close()
            
            # Scriere în fișier
            with open(out_path, "w") as f:
                f.write(fasta_data)
            
            # Numără înregistrările
            num_records = fasta_data.count(">")
            return num_records
        
        # Caz 2: Căutare cu query și descărcare rezultate
        elif query:
            print(f"[info] Căutare query: '{query}' în {db}")
            
            # Pasul 1: esearch pentru a obține lista de ID-uri
            search_handle = Entrez.esearch(
                db=db,
                term=query,
                retmax=retmax,
                usehistory="y"
            )
            search_results = Entrez.read(search_handle)
            search_handle.close()
            
            id_list = search_results.get("IdList", [])
            count = int(search_results.get("Count", 0))
            
            print(f"[info] Găsite {count} rezultate, se descarcă primele {len(id_list)}")
            
            if not id_list:
                print("[warning] Nu s-au găsit rezultate pentru acest query!")
                # Creează fișier gol
                out_path.touch()
                return 0
            
            # Pasul 2: efetch pentru a descărca secvențele
            fetch_handle = Entrez.efetch(
                db=db,
                id=",".join(id_list),
                rettype="fasta",
                retmode="text"
            )
            fasta_data = fetch_handle.read()
            fetch_handle.close()
            
            # Scriere în fișier
            with open(out_path, "w") as f:
                f.write(fasta_data)
            
            return len(id_list)
        
        else:
            print("[error] Trebuie să furnizați --query sau --accession!", file=sys.stderr)
            sys.exit(1)
    
    except Exception as e:
        print(f"[error] Eroare la descărcarea din NCBI: {e}", file=sys.stderr)
        sys.exit(1)


def main():
    ap = argparse.ArgumentParser(
        description="Descarcă secvențe FASTA din NCBI și calculează conținutul GC",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Exemple de utilizare:

  # Descărcare prin query
  python ex01_multifasta_gc.py --email student@example.com \
    --query "TP53[Gene] AND Homo sapiens[Organism]" \
    --retmax 3 \
    --out data/work/olezhka0809/lab01/my_tp53.fa

  # Descărcare prin accession
  python ex01_multifasta_gc.py --email student@example.com \
    --accession NM_000546 \
    --out data/work/olezhka0809/lab01/nm000546.fa

  # Descărcare proteine
  python ex01_multifasta_gc.py --email student@example.com \
    --query "TP53[Gene] AND Homo sapiens[Organism]" \
    --db protein \
    --retmax 5 \
    --out data/work/olezhka0809/lab01/tp53_proteins.fa
        """
    )
    
    ap.add_argument("--email", required=True, help="Email obligatoriu pentru NCBI Entrez")
    ap.add_argument("--api_key", help="NCBI API key (opțional)")
    ap.add_argument("--query", help="Ex: 'TP53[Gene] AND Homo sapiens[Organism]'")
    ap.add_argument("--accession", help="Ex: NM_000546")
    ap.add_argument("--db", default="nuccore", choices=["nuccore", "protein"])
    ap.add_argument("--retmax", type=int, default=3)
    ap.add_argument("--out", required=True, help="Fișier FASTA de ieșire")
    args = ap.parse_args()

    # Validare: trebuie să avem query SAU accession
    if not args.query and not args.accession:
        print("[error] Trebuie să furnizați --query sau --accession!", file=sys.stderr)
        ap.print_help()
        sys.exit(1)

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    print("=" * 80)
    print("NCBI FASTA DOWNLOADER & GC CALCULATOR")
    print("=" * 80)

    # Pasul 1: Descărcare din NCBI
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

    # Pasul 2: Citire fișier FASTA local
    if not out_path.exists() or out_path.stat().st_size == 0:
        print("[warning] Fișierul este gol sau nu există!")
        return

    print("\n" + "=" * 80)
    print("CALCUL CONȚINUT GC")
    print("=" * 80)

    records = SeqIO.parse(out_path, "fasta")

    # Pasul 3: Calculare și afișare GC pentru fiecare secvență
    total_records = 0
    total_gc = 0.0
    
    for rec in records:
        gc = gc_fraction(str(rec.seq))
        total_gc += gc
        total_records += 1
        
        print(f"{rec.id}\tGC={gc:.3f}\tLen={len(rec.seq)}")
    
    # Statistici finale
    if total_records > 0:
        avg_gc = total_gc / total_records
        print("\n" + "-" * 80)
        print(f"Total secvențe: {total_records}")
        print(f"GC mediu: {avg_gc:.3f}")
        print("=" * 80)
    else:
        print("[warning] Nu s-au găsit secvențe în fișier!")


if __name__ == "__main__":
    main()