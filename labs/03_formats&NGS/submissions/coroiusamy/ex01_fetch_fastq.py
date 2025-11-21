#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Exercițiu 03 — Descărcare FASTQ (student-owned)

Obiectiv:
- Alegeți un accession TP53-related (ex. SRR..., ERR...) și DESCĂRCAȚI un fișier FASTQ.
- Salvați in  data/work/<handle>/lab03/your_reads.fastq.gz

Cerințe minime:
- Scriptul trebuie să accepte un accession (de ex. prin arg linie de comandă).
- Scriptul descarcă cel puțin un FASTQ (un singur fișier e suficient pentru exercițiu).
- Scriptul afișează pe stdout calea fișierului descărcat.

Recomandat :
- Suportați .fastq sau .fastq.gz.

NOTĂ:
- Nu contează biblioteca aleasă (requests/urllib/etc.), dar evitați pachete grele.
"""

import urllib.request
import shutil
from pathlib import Path
import sys
import argparse  

def fetch_ena_fastq(accession: str, out_file: Path) -> bool:
    """
    Descarcă un fișier FASTQ de la ENA folosind API-ul filereport.
    """
    
    # Construim URL-ul API ENA
    ena_url = (
        f"https://www.ebi.ac.uk/ena/portal/api/filereport?"
        f"accession={accession}&result=read_run&"
        f"fields=fastq_ftp&format=tsv"
    )

    # mesaje progres
    print(f"[Info] Interogare ENA API pentru: {accession}...", file=sys.stderr)
    
    try:
        with urllib.request.urlopen(ena_url) as response:
            data = response.read().decode('utf-8')
    except Exception as e:
        print(f"[Eroare] ENA API request eșuat: {e}", file=sys.stderr)
        return False

    # parsare TSV 
    lines = data.strip().split('\n')
    if len(lines) < 2: 
        print(f"[Eroare] Răspuns ENA invalid (nu s-au găsit date) pentru {accession}", file=sys.stderr)
        return False

    try:
        header = lines[0].split('\t')
        values = lines[1].split('\t')
        
        # coloana fastq_ftp
        col_index = header.index('fastq_ftp')
        fastq_links_str = values[col_index]
        
        if not fastq_links_str:
            print(f"[Eroare] ENA nu a returnat un link 'fastq_ftp' pentru {accession}", file=sys.stderr)
            return False

        fastq_ftp_link = fastq_links_str.split(';')[0]
        
        # FTP -> HTTP 
        if fastq_ftp_link.startswith("ftp://"):
            http_link = "https://" + fastq_ftp_link[6:] # Taie primii 6 caractere ("ftp://")
        else:
            http_link = "https://" + fastq_ftp_link # Adaugă direct

    except (ValueError, IndexError) as e:
        print(f"[Eroare] Nu am putut parsa răspunsul ENA: {e}\nData: {data}", file=sys.stderr)
        return False

    # descarcare fisier
    print(f"[Info] Se descarcă de la: {http_link}", file=sys.stderr)
    
    try:
        with urllib.request.urlopen(http_link) as response_dl, open(out_file, 'wb') as f_out:
            shutil.copyfileobj(response_dl, f_out)
        return True # Succes
    except Exception as e:
        print(f"[Eroare] Descărcarea fișierului a eșuat: {e}", file=sys.stderr)
        return False


def main():
    # TODO: citiți accession-ul (ex. sys.argv)
    # parser-ul de argumente
    parser = argparse.ArgumentParser(description="Descarcă un fișier FASTQ de la ENA.")
    
    # argumentele cerute
    parser.add_argument("--accession", required=True, help="ID-ul Run ENA/SRA (ex: ERR179724)")
    parser.add_argument("--handle", required=True, help="Handle-ul tău GitHub (pt calea de output)")
    
    # Parsare argumente
    args = parser.parse_args()

    # TODO: descărcați fișierul în Locația ALEASĂ DE VOI
    # Construim calea de output
    out_path = Path(f"data/work/{args.handle}/lab03/your_reads.fastq.gz")
    
    # folder output
    out_path.parent.mkdir(parents=True, exist_ok=True)
    
    # TODO: interogați sursa (ENA/SRA) pentru link FASTQ
    # TODO: descărcați fișierul
    success = fetch_ena_fastq(args.accession, out_path)

    # TODO: print("Downloaded:", <cale_fisier>)
    # verificare descarcare
    if success:
        # cale fisier
        print(f"Downloaded: {out_path}")
    else:
        # mesaj eroare
        print(f"[Eroare] Descărcarea a eșuat pentru {args.accession}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()