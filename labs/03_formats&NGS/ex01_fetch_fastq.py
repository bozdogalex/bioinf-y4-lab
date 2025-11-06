"""
Exercițiu 03 — Descărcare FASTQ (student-owned)

Obiectiv:
- Alegeți un accession TP53-related (ex. SRR..., ERR...) și DESCĂRCAȚI un fișier FASTQ.
- Salvați in  data/work/<handle>/lab03/your_reads.fastq.gz

Cerințe minime:
- Scriptul trebuie să accepte un accession (de ex. prin arg linie de comandă).
- Scriptul descarcă cel puțin un FASTQ (un singur fișier e suficient pentru exercițiu).
- Scriptul afișează pe stdout calea fișierului descărcat.

Recomandat :
- Suportați .fastq sau .fastq.gz.

NOTĂ:
- Nu contează biblioteca aleasă (requests/urllib/etc.), dar evitați pachete grele.
"""

import argparse
from pathlib import Path
import sys
import os
import requests

def main():
    # TODO: citiți accession-ul (ex. sys.argv)
    # TODO: interogați sursa (ENA/SRA) pentru link FASTQ
    # TODO: descărcați fișierul în Locația ALEASĂ DE VOI
    # TODO: print("Downloaded:", <cale_fisier>)
    ap = argparse.ArgumentParser()
    ap.add_argument("--accession", help="Ex: NM_000546")
    args = ap.parse_args()
    accession = args.accession
    
    out_dir = Path("data/work/marabonea16/lab03")
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"{accession}_reads.fastq.gz"

    ena_meta = f"https://www.ebi.ac.uk/ena/portal/api/reads?accession={accession}&result=read_run&fields=fastq_ftp,fastq_http"
    


if __name__ == "__main__":
    main()
