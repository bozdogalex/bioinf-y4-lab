"""
Exercițiu 03 — Descărcare FASTQ (student-owned)
python labs/03_formats\&NGS/ex01_fetch_fastq.py --accession SRR390728 --out data/work/marabonea16/lab03/your_reads.fastq


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

from anyio import Path
import argparse
from pathlib import Path
import urllib.request
import gzip
import shutil


def main():
    # TODO: citiți accession-ul (ex. sys.argv)
    ap = argparse.ArgumentParser()
    ap.add_argument("--accession", required=True, help="Ex: SRR390728")
    ap.add_argument("--out", required=True, help="Fișier FASTQ de ieșire")
    args = ap.parse_args()
    
    accession = args.accession
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    # TODO: interogați sursa (ENA/SRA) pentru link FASTQ
    ena_url = f"https://www.ebi.ac.uk/ena/portal/api/filereport?accession={accession}&result=read_run&fields=run_accession,fastq_ftp&format=tsv"
    with urllib.request.urlopen(ena_url) as response:
            data = response.read().decode('utf-8')
            lines = data.strip().split('\n')
    fastq_links = lines[1].split('\t')[1]
    fastq_link = "ftp://" + fastq_links.split(";")[0]

    print(f"Fetching FASTQ from: {fastq_link}")

    # TODO: descărcați fișierul în Locația ALEASĂ DE VOI
    with urllib.request.urlopen(fastq_link) as response:
        if fastq_link.endswith(".gz"):
            with gzip.open(response, 'rb') as f_in:
                with open(out_path, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
        else:
            shutil.copyfileobj(response, f_out)


    # TODO: print("Downloaded:", <cale_fisier>)
    print("Downloaded:", out_path)
    pass


if __name__ == "__main__":
    main()
