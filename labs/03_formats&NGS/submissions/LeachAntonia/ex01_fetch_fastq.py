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
import sys
import os
import requests

def main():
    # TODO: citiți accession-ul (ex. sys.argv)
    if len(sys.argv) < 2:
        print("Usage: python ex01_fetch_fastq.py <ACCESSION>")
        sys.exit(1)

    accession = sys.argv[1].strip()

    # TODO: interogați sursa (ENA/SRA) pentru link FASTQ
    ena_api = "https://www.ebi.ac.uk/ena/portal/api/filereport"
    params = {
        "accession": accession,
        "result": "read_run",
        "fields": "fastq_ftp",
        "format": "tsv"
    }

    response = requests.get(ena_api, params=params)
    if response.status_code != 200:
        print(f"Accession '{accession}' not found.")
        sys.exit(1)

    lines = response.text.strip().split("\n")
    if len(lines) < 2:
        print(f"FASTQ for '{accession}' not found.")
        sys.exit(1)

    # TODO: descărcați fișierul în Locația ALEASĂ DE VOI
    ftp_links = lines[1].split("\t")[-1].split(";")
    first_link = ftp_links[0]
    if not first_link.startswith("ftp://"):
        first_link = "https://" + first_link.replace("ftp://", "ftp.")

    out_dir = os.path.expanduser("../../../../data/work/LeachAntonia/lab03")
    os.makedirs(out_dir, exist_ok=True)

    file_name = os.path.basename(first_link)
    out_path = os.path.join(out_dir, file_name)

    print(f"Downloading {accession} -> {out_path}")

    with requests.get(first_link, stream=True) as r:
        r.raise_for_status()
        with open(out_path, "wb") as f:
            for chunk in r.iter_content(chunk_size=8192):
                f.write(chunk)

    # TODO: print("Downloaded:", <cale_fisier>)
    print("Downloaded:", out_path)


if __name__ == "__main__":
    main()