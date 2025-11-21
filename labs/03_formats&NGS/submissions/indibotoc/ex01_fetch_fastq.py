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
import os, sys, requests

def main():
    # TODO: citiți accession-ul (ex. sys.argv)
    # TODO: interogați sursa (ENA/SRA) pentru link FASTQ
    handle = "student1"
    out_dir = f"data/work/indibotoc/lab03"
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, "your_reads.fastq.gz")

    # TODO: descărcați fișierul în Locația ALEASĂ DE VOI
    ena_api = f"https://www.ebi.ac.uk/ena/portal/api/filereport"
    params = {
        "accession": accession,
        "result": "read_run",
        "fields": "fastq_ftp",
        "format": "tsv"
    }
    response = requests.get(ena_api, params=params)
    response.raise_for_status()

    lines = response.text.strip().split("\n")
    if len(lines) < 2:
        print("Eroare: accession invalid sau fără FASTQ disponibil.")
        sys.exit(1)

    fastq_url = "https://" + lines[1].split("\t")[1].split(";")[0]  # doar primul fișier
    print(f"Downloading from: {fastq_url}")

    with requests.get(fastq_url, stream=True) as r:
        r.raise_for_status()
        with open(out_path, "wb") as f:
            for chunk in r.iter_content(chunk_size=8192):
                if chunk:
                    f.write(chunk)

    # TODO: print("Downloaded:", <cale_fisier>)
    print("Downloaded:", os.path.abspath(out_path))


if __name__ == "__main__":
    main()