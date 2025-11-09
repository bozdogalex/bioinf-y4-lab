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
import urllib.request
import os

def main():
    # TODO: citiți accession-ul (ex. sys.argv)
    if len(sys.argv) < 2:
        print("Usage: python script.py <accession>")
        return

    accession = sys.argv[1]
    api_url = f"https://www.ebi.ac.uk/ena/portal/api/filereport?accession={accession}&result=read_run&fields=fastq_ftp"

    # interogați sursa (ENA/SRA) pentru link FASTQ
    try:
        with urllib.request.urlopen(api_url) as response:
            data = response.read().decode('utf-8')
    except Exception as e:
        print(f"Error querying API: {e}")
        return

    lines = data.splitlines()
    if len(lines) < 2:
        print("No data found")
        return

    fields = lines[1].split('\t')
    if len(fields) < 2:
        print("Invalid data format")
        return

    fastq_ftps = fields[1].split(';')
    if not fastq_ftps or not fastq_ftps[0]:
        print("No FASTQ found")
        return

    first_fastq_path = fastq_ftps[0]
    first_fastq = "ftp://" + first_fastq_path

    # descărcați fișierul în Locația ALEASĂ DE VOI
    handle = "banmepls"
    dir_path = f"data/work/{handle}/lab03"
    os.makedirs(dir_path, exist_ok=True)

    if first_fastq_path.endswith('.fastq.gz'):
        file_name = "your_reads.fastq.gz"
    elif first_fastq_path.endswith('.fastq'):
        file_name = "your_reads.fastq"
    else:
        print("Unsupported file type")
        return

    file_path = os.path.join(dir_path, file_name)

    # print("Downloaded:", <cale_fisier>)
    try:
        urllib.request.urlretrieve(first_fastq, file_path)
        print("Downloaded:", file_path)
    except Exception as e:
        print(f"Error downloading file: {e}")


if __name__ == "__main__":
    main()
