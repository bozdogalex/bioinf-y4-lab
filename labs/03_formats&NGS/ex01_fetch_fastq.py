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
import requests

def get_fastq_url(accession):
    ena_api = "https://www.ebi.ac.uk/ena/portal/api/filereport"
    params = {
        "accession": accession,
        "result": "read_run",
        "fields": "fastq_ftp",
        "format": "tsv"
    }
    response = requests.get(ena_api, params=params)
    lines = response.text.strip().split('\n')
    if len(lines) < 2:
        raise Exception("No FASTQ found.")
    ftp_links = lines[1].split('\t')[1].split(';')
    http_link = "https://" + ftp_links[0]
    return http_link


def download_fastq(url, output_path):
    r = requests.get(url, stream=True)
    with open(output_path, 'wb') as f:
        for chunk in r.iter_content(chunk_size=8192):
            f.write(chunk)

if __name__ == "__main__":
    accession = sys.argv[1]
    urls = get_fastq_url(accession)
    print(f"urls:", urls)
    output_path = "/workspaces/bioinf-y4-lab/data/work/MariusJalba/lab03/your_reads.fastq.gz"
    download_fastq(urls, output_path)
    print("your_reads.fastq.gz")

