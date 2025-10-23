#!/usr/bin/env python3
import sys
import os
import requests

def main():
    if len(sys.argv) < 2:
        print("Usage: python download_fastq.py <ACCESSION>") #SRR390728
        sys.exit(1)

    accession = sys.argv[1].strip()

    ena_url = f"https://www.ebi.ac.uk/ena/portal/api/filereport"
    params = {
        "accession": accession,
        "result": "read_run",
        "fields": "fastq_ftp",
        "format": "tsv"
    }

    print(f"[INFO] Fetching FASTQ links for {accession}...")

    response = requests.get(ena_url, params=params)
    if response.status_code != 200:
        print(f"[ERROR] ENA request failed with code {response.status_code}")
        sys.exit(1)

    lines = response.text.strip().split("\n")
    if len(lines) < 2:
        print(f"[ERROR] No FASTQ links found for accession {accession}")
        sys.exit(1)

    fastq_links = lines[1].split("\t")[1]
    fastq_urls = ["https://" + link for link in fastq_links.split(";") if link]
    fastq_url = fastq_urls[0]

    print(f"[INFO] Downloading {fastq_url} ...")

    output_dir = os.path.join("data", "work", "StanaAndrei", "sample")
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, "your_reads.fastq.gz")

    with requests.get(fastq_url, stream=True) as r:
        r.raise_for_status()
        with open(output_path, "wb") as f:
            for chunk in r.iter_content(chunk_size=8192):
                if chunk:
                    f.write(chunk)

    print("Downloaded:", os.path.abspath(output_path))


if __name__ == "__main__":
    main()
