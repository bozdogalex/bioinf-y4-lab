#!/usr/bin/env python3
"""
Usage:
    python labs/03_formats\&NGS/submissions/Razvann19/ex01_fetch_fastq.py --accession SRR16599692 \
  --out data/work/Razvann19/lab03/your_reads.fastq.gz
"""

import argparse
import urllib.request
import json
from pathlib import Path


def fetch_fastq_url(accession: str) -> str:
    """Returnează primul link FASTQ de la ENA pentru accesionul dat."""

    api = f"https://www.ebi.ac.uk/ena/portal/api/filereport"
    params = (
        f"?accession={accession}"
        f"&result=read_run"
        f"&fields=fastq_ftp"
        f"&format=json"
    )

    url = api + params
    print(f"[INFO] Querying ENA: {url}")

    with urllib.request.urlopen(url) as response:
        data = json.loads(response.read())

    if not data or "fastq_ftp" not in data[0] or not data[0]["fastq_ftp"]:
        raise ValueError(f"[ERR] Nu s-au găsit linkuri FASTQ pentru {accession}")

    # ENA returnează linkuri ftp, le convertim în https
    fastq_links = data[0]["fastq_ftp"].split(";")
    first_link = "https://" + fastq_links[0]

    print(f"[INFO] Found FASTQ: {first_link}")
    return first_link


def download_file(url: str, out_path: Path):
    """Descarcă fișier prin streaming."""
    print(f"[INFO] Downloading to: {out_path}")

    out_path.parent.mkdir(parents=True, exist_ok=True)

    with urllib.request.urlopen(url) as response, open(out_path, "wb") as f:
        f.write(response.read())

    print(f"[DONE] Downloaded: {out_path}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--accession", required=True, help="SRR/ERR accession ID")
    ap.add_argument("--out", required=True, help="Output FASTQ path (.fastq.gz)")
    args = ap.parse_args()

    acc = args.accession
    out = Path(args.out)

    fastq_url = fetch_fastq_url(acc)
    download_file(fastq_url, out)

    print("Downloaded:", out)


if __name__ == "__main__":
    main()
