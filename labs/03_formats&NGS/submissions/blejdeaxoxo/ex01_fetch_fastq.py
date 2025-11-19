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
import urllib.request
import urllib.parse
import csv
import io
import gzip
import shutil
import argparse

def query_ena_fastq(accession):
    """Return first fastq URL (ftp or http) from ENA filereport for accession, or None."""
    q = (
        "https://www.ebi.ac.uk/ena/portal/api/filereport"
        "?accession={acc}&result=read_run&fields=fastq_ftp&format=tsv"
    ).format(acc=urllib.parse.quote(accession))
    req = urllib.request.Request(q, headers={"User-Agent": "python-urllib/3"})
    with urllib.request.urlopen(req, timeout=30) as resp:
        txt = resp.read().decode("utf-8")
    if not txt.strip():
        return None
    rdr = csv.DictReader(io.StringIO(txt), delimiter="\t")
    for row in rdr:
        val = row.get("fastq_ftp") or ""
        if not val:
            continue
        parts = [p.strip() for p in val.split(";") if p.strip()]
        if parts:
            first = parts[0]
            if first.startswith("ftp://") or first.startswith("http://") or first.startswith("https://"):
                return first
            return "https://" + first
    return None

def download_stream_to_path(url, out_path, compress_if_needed=True):
    req = urllib.request.Request(url, headers={"User-Agent": "python-urllib/3"})
    with urllib.request.urlopen(req, timeout=60) as resp:
        tmp_path = out_path + ".tmp"
        if compress_if_needed and not out_path.endswith(".gz"):
            with gzip.open(tmp_path, "wb") as f:
                shutil.copyfileobj(resp, f)
        else:
            with open(tmp_path, "wb") as f:
                shutil.copyfileobj(resp, f)
        os.replace(tmp_path, out_path)

def main():
    parser = argparse.ArgumentParser(description="Download one FASTQ from ENA for given accession")
    parser.add_argument("accession", help="Run accession (e.g. SRR..., ERR...)")
    parser.add_argument("--handle", help="Your handle (used for output path)", default=os.environ.get("USER", "me"))
    args = parser.parse_args()

    accession = args.accession
    handle = args.handle

    url = query_ena_fastq(accession)
    if not url:
        print(f"ERROR: No FASTQ URL found for accession {accession}", file=sys.stderr)
        sys.exit(2)

    out_dir = os.path.join("data", "work", handle, "lab03")
    os.makedirs(out_dir, exist_ok=True)
    fname = os.path.basename(urllib.parse.urlparse(url).path)
    if not fname:
        fname = f"{accession}.fastq.gz"
    if fname.endswith(".gz"):
        out_name = fname
    else:
        out_name = fname + ".gz"
    out_path = os.path.join(out_dir, out_name)

    try:
        download_stream_to_path(url, out_path, compress_if_needed=not fname.endswith(".gz"))
    except Exception as e:
        print(f"ERROR: download failed: {e}", file=sys.stderr)
        sys.exit(3)

    print("Downloaded:", out_path)



if __name__ == "__main__":
    main()
