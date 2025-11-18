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

def main():
    import sys
    import urllib.request
    from pathlib import Path

    if len(sys.argv) < 2:
        print("Utilizare: python ex01_fetch_fastq.py <accession>")
        sys.exit(1)

    acc = sys.argv[1]
    url = f"https://www.ebi.ac.uk/ena/portal/api/filereport?accession={acc}&result=read_run&fields=run_accession,fastq_ftp&format=tsv"

    try:
        with urllib.request.urlopen(url) as r:
            txt = r.read().decode("utf-8").strip().splitlines()
    except Exception as e:
        print(f"Eroare: {e}", file=sys.stderr)
        sys.exit(1)

    if len(txt) < 2:
        print(f"Eroare: nu am găsit rânduri pentru {acc}", file=sys.stderr)
        sys.exit(1)

    header = txt[0].split("\t")
    row = txt[1].split("\t")
    try:
        i = header.index("fastq_ftp")
        field = row[i].strip()
    except Exception:
        print("Eroare: câmpul fastq_ftp lipsește", file=sys.stderr)
        sys.exit(1)

    if not field:
        print(f"Eroare: fără link FASTQ pentru {acc}", file=sys.stderr)
        sys.exit(1)

    first = field.split(";")[0].strip()
    if first.startswith("ftp://"):
        fastq = first.replace("ftp://", "https://", 1)
    elif first.startswith("ftp.sra.ebi.ac.uk"):
        fastq = "https://" + first
    else:
        fastq = first

    dest_dir = Path("data/work/rauldavid35/lab03")
    dest_dir.mkdir(parents=True, exist_ok=True)
    filename = fastq.split("/")[-1]
    outpath = dest_dir / filename

    try:
        print(f"Descărcam {acc}...", file=sys.stderr)
        urllib.request.urlretrieve(fastq, outpath)
    except Exception as e:
        print(f"Eroare la descărcare: {e}", file=sys.stderr)
        sys.exit(1)

    print(f"Downloaded: {outpath}")


if __name__ == "__main__":
    main()
