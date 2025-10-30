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
    import os
    from pathlib import Path
    
    # TODO: citiți accession-ul (ex. sys.argv)
    if len(sys.argv) < 2:
        print("Utilizare: python ex01_fetch_fastq.py <accession>")
        print("Exemplu: python ex01_fetch_fastq.py SRR123456")
        sys.exit(1)
    
    accession = sys.argv[1]
    
    # TODO: interogați sursa (ENA/SRA) pentru link FASTQ
    # Construim URL-ul pentru ENA (European Nucleotide Archive)
    ena_url = f"https://www.ebi.ac.uk/ena/portal/api/filereport?accession={accession}&result=read_run&fields=run_accession,fastq_ftp&format=tsv"
    
    try:
        # Facem request pentru a obține link-ul FASTQ
        with urllib.request.urlopen(ena_url) as response:
            data = response.read().decode('utf-8')
            lines = data.strip().split('\n')
            
            if len(lines) < 2:
                print(f"Eroare: Accession '{accession}' nu a fost găsit.")
                sys.exit(1)
            
            # Parsam răspunsul
            fastq_links = lines[1].split('\t')[1]  # Al doilea câmp conține link-urile
            
            if not fastq_links:
                print(f"Eroare: Nu s-a găsit FASTQ pentru '{accession}'.")
                sys.exit(1)
            
            # Preluam primul link (poate fi mai multe separate prin ';')
            fastq_url = f"ftp://{fastq_links.split(';')[0]}"
            
            # TODO: descărcați fișierul în Locația ALEASĂ DE VOI
            # Creem directorul de destinație
            dest_dir = Path(f"../../../../data/work/AlexTGoCreative/lab03")
            dest_dir.mkdir(parents=True, exist_ok=True)
            
            # Extragem numele fișierului din URL
            filename = fastq_url.split('/')[-1]
            filepath = dest_dir / filename
            
            # Descărcam fișierul
            print(f"Descărcam {accession}...", file=sys.stderr)
            urllib.request.urlretrieve(fastq_url, filepath)
            
            # TODO: print("Downloaded:", <cale_fisier>)
            print(f"Downloaded: {filepath}")
            
    except Exception as e:
        print(f"Eroare la descărcare: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
