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

"""

Exemplu folosire:
    python ex03_download_fastq.py SRR000001
    python ex03_download_fastq.py ERR000001 --output data/work/olezhka0809/lab03/
"""

import sys
import os
import requests
from pathlib import Path


def get_fastq_url(accession):
    """
    Obține URL-ul pentru download FASTQ de la ENA (European Nucleotide Archive).
    ENA oferă un API simplu pentru a obține link-uri directe.
    """
    # ENA API pentru metadata
    ena_api = f"https://www.ebi.ac.uk/ena/portal/api/filereport"
    params = {
        "accession": accession,
        "result": "read_run",
        "fields": "run_accession,fastq_ftp,fastq_aspera,fastq_galaxy",
        "format": "json",
        "download": "true"
    }
    
    try:
        response = requests.get(ena_api, params=params, timeout=10)
        response.raise_for_status()
        data = response.json()
        
        if not data:
            print(f"Eroare: Accession '{accession}' nu a fost găsit în ENA.")
            return None
        
        # Extragem primul link FTP disponibil
        fastq_ftp = data[0].get("fastq_ftp", "")
        if not fastq_ftp:
            print(f"Eroare: Nu există fișiere FASTQ pentru {accession}")
            return None
        
        # fastq_ftp poate conține multiple fișiere separate prin ';'
        # Luăm primul fișier (pentru paired-end ar fi 2 fișiere)
        first_file = fastq_ftp.split(";")[0]
        
        # Convertim FTP la HTTP pentru download mai simplu
        url = f"http://{first_file}"
        return url
        
    except requests.exceptions.RequestException as e:
        print(f"Eroare la interogarea ENA: {e}")
        return None


def download_fastq(url, output_path):
    """
    Descarcă fișierul FASTQ de la URL specificat.
    """
    try:
        print(f"Descărcare din: {url}")
        print("Acest proces poate dura câteva minute...")
        
        # Download cu progress
        response = requests.get(url, stream=True, timeout=30)
        response.raise_for_status()
        
        total_size = int(response.headers.get('content-length', 0))
        downloaded = 0
        
        with open(output_path, 'wb') as f:
            for chunk in response.iter_content(chunk_size=8192):
                if chunk:
                    f.write(chunk)
                    downloaded += len(chunk)
                    # Progress simplu
                    if total_size > 0:
                        percent = (downloaded / total_size) * 100
                        print(f"\rProgress: {percent:.1f}%", end="", flush=True)
        
        print()  # newline după progress
        return True
        
    except requests.exceptions.RequestException as e:
        print(f"Eroare la download: {e}")
        return False


def main():
    # Parsare argumente
    if len(sys.argv) < 2:
        print("Folosire: python ex03_download_fastq.py <ACCESSION> [--output <DIR>]")
        print("Exemplu: python ex03_download_fastq.py SRR000001")
        print("\nAccessions TP53-related sugerate pentru testare:")
        print("  - SRR000001 (mic, pentru testare rapidă)")
        print("  - ERR000001")
        sys.exit(1)
    
    accession = sys.argv[1]
    
    # Directorul de output implicit
    output_dir = "data/work/olezhka0809/lab03"
    
    # Verifică dacă utilizatorul a specificat alt director
    if "--output" in sys.argv:
        idx = sys.argv.index("--output")
        if idx + 1 < len(sys.argv):
            output_dir = sys.argv[idx + 1]
    
    # Creează directorul dacă nu există
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    
    # Obține URL-ul pentru download
    print(f"Căutare accession: {accession}")
    url = get_fastq_url(accession)
    
    if not url:
        print("Nu s-a putut obține URL-ul pentru download.")
        sys.exit(1)
    
    # Determină numele fișierului
    filename = url.split("/")[-1]
    output_path = os.path.join(output_dir, filename)
    
    # Verifică dacă fișierul există deja
    if os.path.exists(output_path):
        print(f"Fișierul există deja: {output_path}")
        response = input("Doriți să-l descărcați din nou? (y/n): ")
        if response.lower() != 'y':
            print(f"Downloaded: {output_path}")
            return
    
    # Descarcă fișierul
    success = download_fastq(url, output_path)
    
    if success:
        # Verifică dimensiunea fișierului descărcat
        file_size = os.path.getsize(output_path)
        file_size_mb = file_size / (1024 * 1024)
        
        print(f"\nDescărcare completă!")
        print(f"  Dimensiune: {file_size_mb:.2f} MB")
        print(f"Downloaded: {output_path}")
    else:
        print("Descărcarea a eșuat.")
        sys.exit(1)


if __name__ == "__main__":
    main()
