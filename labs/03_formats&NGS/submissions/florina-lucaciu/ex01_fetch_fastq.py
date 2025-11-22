#!/usr/bin/env python3
"""
Exercițiu 03 — Descărcare FASTQ (student-owned)
Scriptul acceptă un accession (ex. ERR10189843), interoghează ENA
pentru link-ul FASTQ și îl descarcă în locația specificată.
"""

#pentru rulare: python labs/'03_formats&NGS'/submissions/florina-lucaciu/ex01_fetch_fastq.py SRR32749910

import sys
import os
import urllib.request
import urllib.error


USER_HANDLE = "florina-lucaciu"



def get_ena_fastq_url(accession: str) -> str:
    """
    Interoghează API-ul ENA pentru a obține link-ul de descărcare FTP
    pentru un accession (Run) dat.
    """
    print(f"Se interoghează ENA API pentru {accession}...")
    
    api_url = (
        f"https://www.ebi.ac.uk/ena/portal/api/filereport"
        f"?accession={accession}"
        f"&result=read_run"
        f"&fields=fastq_ftp"
        f"&format=tsv"
    )

    try:
        with urllib.request.urlopen(api_url) as response:
            data = response.read().decode('utf-8').strip()
            
            lines = data.split('\n')
            
            if len(lines) < 2:
                raise ValueError(f"Nu s-au gasit date valide pentru {accession}. Raspuns API: {data}")
            
            header = lines[0].split('\t')
            data_line = lines[1].split('\t')
            
            try:
                ftp_index = header.index('fastq_ftp')
            except ValueError:
                raise ValueError("Raspunsul API ENA nu contine coloana 'fastq_ftp'.")

            ftp_paths_str = data_line[ftp_index]
            
            if not ftp_paths_str:
                raise ValueError(f"Niciun fisier FASTQ (FTP) nu a fost găsit pentru {accession}.")

            first_ftp_path = ftp_paths_str.split(';')[0]
            

            # Calea din API (first_ftp_path) de forma 'ftp.sra.ebi.ac.uk/vol1/...'
            full_ftp_url = f"ftp://{first_ftp_path}"
            
            return full_ftp_url

    except urllib.error.URLError as e:
        raise ConnectionError(f"Eroare de retea la contactarea ENA API: {e.reason}")
    except Exception as e:
        raise RuntimeError(f"Eroare la procesarea raspunsului ENA: {e}")

def download_file(url: str, save_path: str):
    """
    Descarca un fisier de la un URL si il salveaza in calea specificata.
    """
    print(f"Incepe descarcarea de la: {url}")
    print(f"Se salveaza in: {save_path}")
    
    try:

        output_dir = os.path.dirname(save_path)
        if output_dir: 
            os.makedirs(output_dir, exist_ok=True)
        
        # Descarcă fișierul
        urllib.request.urlretrieve(url, save_path)
        
        print("\nDescarcare finalizată cu succes.")
        
    except urllib.error.URLError as e:
        raise ConnectionError(f"Eroare la descarcarea fisierului: {e.reason}")
    except OSError as e:
        raise OSError(f"Eroare la scrierea fisierului pe disc: {e}")
    except Exception as e:
        raise RuntimeError(f"O eroare necunoscuta a aparut la descarcare: {e}")


def main():
    # TODO: citiți accession-ul (ex. sys.argv)
    if len(sys.argv) < 2:
        print(f"Eroare: Trebuie sa furnizati un accession (ex. SRR... sau ERR...).", file=sys.stderr)
        print(f"Utilizare: python3 {sys.argv[0]} <ACCESSION>", file=sys.stderr)
        sys.exit(1)
        
    accession = sys.argv[1]

    
    if USER_HANDLE == "<handle>":
        print("EROARE: Editati scriptul si modificati variabila USER_HANDLE.", file=sys.stderr)
        sys.exit(1)

    try:
        # TODO: interogați sursa (ENA/SRA) pentru link FASTQ
        fastq_url = get_ena_fastq_url(accession)
        
        
        output_dir = os.path.join("data", "work", USER_HANDLE, "lab03")
        output_filename = "your_reads.fastq.gz" 
        output_path = os.path.join(output_dir, output_filename)

        # TODO: descărcați fișierul
        download_file(fastq_url, output_path)

        # TODO: print("Downloaded:", <cale_fisier>)
        print(f"Downloaded: {output_path}")

    except Exception as e:
        print(f"\nA aparut o eroare: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()