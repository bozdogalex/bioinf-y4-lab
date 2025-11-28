#!/usr/bin/env python3
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
import requests
import argparse
import time


def get_fastq_links(accession):
    """Obține link-urile FASTQ de la ENA API pentru un accession dat"""
    url = "https://www.ebi.ac.uk/ena/portal/api/filereport"
    params = {
        'accession': accession,
        'result': 'read_run',
        'fields': 'fastq_ftp,sample_accession,instrument_platform',
        'format': 'json',
        'download': 'false'
    }
    
    headers = {
        'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36'
    }
    
    try:
        print(f"🔍 Interoghez ENA API pentru {accession}...")
        response = requests.get(url, params=params, headers=headers, timeout=30)
        
        # Afișează detalii pentru debugging
        print(f"📡 Status code: {response.status_code}")
        print(f"📡 URL: {response.url}")
        
        if response.status_code != 200:
            raise RuntimeError(f"EROARE HTTP {response.status_code}: {response.text}")
        
        # Verifică conținutul răspunsului
        if not response.text.strip():
            raise ValueError("Răspuns gol de la server")
        
        # Încearcă să parseze JSON
        try:
            data = response.json()
        except ValueError as e:
            print(f"❌ Răspuns primit (primele 500 caractere): {response.text[:500]}")
            raise ValueError(f"Nu pot parsa răspunsul JSON: {e}")
        
        if not data:
            raise ValueError(f"Nu s-au găsit date pentru accession: {accession}")
        
        # Extrage link-urile FTP din răspuns
        ftp_links = data[0].get('fastq_ftp', '')
        if not ftp_links:
            raise ValueError(f"Nu s-au găsit link-uri FASTQ pentru: {accession}")
        
        # Afișează informații despre sample
        sample_acc = data[0].get('sample_accession', 'N/A')
        platform = data[0].get('instrument_platform', 'N/A')
        print(f"📊 Sample: {sample_acc}, Platform: {platform}")
        
        # Transformă în link-uri HTTP (ENA permite acest lucru)
        http_links = []
        for ftp_link in ftp_links.split(';'):
            if ftp_link.strip():
                # Folosește FTP direct sau transformă în HTTP
                http_link = f"https://{ftp_link.strip()}"
                http_links.append(http_link)
        
        print(f"🔗 Am găsit {len(http_links)} fișier(e) FASTQ")
        for link in http_links:
            print(f"   - {os.path.basename(link)}")
        
        return http_links
    
    except requests.exceptions.Timeout:
        raise RuntimeError("Timeout la interogarea ENA API")
    except requests.exceptions.ConnectionError:
        raise RuntimeError("Eroare de conexiune la ENA API")
    except requests.exceptions.RequestException as e:
        raise RuntimeError(f"Eroare la interogarea ENA API: {e}")


def download_file(url, output_path):
    """Descarcă un fișier de la URL la calea specificată"""
    try:
        print(f"⬇️  Descarc {os.path.basename(url)}...")
        
        headers = {
            'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36'
        }
        
        response = requests.get(url, stream=True, headers=headers, timeout=60)
        response.raise_for_status()
        
        total_size = int(response.headers.get('content-length', 0))
        downloaded_size = 0
        
        with open(output_path, 'wb') as f:
            for chunk in response.iter_content(chunk_size=8192):
                if chunk:
                    f.write(chunk)
                    downloaded_size += len(chunk)
                    
                    # Afișează progres pentru fișiere mari
                    if total_size > 0:
                        percent = (downloaded_size / total_size) * 100
                        print(f"   📦 Progres: {percent:.1f}%", end='\r')
        
        print()  # Linie nouă după progres
        return output_path
    
    except requests.exceptions.RequestException as e:
        # Șterge fișierul parțial descărcat
        if os.path.exists(output_path):
            os.remove(output_path)
        raise RuntimeError(f"Eroare la descărcare: {e}")


def main():
    parser = argparse.ArgumentParser(description='Descarcă fișier FASTQ pentru un accession')
    parser.add_argument('accession', help='Accession-ul (ex: SRR..., ERR...)')
    parser.add_argument('--output', '-o', help='Calea de output (opțional)')
    
    args = parser.parse_args()
    
    # Calea de output implicită
    if args.output:
        output_path = args.output
    else:
        # Creați directorul dacă nu există
        output_dir = f"/workspaces/bioinf-y4-lab/labs/01_intro&databases/data/work/Botoaca-Florentina-Veronica/lab01"
        os.makedirs(output_dir, exist_ok=True)
        output_path = os.path.join(output_dir, "your_reads.fastq.gz")
    
    try:
        # Obține link-urile FASTQ
        fastq_links = get_fastq_links(args.accession)
        
        if not fastq_links:
            print("❌ Nu s-au găsit link-uri FASTQ")
            sys.exit(1)
        
        # Descarcă primul fișier FASTQ (cerința minimă)
        first_link = fastq_links[0]
        downloaded_file = download_file(first_link, output_path)
        
        # Verifică dimensiunea fișierului
        file_size = os.path.getsize(downloaded_file)
        file_size_mb = file_size / (1024 * 1024)
        
        print(f"✅ Descărcat cu succes: {downloaded_file}")
        print(f"📊 Dimensiune: {file_size_mb:.2f} MB ({file_size} bytes)")
        
    except Exception as e:
        print(f"❌ Eroare: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()