"""
Exercițiu 03 — Descărcare FASTQ (student-owned)

Obiectiv:
- Alegeți un accession TP53-related (ex. SRR..., ERR...) și DESCĂRCAȚI un fișier FASTQ.
- Salvați in  data/work/<handle>/lab03/your_reads.fastq.gz

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

BASE_DIR = "data/work/filippaulescu/lab03"

def get_fastq_url(accession):
    """
    Interogheaza ENA (European Nucleotide Archive) pentru a obtine URL-ul FASTQ.
    ENA ofera un API simplu de cautare care returneaza metadate in format TSV.
    """
    # Campurile necesare: "fastq_ftp" (link-ul de descarcare)
    ENA_URL = (
        f"https://www.ebi.ac.uk/ena/portal/api/filereport?accession={accession}&"
        f"result=read_run&fields=fastq_ftp&download=false"
    )
    
    try:
        with urllib.request.urlopen(ENA_URL) as response:
            data = response.read().decode('utf-8').strip()
    except urllib.error.URLError as e:
        print(f"Eroare la interogarea ENA: {e}", file=sys.stderr)
        return None

    # ENA returneaza un header (fastq_ftp) si apoi link-ul pe linia urmatoare
    lines = data.split('\n')
    if len(lines) < 2:
        print(f"Eroare: Nu s-au gasit link-uri FASTQ pentru accession-ul {accession}.", file=sys.stderr)
        return None

    # CORECTIE PENTRU TABULATIE: Separam pe tabulator (TSV) si apoi pe punct si virgula (paired-end)
    # Ia ultimul câmp (link-ul FTP) și apoi îl separă pe ';'
    fastq_field = lines[1].split('\t')[-1] 
    ftp_links = fastq_field.split(';')

    # Pentru cerinta 'un singur fisier e suficient', alegem primul link.
    if ftp_links and ftp_links[0]:
        # 1. Ne asiguram ca e un link complet, daca nu, ii adaugam prefixul 'ftp.sra.ebi.ac.uk/'
        full_ftp_link = ftp_links[0]
        if not full_ftp_link.startswith('ftp'):
             full_ftp_link = f"ftp.sra.ebi.ac.uk/{full_ftp_link}"
             
        # 2. Inlocuim protocolul ftp:// cu https://
        http_link = full_ftp_link.replace("ftp://", "https://")
        
        # 3. Ne asiguram ca protocolul este prezent, in cazul in care full_ftp_link nu avea 'ftp://'
        # Aceasta este corectia cheie: adauga "https://" daca link-ul inca nu incepe cu el
        if not http_link.startswith('https://'):
            http_link = f"https://{http_link}"
            
        # Inlocuim 'ftp.sra.ebi.ac.uk' cu 'ftp.ebi.ac.uk' (gazda mai stabila pentru HTTPS)
        http_link = http_link.replace("ftp.sra.ebi.ac.uk", "ftp.ebi.ac.uk")
        
        return http_link
    
    return None
def main():
    # TODO: citiți accession-ul (ex. sys.argv)
    # TODO: interogați sursa (ENA/SRA) pentru link FASTQ
    # TODO: descărcați fișierul în Locația ALEASĂ DE VOI
    # TODO: print("Downloaded:", <cale_fisier>)
    if len(sys.argv) < 2:
        print("Utilizare: python script.py <accession_SRR_ERR>", file=sys.stderr)
        sys.exit(1)

    accession = sys.argv[1]

    # 1. Interogare ENA pentru URL-ul FASTQ
    fastq_url = get_fastq_url(accession)
    
    if not fastq_url:
        print(f"Nu s-a putut obtine URL-ul pentru {accession}.", file=sys.stderr)
        sys.exit(1)

    # 2. Setare cale de salvare (data/work/<handle>/lab03/your_reads.fastq.gz)
    os.makedirs(BASE_DIR, exist_ok=True)
    
    # Numele fisierului de pe server (ex: SRR..._1.fastq.gz)
    remote_filename = fastq_url.split('/')[-1]
    
    # Numele final al fisierului conform cerintei
    local_path = os.path.join(BASE_DIR, "your_reads.fastq.gz")
    
    # 3. Descarcare fisier
    print(f"Incepere descarcare {remote_filename} de la {fastq_url}...", file=sys.stderr)
    
    try:
        # Folosim urlretrieve pentru descarcare directa intr-un fisier
        urllib.request.urlretrieve(fastq_url, local_path)
        
        # 4. Afisare calea fisierului descarcat
        print("Downloaded:", os.path.abspath(local_path))
        
    except urllib.error.URLError as e:
        print(f"Eroare la descarcarea fisierului: {e}", file=sys.stderr)
        # Sterge fisierul incomplet, daca exista
        if os.path.exists(local_path):
            os.remove(local_path)
        sys.exit(1)
    pass


if __name__ == "__main__":
    main()