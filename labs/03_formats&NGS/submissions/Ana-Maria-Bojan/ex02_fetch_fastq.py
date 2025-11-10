
import sys
import os
import urllib.request

def main():
    # TODO: citiți accession-ul (ex. sys.argv)
    # TODO: interogați sursa (ENA/SRA) pentru link FASTQ
    # TODO: descărcați fișierul în Locația ALEASĂ DE VOI
    # TODO: print("Downloaded:", <cale_fisier>)
        
    accession = sys.argv[1]   

    out_dir = f"/workspaces/bioinf-y4-lab/data/work/Ana-Maria-Bojan/lab03"
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, "your_reads.fastq.gz")

    api_url = f"https://www.ebi.ac.uk/ena/portal/api/filereport?accession={accession}&result=read_run&fields=fastq_ftp&format=tsv"

    with urllib.request.urlopen(api_url) as response:
        text = response.read().decode("utf-8")

    lines = text.strip().split("\n")
    if len(lines) < 2:
        print("nu s-au gasit linkuri FASTQ ")
        sys.exit(1)

    ftp_link = lines[1].split("\t")[-1].split(";")[0]  
    https_link = "https://" + ftp_link  

    print(f"Downloading {https_link} ...")
    urllib.request.urlretrieve(https_link, out_path)

    print(f"Downloaded: {out_path}")


if __name__ == "__main__":
    main()