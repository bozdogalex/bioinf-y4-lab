import os
import requests

handle = "MadalinaNegru"
outdir = f"data/work/{handle}/lab03"
outfile = os.path.join(outdir, "your_reads.fastq.gz")

os.makedirs(outdir, exist_ok=True)

# URL direct către FASTQ
url = "https://www.ebi.ac.uk/ena/browser/view/SRR17944414"

print(f"Descărcare în {outfile} ...")
r = requests.get(url, stream=True)
if r.status_code == 200:
    with open(outfile, "wb") as f:
        for chunk in r.iter_content(chunk_size=8192):
            f.write(chunk)
    print("Descărcare completă!")
else:
    print("Eroare la descărcare:", r.status_code)

