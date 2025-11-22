import os
import requests

# 3 secvențe TP53 (Human, Mouse, Human Variant)
FASTA_ACCESSIONS = [
    "AF209148",   # Human TP53
    "AF151353",   # Mouse TP53
    "AF209133"    # Human TP53 variant
]

def download_fasta(accession):
    """Download FASTA from ENA."""
    url = f"https://www.ebi.ac.uk/ena/browser/api/fasta/{accession}?download=true"
    r = requests.get(url)

    if r.status_code != 200:
        raise RuntimeError(f"Could not download FASTA for accession {accession}")

    return r.text


def main():
    out_dir = os.path.expanduser("/workspaces/bioinf-y4-lab/data/work/LeachAntonia/lab04")
    os.makedirs(out_dir, exist_ok=True)

    out_path = os.path.join(out_dir, "tp53_sequences.fasta")

    print(f"Creating multi-FASTA -> {out_path}")

    with open(out_path, "w") as f:
        for acc in FASTA_ACCESSIONS:
            print(f"  - downloading {acc}")
            seq = download_fasta(acc)
            f.write(seq.strip() + "\n\n")

    print("\nMulti-FASTA created:", out_path)


if __name__ == "__main__":
    main()