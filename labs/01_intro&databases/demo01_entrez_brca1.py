"""Inspect a versioned BRCA1 transcript offline; optionally refresh from NCBI."""
import argparse
from io import StringIO
import os
from pathlib import Path
import socket

from Bio import Entrez, SeqIO
from Bio.SeqUtils import gc_fraction

ACCESSION = "NM_007294.4"
CACHE = Path(__file__).resolve().parent / "data" / "brca1.gb"


def validate_record(record):
    """Reject chromosome/contig records and incomplete sequence data."""
    if (record.id != ACCESSION or not record.seq.defined
            or not len(record) or "BRCA1" not in record.description
            or record.annotations.get("organism") != "Homo sapiens"
            or record.annotations.get("molecule_type") != "mRNA"):
        raise ValueError(f"Expected a complete human BRCA1 transcript {ACCESSION}")
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--refresh", action="store_true", help="Download the fixed accession")
    parser.add_argument("--email", default=os.environ.get("NCBI_EMAIL"))
    parser.add_argument("--out", type=Path, default=Path("data/work/demo/lab01/brca1.gb"))
    args = parser.parse_args()
    if args.refresh:
        if not args.email:
            parser.error("Set NCBI_EMAIL or --email to your real email; omit --refresh for offline use.")
        Entrez.email = args.email
        Entrez.tool = "bioinf_y4_lab"
        Entrez.api_key = os.environ.get("NCBI_API_KEY")
        Entrez.max_tries = 2
        Entrez.sleep_between_tries = 1
        socket.setdefaulttimeout(20)
        try:
            with Entrez.efetch(db="nuccore", id=ACCESSION, rettype="gb", retmode="text") as handle:
                text = handle.read()
            record = validate_record(SeqIO.read(StringIO(text), "genbank"))
        except (OSError, ValueError) as exc:
            raise SystemExit(f"NCBI retrieval failed: {exc}. Retry without --refresh for cached data.")
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(text, encoding="utf-8")
        source = str(args.out)
    else:
        record = validate_record(SeqIO.read(CACHE, "genbank"))
        source = f"bundled offline cache: {CACHE.name}"
    print("Source:", source)
    print("ID:", record.id)
    print("Title:", record.description)
    print("Length:", len(record), "nt")
    print("GC fraction:", round(gc_fraction(record.seq), 3))
    print("First 50 nt:", record.seq[:50])
    print("This is an mRNA transcript represented using T, not the complete genomic locus.")


if __name__ == "__main__":
    main()
