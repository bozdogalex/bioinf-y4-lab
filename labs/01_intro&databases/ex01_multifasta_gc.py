#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Exercise (Lab 1): FASTA Download + GC Content Calculation

Purpose:
  1) Download a FASTA file from NCBI (nucleotide or protein).
  2) Save the file locally in data/work/<handle>/lab01/ (DO NOT upload it to git).
  3) Compute the GC fraction for each record in the file.

Instructions:
  - Run the script with the required arguments (examples):
      python ex01_multifasta_gc.py --email student@example.com \
        --query "TP53[Gene] AND Homo sapiens[Organism]" \
        --retmax 3 \
        --out data/work/<handle>/lab01/my_tp53.fa

      python ex01_multifasta_gc.py --email student@example.com \
        --accession NM_000546 \
        --out data/work/<handle>/lab01/nm000546.fa

  - Steps to complete:
    1) Configure Entrez with email (and optional api_key).
    2) If you receive an accession → download that record using efetch.
    3) If you receive a query → perform esearch for IdList, then efetch those IDs.
    4) Write the results to the file specified by --out.
    5) Read the local FASTA file and calculate GC for each sequence.
    6) Display results on screen: <id>\tGC=<value with 3 decimals>.
"""

import argparse
from pathlib import Path
import sys

from Bio import SeqIO
# from Bio import Entrez  # TODO: uncomment and use for download


def gc_fraction(seq: str) -> float:
    """GC fraction for a sequence; robust to lower/upper case and non-ATGC."""
    s = seq.upper()
    atgc = [c for c in s if c in ("A", "T", "G", "C")]
    if not atgc:
        return 0.0
    g = atgc.count("G")
    c = atgc.count("C")
    return (g + c) / float(len(atgc))


def download_fasta(email: str, out_path: Path, query: str = None,
                   accession: str = None, db: str = "nuccore",
                   retmax: int = 3, api_key: str = None) -> int:
    """
    TODO: Implement downloading from NCBI.
    Steps:
      - Configure Entrez with email (and optional api_key).
      - If accession is provided: download that record.
      - Otherwise, if query is provided: perform esearch -> list of IDs, then efetch.
      - Write results to out_path.
      - Return the number of records written.
    """
    raise NotImplementedError("TODO: implement download using Entrez")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--email", required=True, help="Required email for NCBI Entrez")
    ap.add_argument("--api_key", help="NCBI API key (optional)")
    ap.add_argument("--query", help='Ex: "TP53[Gene] AND Homo sapiens[Organism]"')
    ap.add_argument("--accession", help="Ex: NM_000546")
    ap.add_argument("--db", default="nuccore", choices=["nuccore", "protein"])
    ap.add_argument("--retmax", type=int, default=3)
    ap.add_argument("--out", required=True, help="Output FASTA file")
    args = ap.parse_args()

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    # TODO: Call the download_fasta(...) function and save the results
    # n = download_fasta(args.email, out_path, query=args.query,
    #                    accession=args.accession, db=args.db,
    #                    retmax=args.retmax, api_key=args.api_key)
    # print(f"[ok] Wrote {n} records to: {out_path}")

    # TODO: Read the FASTA file with SeqIO.parse
    # records = ...

    # TODO: Compute GC for each sequence and print the results
    # for rec in records:
    #     print(f"{rec.id}\tGC={value:.3f}")


if __name__ == "__main__":
    main()