#!/usr/bin/env python3
"""
Exercise 01: Multi-FASTA GC content calculator with Entrez download
Complete the TODO sections to fetch sequences and calculate GC content.
"""

import argparse
from Bio import Entrez, SeqIO
from pathlib import Path


def fetch_by_query(email: str, query: str, retmax: int, output_file: str):
    """Fetch sequences using an Entrez query and save to FASTA."""
    Entrez.email = email
    
    # TODO: Search nucleotide database
    print(f"Searching for: {query}")
    handle = Entrez.esearch(db="nucleotide", term=query, retmax=retmax)
    record = Entrez.read(handle)
    handle.close()
    
    id_list = record["IdList"]
    print(f"Found {len(id_list)} sequences")
    
    if not id_list:
        print("No sequences found!")
        return
    
    # TODO: Fetch sequences in FASTA format
    print(f"Fetching {len(id_list)} sequences...")
    handle = Entrez.efetch(db="nucleotide", id=id_list, rettype="fasta", retmode="text")
    sequences = handle.read()
    handle.close()
    
    # Save to file
    Path(output_file).parent.mkdir(parents=True, exist_ok=True)
    with open(output_file, 'w') as f:
        f.write(sequences)
    print(f"Saved to {output_file}")


def fetch_by_accession(email: str, accession: str, output_file: str):
    """Fetch a single sequence by accession number."""
    Entrez.email = email
    
    # TODO: Fetch by accession
    print(f"Fetching accession: {accession}")
    handle = Entrez.efetch(db="nucleotide", id=accession, rettype="fasta", retmode="text")
    sequence = handle.read()
    handle.close()
    
    # Save to file
    Path(output_file).parent.mkdir(parents=True, exist_ok=True)
    with open(output_file, 'w') as f:
        f.write(sequence)
    print(f"Saved to {output_file}")


def calculate_gc_content(fasta_file: str):
    """Calculate GC content for all sequences in a FASTA file."""
    print(f"\nCalculating GC content for {fasta_file}")
    print("-" * 60)
    
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq = str(record.seq).upper()
        
        # TODO: Calculate GC content
        gc_count = seq.count('G') + seq.count('C')
        total = len(seq)
        gc_fraction = gc_count / total if total > 0 else 0
        
        print(f"ID: {record.id}")
        print(f"Description: {record.description}")
        print(f"Length: {total} bp")
        print(f"GC content: {gc_fraction:.4f} ({gc_fraction*100:.2f}%)")
        print("-" * 60)


def main():
    parser = argparse.ArgumentParser(description="Fetch sequences and calculate GC content")
    parser.add_argument("--email", required=True, help="Your email for Entrez")
    parser.add_argument("--query", help="Entrez search query")
    parser.add_argument("--accession", help="Single accession number")
    parser.add_argument("--retmax", type=int, default=5, help="Max results for query")
    parser.add_argument("--out", required=True, help="Output FASTA file")
    
    args = parser.parse_args()
    
    # Fetch sequences
    if args.query:
        fetch_by_query(args.email, args.query, args.retmax, args.out)
    elif args.accession:
        fetch_by_accession(args.email, args.accession, args.out)
    else:
        print("ERROR: Provide either --query or --accession")
        return
    
    # Calculate GC content
    calculate_gc_content(args.out)


if __name__ == "__main__":
    main()