#!/usr/bin/env python3
"""
VCF to PubMed Query Script
Parses a VCF file, extracts variants, and queries PubMed for related literature.
"""

from Bio import Entrez
import sys
from pathlib import Path

# Configure Entrez
Entrez.email = "alexandru.tulbure@student.upt.ro"  
HANDLE = "AlexTGoCreative"


def parse_vcf(vcf_file):
    """
    Parse a VCF file and extract variant information.
    
    Args:
        vcf_file (str): Path to the VCF file
        
    Returns:
        list: List of variant dictionaries
    """
    variants = []
    
    try:
        with open(vcf_file, 'r') as f:
            for line in f:
                line = line.strip()
                
                # Skip empty lines and header lines
                if not line or line.startswith('##'):
                    continue
                
                # Column header line
                if line.startswith('#CHROM'):
                    continue
                
                # Parse variant line
                # Split by whitespace (handles both tabs and multiple spaces)
                fields = line.split()
                if len(fields) < 5:
                    continue
                
                variant = {
                    'chrom': fields[0],
                    'pos': fields[1],
                    'id': fields[2],
                    'ref': fields[3],
                    'alt': fields[4],
                    'qual': fields[5] if len(fields) > 5 else '.',
                    'filter': fields[6] if len(fields) > 6 else '.',
                    'info': fields[7] if len(fields) > 7 else '.'
                }
                
                variants.append(variant)
        
        print(f"Parsed {len(variants)} variants from VCF file")
        
    except Exception as e:
        print(f"Error parsing VCF file: {e}")
        sys.exit(1)
    
    return variants


def search_pubmed_variant(query, max_results=5):
    """
    Search PubMed for a specific query and return article count and IDs.
    
    Args:
        query (str): Search query string
        max_results (int): Maximum number of results to return
        
    Returns:
        tuple: (count, pmid_list)
    """
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results)
        record = Entrez.read(handle)
        handle.close()
        count = int(record.get("Count", 0))
        pmid_list = record.get("IdList", [])
        return count, pmid_list
    except Exception as e:
        print(f"Error searching PubMed for '{query}': {e}")
        return 0, []


def fetch_article_summary(pmid_list):
    """
    Fetch article summaries for a list of PubMed IDs.
    
    Args:
        pmid_list (list): List of PubMed IDs
        
    Returns:
        list: List of article summary dictionaries
    """
    if not pmid_list:
        return []
    
    summaries = []
    
    try:
        # Fetch article summaries
        handle = Entrez.esummary(db="pubmed", id=",".join(pmid_list))
        records = Entrez.read(handle)
        handle.close()
        
        for record in records:
            if isinstance(record, dict):
                summary = {
                    'pmid': record.get('Id', 'N/A'),
                    'title': record.get('Title', 'Title not available'),
                    'authors': record.get('AuthorList', []),
                    'source': record.get('Source', 'N/A'),
                    'pubdate': record.get('PubDate', 'N/A')
                }
                summaries.append(summary)
        
    except Exception as e:
        print(f"Error fetching article summaries: {e}")
    
    return summaries


def query_variant(variant):
    """
    Query PubMed for a specific variant.
    
    Args:
        variant (dict): Variant dictionary
        
    Returns:
        dict: Query results with count and articles
    """
    result = {
        'variant': variant,
        'query': '',
        'count': 0,
        'articles': []
    }
    
    # Determine query based on variant ID (primary strategy)
    if variant['id'] != '.' and variant['id'].startswith('rs'):
        # Use rsID if available
        query = variant['id']
        result['query'] = query
        print(f"\nQuerying PubMed for rsID: {query}")
    else:
        # Use chromosomal position and TP53
        query = f"chr{variant['chrom']}:{variant['pos']} AND TP53"
        result['query'] = query
        print(f"\nQuerying PubMed for position: {query}")
    
    # Search PubMed
    count, pmid_list = search_pubmed_variant(query, max_results=5)
    result['count'] = count
    
    print(f"  Found {count} articles")
    
    # If no results, try fallback with mutation name from INFO field
    if count == 0 and 'NOTE=' in variant['info']:
        note_value = variant['info'].split('NOTE=')[1].split(';')[0]
        # Extract mutation notation (e.g., TP53_R273H -> R273H)
        if '_' in note_value:
            mutation_name = note_value.split('_')[1]
            query = f"TP53 {mutation_name}"
            result['query'] = f"{result['query']} (fallback: {query})"
            print(f"  No results. Trying fallback query: {query}")
            
            count, pmid_list = search_pubmed_variant(query, max_results=5)
            result['count'] = count
            print(f"  Found {count} articles with fallback query")
    
    # Fetch article summaries if results found
    if pmid_list:
        print(f"  Fetching details for {len(pmid_list)} articles...")
        result['articles'] = fetch_article_summary(pmid_list)
    
    return result


def save_results(results, output_file):
    """
    Save variant query results to a text file.
    
    Args:
        results (list): List of result dictionaries
        output_file (str): Output filename
    """
    try:
        with open(output_file, 'w', encoding='utf-8') as f:
            f.write("=" * 80 + "\n")
            f.write("VCF Variant → PubMed Query Results\n")
            f.write("=" * 80 + "\n\n")
            
            for i, result in enumerate(results, 1):
                variant = result['variant']
                
                f.write(f"Variant {i}\n")
                f.write("-" * 80 + "\n")
                f.write(f"Chromosome: {variant['chrom']}\n")
                f.write(f"Position: {variant['pos']}\n")
                f.write(f"ID: {variant['id']}\n")
                f.write(f"Reference: {variant['ref']}\n")
                f.write(f"Alternate: {variant['alt']}\n")
                f.write(f"Quality: {variant['qual']}\n")
                f.write(f"Filter: {variant['filter']}\n\n")
                
                f.write(f"PubMed Query: {result['query']}\n")
                f.write(f"Total Articles Found: {result['count']}\n\n")
                
                if result['articles']:
                    f.write("Top Articles:\n")
                    f.write("-" * 80 + "\n")
                    
                    for j, article in enumerate(result['articles'], 1):
                        f.write(f"\n{j}. PMID: {article['pmid']}\n")
                        f.write(f"   Title: {article['title']}\n")
                        
                        # Format authors
                        if isinstance(article['authors'], list) and article['authors']:
                            authors_str = ", ".join(article['authors'][:5])
                            if len(article['authors']) > 5:
                                authors_str += " et al."
                            f.write(f"   Authors: {authors_str}\n")
                        
                        f.write(f"   Source: {article['source']}\n")
                        f.write(f"   Date: {article['pubdate']}\n")
                else:
                    f.write("No articles found for this variant.\n")
                
                f.write("\n" + "=" * 80 + "\n\n")
        
        print(f"\nResults saved to {output_file}")
        
    except Exception as e:
        print(f"Error saving results: {e}")


def main():
    """Main function to execute VCF to PubMed query."""
    # Define file paths
    vcf_file = "test_AlexTGoCreative.vcf"
    output_file = f"variants_{HANDLE}.txt"
    
    # Check if VCF file exists
    if not Path(vcf_file).exists():
        print(f"Error: VCF file not found at {vcf_file}")
        sys.exit(1)
    
    print("=" * 80)
    print("VCF to PubMed Query Analysis")
    print("=" * 80)
    print(f"Input VCF: {vcf_file}")
    print(f"Output file: {output_file}")
    print("=" * 80)
    
    # Parse VCF file
    print("\nParsing VCF file...")
    variants = parse_vcf(vcf_file)
    
    if not variants:
        print("No variants found in VCF file.")
        sys.exit(1)
    
    # Select variants to query (at least 2, up to 5)
    num_variants = min(5, len(variants))
    selected_variants = variants[:num_variants]
    
    print(f"\nSelected {len(selected_variants)} variants for PubMed query")
    
    # Query PubMed for each variant
    results = []
    for variant in selected_variants:
        result = query_variant(variant)
        results.append(result)
    
    # Save results
    print("\nSaving results...")
    save_results(results, output_file)
    
    print("\n" + "=" * 80)
    print("Analysis complete!")
    print("=" * 80)


if __name__ == "__main__":
    main()
