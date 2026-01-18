"""
Fetch TP53 protein sequences from multiple species for phylogenetic analysis
"""

from Bio import Entrez, SeqIO
from pathlib import Path

# Configure email for Entrez
Entrez.email = "alexTGoCreative@example.com"

def fetch_tp53_sequences():
    """Fetch TP53 protein sequences from different species"""
    
    # Search for TP53 protein sequences from various species
    species_queries = [
        "TP53[Gene] AND Homo sapiens[Organism] AND p53[Protein]",
        "TP53[Gene] AND Mus musculus[Organism]",
        "TP53[Gene] AND Rattus norvegicus[Organism]",
        "TP53[Gene] AND Danio rerio[Organism]",
        "TP53[Gene] AND Xenopus tropicalis[Organism]",
        "TP53[Gene] AND Gallus gallus[Organism]",
        "TP53[Gene] AND Bos taurus[Organism]",
        "TP53[Gene] AND Canis familiaris[Organism]",
        "TP53[Gene] AND Sus scrofa[Organism]",
        "TP53[Gene] AND Macaca mulatta[Organism]",
    ]
    
    sequences = []
    
    for query in species_queries:
        print(f"Searching: {query}")
        try:
            # Search for protein IDs
            search_handle = Entrez.esearch(db="protein", term=query, retmax=1)
            search_results = Entrez.read(search_handle)
            search_handle.close()
            
            if search_results["IdList"]:
                protein_id = search_results["IdList"][0]
                
                # Fetch the sequence
                fetch_handle = Entrez.efetch(db="protein", id=protein_id, rettype="fasta", retmode="text")
                record = SeqIO.read(fetch_handle, "fasta")
                fetch_handle.close()
                
                sequences.append(record)
                print(f"  ✓ Fetched: {record.id} ({record.description[:60]}...)")
            else:
                print(f"  ✗ No sequences found")
                
        except Exception as e:
            print(f"  ✗ Error: {e}")
            continue
    
    return sequences

if __name__ == "__main__":
    print("Fetching TP53 protein sequences from multiple species...\n")
    
    sequences = fetch_tp53_sequences()
    
    print(f"\nTotal sequences fetched: {len(sequences)}")
    
    # Save to multi-FASTA file
    output_path = Path(__file__).parent / "tp53_multi_species.fasta"
    SeqIO.write(sequences, output_path, "fasta")
    print(f"Saved to: {output_path}")
