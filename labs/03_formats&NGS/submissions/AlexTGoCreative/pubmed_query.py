#!/usr/bin/env python3
"""
PubMed Query Script
Searches PubMed for 'TP53 AND cancer' and retrieves article information.
"""

from Bio import Entrez
import sys

# Configure Entrez
Entrez.email = "alexandru.tulbure@student.upt.ro"
HANDLE = "AlexTGoCreative"


def search_pubmed(query, max_results=5):
    """
    Search PubMed with a query and return article IDs.
    
    Args:
        query (str): Search query string
        max_results (int): Maximum number of results to return
        
    Returns:
        list: List of PubMed IDs
    """
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results)
        record = Entrez.read(handle)
        handle.close()
        return record["IdList"]
    except Exception as e:
        print(f"Error searching PubMed: {e}")
        return []


def fetch_article_details(pmid_list):
    """
    Fetch article details for a list of PubMed IDs.
    
    Args:
        pmid_list (list): List of PubMed IDs
        
    Returns:
        list: List of article dictionaries with title, authors, and abstract
    """
    articles = []
    
    try:
        # Fetch article details
        handle = Entrez.efetch(db="pubmed", id=pmid_list, retmode="xml")
        records = Entrez.read(handle)
        handle.close()
        
        for record in records['PubmedArticle']:
            article = {}
            
            # Extract title
            try:
                article['title'] = record['MedlineCitation']['Article']['ArticleTitle']
            except KeyError:
                article['title'] = "Title not available"
            
            # Extract authors
            try:
                author_list = record['MedlineCitation']['Article']['AuthorList']
                authors = []
                for author in author_list:
                    if 'LastName' in author and 'Initials' in author:
                        authors.append(f"{author['LastName']} {author['Initials']}")
                    elif 'CollectiveName' in author:
                        authors.append(author['CollectiveName'])
                article['authors'] = ", ".join(authors) if authors else "Authors not available"
            except (KeyError, TypeError):
                article['authors'] = "Authors not available"
            
            # Extract abstract
            try:
                abstract = record['MedlineCitation']['Article']['Abstract']['AbstractText']
                if isinstance(abstract, list):
                    article['abstract'] = " ".join([str(text) for text in abstract])
                else:
                    article['abstract'] = str(abstract)
            except KeyError:
                article['abstract'] = "Abstract not available"
            
            # Extract PMID
            try:
                article['pmid'] = str(record['MedlineCitation']['PMID'])
            except KeyError:
                article['pmid'] = "PMID not available"
            
            articles.append(article)
            
    except Exception as e:
        print(f"Error fetching article details: {e}")
    
    return articles


def save_results(articles, filename):
    """
    Save article information to a text file.
    
    Args:
        articles (list): List of article dictionaries
        filename (str): Output filename
    """
    try:
        with open(filename, 'w', encoding='utf-8') as f:
            f.write("=" * 80 + "\n")
            f.write("PubMed Query Results: TP53 AND cancer\n")
            f.write("=" * 80 + "\n\n")
            
            for i, article in enumerate(articles, 1):
                f.write(f"Article {i}\n")
                f.write("-" * 80 + "\n")
                f.write(f"PMID: {article.get('pmid', 'N/A')}\n\n")
                f.write(f"Title: {article.get('title', 'N/A')}\n\n")
                f.write(f"Authors: {article.get('authors', 'N/A')}\n\n")
                f.write(f"Abstract:\n{article.get('abstract', 'N/A')}\n\n")
                f.write("=" * 80 + "\n\n")
        
        print(f"Results saved to {filename}")
        
    except Exception as e:
        print(f"Error saving results: {e}")


def main():
    """Main function to execute PubMed query and save results."""
    query = "TP53 AND cancer"
    max_results = 5
    output_file = f"pubmed_{HANDLE}.txt"
    
    print(f"Searching PubMed for: '{query}'")
    print(f"Maximum results: {max_results}\n")
    
    # Search PubMed
    pmid_list = search_pubmed(query, max_results)
    
    if not pmid_list:
        print("No results found or error occurred during search.")
        sys.exit(1)
    
    print(f"Found {len(pmid_list)} articles")
    print(f"PubMed IDs: {', '.join(pmid_list)}\n")
    
    # Fetch article details
    print("Fetching article details...")
    articles = fetch_article_details(pmid_list)
    
    if not articles:
        print("No article details retrieved.")
        sys.exit(1)
    
    print(f"Retrieved {len(articles)} articles\n")
    
    # Save results
    save_results(articles, output_file)
    
    print("\nDone!")


if __name__ == "__main__":
    main()
