from __future__ import annotations
import gc
import sys
from pathlib import Path
from typing import Dict, List, Set, Tuple, Iterator

import numpy as np
import pandas as pd
import networkx as nx

# --- CONFIGURATION ---
INPUT_CSV = Path("data/work/StanaAndrei/lab06/gene_expr_matrix.csv")
OUTPUT_DIR = Path("labs/06_networks/submissions/StanaAndrei")
OUTPUT_CSV = OUTPUT_DIR / "modules_StanaAndrei.csv"

# RAM SAFETY SETTINGS
CHUNK_SIZE = 2000          # Read 2000 lines at a time (low memory footprint)
TARGET_GENE_COUNT = 3000   # Keep top 3000 genes (approx 9M cells = 0.07GB RAM)

# ANALYSIS SETTINGS
CORR_METHOD = "spearman"
ADJ_THRESHOLD = 0.7        # Lowered slightly from 0.8 to ensure edges exist
# ---------------------

def get_top_genes_iteratively(path: Path, top_n: int) -> List[str]:
    """
    PASS 1: Iterate through file, calculate variance per chunk, 
    discard data immediately to save RAM.
    """
    print(f"   [Pass 1] Scanning {path} for high-variance genes...")
    
    # Store (gene_name, variance) pairs
    gene_variances: List[Tuple[str, float]] = []
    
    # Create an iterator that reads the file in pieces
    chunk_iter = pd.read_csv(path, index_col=0, chunksize=CHUNK_SIZE)
    
    for i, chunk in enumerate(chunk_iter):
        # 1. Convert to float32 (save memory)
        chunk = chunk.astype(np.float32)
        
        # 2. Log transform
        np.log2(chunk + 1, out=chunk)
        
        # 3. Calculate variance for this chunk of genes
        vars = chunk.var(axis=1)
        
        # 4. Save results
        # zip creates an iterator of (index, value) tuples
        gene_variances.extend(zip(vars.index, vars.values))
        
        # 5. Manual Cleanup
        del chunk, vars
        if i % 5 == 0:
            gc.collect()
            print(f"      ... scanned {len(gene_variances)} genes so far")

    print(f"   [Pass 1] Completed. Total genes scanned: {len(gene_variances)}")
    
    # Sort by variance (descending) and take top N
    # x[1] is the variance value
    gene_variances.sort(key=lambda x: x[1], reverse=True)
    top_genes = [x[0] for x in gene_variances[:top_n]]
    
    return top_genes

def load_selected_genes(path: Path, keep_genes: Set[str]) -> pd.DataFrame:
    """
    PASS 2: Read file again, but ONLY keep rows that are in 'keep_genes'.
    """
    print(f"   [Pass 2] Harvesting top {len(keep_genes)} genes from file...")
    
    filtered_chunks = []
    chunk_iter = pd.read_csv(path, index_col=0, chunksize=CHUNK_SIZE)
    
    for chunk in chunk_iter:
        # Keep only rows that are in our target set
        # intersection is fast
        common_genes = chunk.index.intersection(keep_genes)
        
        if not common_genes.empty:
            # Append only the small slice
            subset = chunk.loc[common_genes].astype(np.float32)
            # Log transform now so we don't have to do it on the big matrix later
            np.log2(subset + 1, out=subset)
            filtered_chunks.append(subset)
            
        del chunk
        
    if not filtered_chunks:
        raise ValueError("Error: No matching genes found in Pass 2. Check CSV integrity.")

    # Combine the small pieces into one DataFrame
    df_final = pd.concat(filtered_chunks)
    return df_final

def generate_edges(corr_df: pd.DataFrame, threshold: float) -> Iterator[Tuple[str, str]]:
    """
    Generator function that yields edges one by one using Numpy.
    Avoids creating a massive list of edges in memory.
    """
    matrix = corr_df.values
    genes = corr_df.index.to_numpy()
    n = len(genes)
    
    # Iterate over the upper triangle only
    # We use numpy indices to find valid spots
    mask = np.triu(np.ones(matrix.shape, dtype=bool), k=1)
    
    # Find indices where condition is met
    # rows, cols will be numpy arrays of indices
    rows, cols = np.where((np.abs(matrix) >= threshold) & mask)
    
    total_edges = len(rows)
    print(f"   -> Found {total_edges} edges. Streaming to graph...")
    
    for r, c in zip(rows, cols):
        yield (genes[r], genes[c])

if __name__ == "__main__":
    gc.enable()
    
    try:
        # 1. PASS 1: Find best genes without loading file
        top_gene_names = get_top_genes_iteratively(INPUT_CSV, TARGET_GENE_COUNT)
        top_gene_set = set(top_gene_names) # Convert to set for O(1) lookups
        
        # 2. PASS 2: Load only those genes
        expr_data = load_selected_genes(INPUT_CSV, top_gene_set)
        print(f"   Data loaded successfully. Shape: {expr_data.shape} (Genes x Probes)")
        
        # 3. Correlation (Safe size: 3000 x 3000)
        print("   Calculating Correlation Matrix...")
        corr = expr_data.T.corr(method=CORR_METHOD)
        
        # Free expression data
        del expr_data
        gc.collect()
        
        # 4. Build Graph using Iterator
        print(f"   Building Graph (Threshold >= {ADJ_THRESHOLD})...")
        G = nx.Graph()
        G.add_nodes_from(corr.index) # Add nodes first (handles isolates)
        
        # Consume the generator
        edge_generator = generate_edges(corr, ADJ_THRESHOLD)
        G.add_edges_from(edge_generator)
        
        # Free correlation matrix
        del corr
        gc.collect()
        
        # Remove isolates to clean up result
        G.remove_nodes_from(list(nx.isolates(G)))
        
        print(f"   Final Graph: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges.")
        
        # 5. Detect Modules
        if G.number_of_nodes() > 0:
            print("   Detecting modules...")
            try:
                from networkx.algorithms.community import louvain_communities
                comms = louvain_communities(G, seed=42)
            except ImportError:
                from networkx.algorithms.community import greedy_modularity_communities
                comms = greedy_modularity_communities(G)
                
            mapping = {}
            for i, c in enumerate(comms, 1):
                for gene in c:
                    mapping[gene] = i
            
            OUTPUT_CSV.parent.mkdir(parents=True, exist_ok=True)
            pd.DataFrame({
                "Gene": list(mapping.keys()), 
                "Module": list(mapping.values())
            }).sort_values(["Module", "Gene"]).to_csv(OUTPUT_CSV, index=False)
            print(f"   SUCCESS! Saved to {OUTPUT_CSV}")
        else:
            print("   Warning: Graph is empty. Try lowering ADJ_THRESHOLD.")

    except Exception as e:
        print(f"\nCRITICAL FAILURE: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)