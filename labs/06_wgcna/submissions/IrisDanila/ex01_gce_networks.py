
from __future__ import annotations
from pathlib import Path
from typing import Dict, Iterable

import numpy as np
import pandas as pd
import networkx as nx


# --------------------------
# Config
# --------------------------
INPUT_CSV = Path("data/work/IrisDanila/lab06/expression_matrix.csv")
OUTPUT_DIR = Path("labs/06_wgcna/submissions/IrisDanila")
OUTPUT_CSV = OUTPUT_DIR / "modules_IrisDanila.csv"

CORR_METHOD = "spearman"   # Spearman correlation (robust to non-linearity)
VARIANCE_THRESHOLD = 0.5   # Filter low-variance genes
ADJ_THRESHOLD = 0.6        # Correlation threshold for edges
USE_ABS_CORR = True        # Use absolute correlation
MAKE_UNDIRECTED = True     # Undirected network (co-expression is symmetric)


def read_expression_matrix(path: Path) -> pd.DataFrame:
    """
    Citește matricea de expresie din CSV.
    Format: rânduri = gene, coloane = samples
    """
    if not path.exists():
        raise FileNotFoundError(
            f"Nu am găsit {path}. Rulați create_expression_data.py mai întâi."
        )
    df = pd.read_csv(path, index_col=0)
    if df.empty:
        raise ValueError("Matricea de expresie este goală.")
    
    print(f"Loaded expression matrix: {df.shape[0]} genes × {df.shape[1]} samples")
    return df


def log_and_filter(df: pd.DataFrame, variance_threshold: float) -> pd.DataFrame:
    """
    Preprocesare:
    - Aplică log2(x+1) transformare
    - Filtrează genele cu varianță scăzută
    """
    print("\n" + "="*70)
    print("PREPROCESSING")
    print("="*70)
    
    # Log transformation
    print("Applying log2(x+1) transformation...")
    df_log = np.log2(df + 1)
    
    # Calculate variance across samples for each gene
    gene_variances = df_log.var(axis=1)
    
    # Filter low-variance genes
    print(f"Filtering genes with variance <= {variance_threshold}...")
    df_filt = df_log.loc[gene_variances > variance_threshold]
    
    n_removed = len(df_log) - len(df_filt)
    print(f"  Removed: {n_removed} genes")
    print(f"  Retained: {len(df_filt)} genes")
    
    return df_filt


def correlation_matrix(df: pd.DataFrame, 
                       method: str = "spearman",
                       use_abs: bool = True) -> pd.DataFrame:
    """
    Calculează matricea de corelație între gene.
    
    Args:
        df: Expression matrix (genes × samples)
        method: 'pearson' or 'spearman'
        use_abs: If True, use absolute correlation values
    
    Returns:
        Correlation matrix (genes × genes)
    """
    print("\n" + "="*70)
    print("CORRELATION MATRIX")
    print("="*70)
    
    print(f"Computing {method.capitalize()} correlation between genes...")
    
    # df is (genes × samples), we want correlation between genes
    # So we transpose and correlate
    corr = df.T.corr(method=method)
    
    if use_abs:
        print("Using absolute correlation values...")
        corr = corr.abs()
    
    print(f"Correlation matrix shape: {corr.shape}")
    print(f"Correlation range: [{corr.min().min():.3f}, {corr.max().max():.3f}]")
    
    return corr


def adjacency_from_correlation(corr: pd.DataFrame,
                               threshold: float,
                               weighted: bool = False) -> pd.DataFrame:
    """
    Construiește matricea de adiacență din corelații.
    
    Args:
        corr: Correlation matrix
        threshold: Minimum correlation for an edge
        weighted: If True, use correlation as edge weight; else binary
    
    Returns:
        Adjacency matrix
    """
    print("\n" + "="*70)
    print("ADJACENCY MATRIX")
    print("="*70)
    
    print(f"Threshold: {threshold}")
    print(f"Weighted: {weighted}")
    
    if weighted:
        A = corr.copy()
        A[A < threshold] = 0
    else:
        A = (corr >= threshold).astype(int)
    
    # Remove self-loops
    np.fill_diagonal(A.values, 0.0)
    
    n_edges = (A > 0).sum().sum() // 2  # Divide by 2 for undirected
    print(f"Edges created: {n_edges}")
    
    return A


def graph_from_adjacency(A: pd.DataFrame, undirected: bool = True) -> nx.Graph:
    """
    Construiește graf NetworkX din matricea de adiacență.
    """
    print("\n" + "="*70)
    print("GRAPH CONSTRUCTION")
    print("="*70)
    
    if undirected:
        G = nx.from_pandas_adjacency(A)
        print("Graph type: Undirected")
    else:
        G = nx.from_pandas_adjacency(A, create_using=nx.DiGraph)
        print("Graph type: Directed")
    
    # Remove isolated nodes
    isolates = list(nx.isolates(G))
    if isolates:
        print(f"Removing {len(isolates)} isolated nodes...")
        G.remove_nodes_from(isolates)
    
    print(f"Final graph: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")
    
    # Network statistics
    if G.number_of_nodes() > 0:
        density = nx.density(G)
        print(f"Network density: {density:.4f}")
        
        if nx.is_connected(G):
            print("Network is connected")
            avg_path_length = nx.average_shortest_path_length(G)
            print(f"Average shortest path length: {avg_path_length:.2f}")
        else:
            components = list(nx.connected_components(G))
            print(f"Network has {len(components)} connected components")
            largest = max(components, key=len)
            print(f"Largest component: {len(largest)} nodes")
    
    return G


def detect_modules_louvain_or_greedy(G: nx.Graph) -> Dict[str, int]:
    """
    Detectează module (comunități) în rețea folosind Louvain sau greedy modularity.
    
    Returns:
        Dictionary mapping gene -> module_id
    """
    print("\n" + "="*70)
    print("MODULE DETECTION")
    print("="*70)
    
    # Try Louvain first (better algorithm)
    try:
        from networkx.algorithms.community import louvain_communities
        print("Using Louvain algorithm...")
        communities_iterable = louvain_communities(G, seed=42)
        communities = [set(c) for c in communities_iterable]
    except (ImportError, AttributeError):
        # Fallback to greedy modularity (always available)
        print("Louvain not available, using greedy modularity...")
        from networkx.algorithms.community import greedy_modularity_communities
        communities_iterable: Iterable[Iterable[str]] = greedy_modularity_communities(G)
        communities = [set(c) for c in communities_iterable]
    
    # Create gene -> module mapping
    mapping: Dict[str, int] = {}
    for midx, comm in enumerate(communities, start=1):
        for gene in comm:
            mapping[gene] = midx
    
    print(f"Detected {len(communities)} modules")
    
    # Module size distribution
    module_sizes = [len(c) for c in communities]
    print(f"Module sizes: min={min(module_sizes)}, max={max(module_sizes)}, mean={np.mean(module_sizes):.1f}")
    
    # Show module composition
    for midx, comm in enumerate(communities[:5], start=1):  # Show first 5 modules
        print(f"  Module {midx}: {len(comm)} genes")
    
    if len(communities) > 5:
        print(f"  ... and {len(communities) - 5} more modules")
    
    return mapping


def save_modules_csv(mapping: Dict[str, int], out_csv: Path) -> None:
    """
    Salvează mapping-ul gene → module în CSV.
    """
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    
    df_modules = (
        pd.DataFrame({
            "Gene": list(mapping.keys()), 
            "Module": list(mapping.values())
        })
        .sort_values(["Module", "Gene"])
    )
    
    df_modules.to_csv(out_csv, index=False)
    
    print(f"\n✓ Gene-to-module mapping saved to: {out_csv}")
    print(f"\nFirst 10 genes:")
    print(df_modules.head(10))


def main():
    print("="*70)
    print("LAB 6 - GENE CO-EXPRESSION NETWORK ANALYSIS")
    print("="*70)
    
    # 1. Load expression data
    expr = read_expression_matrix(INPUT_CSV)
    
    # 2. Preprocess: log transform and filter low-variance genes
    expr_pp = log_and_filter(expr, variance_threshold=VARIANCE_THRESHOLD)
    
    # 3. Calculate correlation matrix
    corr = correlation_matrix(expr_pp, method=CORR_METHOD, use_abs=USE_ABS_CORR)
    
    # 4. Build adjacency matrix
    adj = adjacency_from_correlation(corr, threshold=ADJ_THRESHOLD, weighted=False)
    
    # 5. Create graph
    G = graph_from_adjacency(adj, undirected=MAKE_UNDIRECTED)
    
    # 6. Detect modules (communities)
    gene_to_module = detect_modules_louvain_or_greedy(G)
    
    # 7. Save results
    save_modules_csv(gene_to_module, OUTPUT_CSV)
    
    print("\n" + "="*70)
    print("ANALYSIS COMPLETE!")
    print("="*70)
    print(f"\nResults saved to: {OUTPUT_CSV}")
    print(f"Detected {len(set(gene_to_module.values()))} co-expression modules")
    print("="*70)


if __name__ == "__main__":
    main()
