"""
Exercise 8 — Visualization of Co-Expression Networks + Hub Genes

TODO:
- Load the expression matrix and module mapping from Lab 6
- Rebuild (or load) the adjacency matrix
- Construct the graph from adjacency
- Color nodes by module
- Compute hub genes (top degree)
- Visualize and export the network figure (.png)
- Export hub genes to CSV (optional)
"""

from __future__ import annotations
from pathlib import Path
from typing import Dict, Iterable, Optional

import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt


# --------------------------
# Config — complete with your values
# --------------------------
HANDLE = "AlexTGoCreative"

# Input files
EXPR_CSV = Path(f"../../../../data/sample/expression_data.csv")
MODULES_CSV = Path(f"../../../06_wgcna/submissions/{HANDLE}/modules_{HANDLE}.csv")

# Optional: if you saved adjacency in Lab 6, load it here
PRECOMPUTED_ADJ_CSV: Optional[Path] = None

# Parameters for adjacency reconstruction (if PRECOMPUTED_ADJ_CSV is None)
CORR_METHOD = "spearman"   # choose "pearson" or "spearman"
USE_ABS_CORR = True        # use absolute correlations?
ADJ_THRESHOLD = 0.6        # correlation threshold
WEIGHTED = False           # use weighted adjacency or binary?

# Visualization parameters
SEED = 42
TOPK_HUBS = 10
NODE_BASE_SIZE = 60
EDGE_ALPHA = 0.15

# Outputs
OUT_DIR = Path(f"")
OUT_PNG = OUT_DIR / f"network_{HANDLE}.png"
OUT_HUBS = OUT_DIR / f"hubs_{HANDLE}.csv"


# --------------------------
# Utils
# --------------------------
def ensure_exists(path: Path) -> None:
    """Check that a file exists."""
    if not path.exists():
        raise FileNotFoundError(f"File not found: {path}")
    print(f"✓ Found: {path}")


def read_expression_matrix(path: Path) -> pd.DataFrame:
    """
    Read expression matrix from CSV.
    First column should contain gene names (used as index).
    """
    df = pd.read_csv(path, index_col=0)
    print(f"Loaded expression matrix: {df.shape[0]} genes × {df.shape[1]} samples")
    return df


def read_modules_csv(path: Path) -> Dict[str, int]:
    """
    Read module mapping from CSV.
    Expected columns: Gene, Module
    Returns dict: gene -> module_id
    """
    df = pd.read_csv(path)
    gene2module = dict(zip(df['Gene'], df['Module']))
    print(f"Loaded {len(gene2module)} genes across {len(set(gene2module.values()))} modules")
    return gene2module


def correlation_to_adjacency(expr: pd.DataFrame,
                             method: str,
                             use_abs: bool,
                             threshold: float,
                             weighted: bool) -> pd.DataFrame:
    """
    Compute adjacency matrix from expression data.
    - compute correlation matrix between genes (rows)
    - optionally apply abs()
    - apply threshold to build adjacency
    - remove diagonal
    """
    # Log-transform
    expr_log = np.log2(expr + 1)
    
    # Compute correlation between genes (transpose to have genes as columns for .corr())
    corr = expr_log.T.corr(method=method)
    
    if use_abs:
        corr = corr.abs()
    
    # Build adjacency matrix
    if weighted:
        A = corr.copy()
        A[A < threshold] = 0
    else:
        A = (corr >= threshold).astype(int)
    
    # Remove diagonal (self-loops)
    np.fill_diagonal(A.values, 0)
    
    num_edges = (A > 0).sum().sum() // 2
    print(f"Adjacency matrix: {num_edges} edges above threshold {threshold}")
    return A


def graph_from_adjacency(A: pd.DataFrame) -> nx.Graph:
    """
    Convert adjacency DataFrame to NetworkX graph.
    Remove isolated nodes.
    """
    G = nx.from_pandas_adjacency(A)
    
    # Remove isolated nodes
    isolates = list(nx.isolates(G))
    if isolates:
        G.remove_nodes_from(isolates)
        print(f"Removed {len(isolates)} isolated nodes")
    
    print(f"Graph: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")
    return G


def color_map_from_modules(nodes: Iterable[str], gene2module: Dict[str, int]) -> Dict[str, str]:
    """
    Assign a color to each node based on its module.
    Uses matplotlib 'tab10' palette.
    """
    cmap = plt.get_cmap('tab10')
    color_map = {}
    for node in nodes:
        module_id = gene2module.get(node, 0)
        color_map[node] = cmap(module_id % 10)
    return color_map


def compute_hubs(G: nx.Graph, topk: int) -> pd.DataFrame:
    """
    Compute hub genes based on degree centrality.
    Returns top-k genes with highest degree.
    """
    # Compute degree for all nodes
    degrees = dict(G.degree())
    
    # Sort by degree (descending) and take top-k
    sorted_genes = sorted(degrees.items(), key=lambda x: x[1], reverse=True)
    top_genes = sorted_genes[:topk]
    
    # Create DataFrame
    df = pd.DataFrame(top_genes, columns=['Gene', 'Degree'])
    
    print(f"Top {topk} hub genes:")
    for gene, deg in top_genes:
        print(f"  {gene}: degree={deg}")
    
    return df


# --------------------------
# Main
# --------------------------
if __name__ == "__main__":
    print("=" * 70)
    print("Exercise 8 — Visualization of Co-Expression Networks + Hub Genes")
    print("=" * 70)
    
    # 1. Verify input files exist
    print("\n[1] Checking input files...")
    ensure_exists(EXPR_CSV)
    ensure_exists(MODULES_CSV)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    
    # 2. Load expression matrix and module mapping
    print("\n[2] Loading data...")
    expr = read_expression_matrix(EXPR_CSV)
    gene2module = read_modules_csv(MODULES_CSV)
    
    # 3. Load or reconstruct adjacency
    print("\n[3] Building adjacency matrix...")
    if PRECOMPUTED_ADJ_CSV and PRECOMPUTED_ADJ_CSV.exists():
        print(f"Loading precomputed adjacency from {PRECOMPUTED_ADJ_CSV}")
        A = pd.read_csv(PRECOMPUTED_ADJ_CSV, index_col=0)
        # Filter to genes in modules
        genes_in_modules = set(gene2module.keys())
        A = A.loc[A.index.isin(genes_in_modules), A.columns.isin(genes_in_modules)]
    else:
        print("Computing adjacency from correlations...")
        A = correlation_to_adjacency(expr, CORR_METHOD, USE_ABS_CORR, ADJ_THRESHOLD, WEIGHTED)
    
    # 4. Build graph
    print("\n[4] Building graph...")
    G = graph_from_adjacency(A)
    
    if G.number_of_nodes() == 0:
        print("⚠ Warning: Graph has no nodes after removing isolates!")
        print("Consider lowering the ADJ_THRESHOLD or using different data.")
        exit(1)
    
    # 5. Compute colors by module
    print("\n[5] Assigning colors by module...")
    node_color_map = color_map_from_modules(G.nodes(), gene2module)
    node_colors = [node_color_map[n] for n in G.nodes()]
    
    # 6. Compute hub genes
    print("\n[6] Computing hub genes...")
    hubs_df = compute_hubs(G, TOPK_HUBS)
    hub_nodes = set(hubs_df['Gene'].tolist())
    node_sizes = [NODE_BASE_SIZE * 3 if n in hub_nodes else NODE_BASE_SIZE for n in G.nodes()]
    
    # 7. Compute layout and draw graph
    print("\n[7] Drawing network...")
    pos = nx.spring_layout(G, seed=SEED)
    
    plt.figure(figsize=(12, 10))
    
    # Draw edges
    nx.draw_networkx_edges(G, pos, alpha=EDGE_ALPHA, width=0.5)
    
    # Draw nodes (colored by module)
    nx.draw_networkx_nodes(
        G, pos,
        node_color=node_colors,
        node_size=node_sizes,
        alpha=0.8
    )
    
    # Draw labels only for hubs
    hub_labels = {n: n for n in G.nodes() if n in hub_nodes}
    nx.draw_networkx_labels(G, pos, labels=hub_labels, font_size=9, font_weight='bold')
    
    plt.title(f"Co-Expression Network ({G.number_of_nodes()} genes, {G.number_of_edges()} edges)",
              fontsize=14, fontweight='bold')
    plt.axis('off')
    plt.tight_layout()
    
    # 8. Save network figure
    print("\n[8] Saving figure...")
    plt.savefig(OUT_PNG, dpi=200, bbox_inches='tight')
    print(f"✓ Saved network visualization to: {OUT_PNG}")
    
    # 9. Save hub genes to CSV
    print("\n[9] Saving hub genes...")
    hubs_df.to_csv(OUT_HUBS, index=False)
    print(f"✓ Saved hub genes to: {OUT_HUBS}")
    
    print("\n" + "=" * 70)
    print("✓ Visualization complete!")
    print("=" * 70)
    
    plt.show()
