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
import matplotlib.colors as mcolors
import gc

# --------------------------
# Config — complete with your values
# --------------------------
HANDLE = "coroiusamy"

# Input files
EXPR_CSV = Path(f"data/work/{HANDLE}/lab06/expression_matrix.csv")
MODULES_CSV = Path(f"/workspaces/bioinf-y4-lab/labs/07_network_viz/submissions/{HANDLE}/modules_coroiusamy.csv")

PRECOMPUTED_ADJ_CSV: Optional[Path] = None

# Parameters for adjacency reconstruction 
CORR_METHOD = "spearman"   # TODO: choose "pearson" or "spearman"
USE_ABS_CORR = True        # TODO: use absolute correlations?
ADJ_THRESHOLD = 0.70       # TODO: correlation threshold
WEIGHTED = False           # TODO: use weighted adjacency or binary?

# Visualization parameters
SEED = 42
TOPK_HUBS = 10
NODE_BASE_SIZE = 60
EDGE_ALPHA = 0.15
MAX_GENES = 1500  # to prevent memory crash

# Outputs
OUT_DIR = Path(f"labs/07_network_viz/submissions/{HANDLE}")
OUT_PNG = OUT_DIR / f"network_{HANDLE}.png"
OUT_HUBS = OUT_DIR / f"hubs_{HANDLE}.csv"


# --------------------------
# Utils
# --------------------------
def ensure_exists(path: Path) -> None:
    """TODO: check that a file exists."""
    if not path.exists():
        raise FileNotFoundError(f"Input file not found: {path}")
    print(f"[OK] Found {path}")


def read_expression_matrix(path: Path) -> pd.DataFrame:
    """
    TODO:
    - read CSV
    - set index to gene names
    """
    print(f"Reading expression matrix from {path}...")
    df = pd.read_csv(path, index_col=0)
    df = df.select_dtypes(include=[np.number])
    return df


def read_modules_csv(path: Path) -> Dict[str, int]:
    """
    TODO:
    - read CSV with columns: Gene, Module
    - return dict: gene -> module_id
    """
    print(f"Reading modules from {path}...")
    df = pd.read_csv(path)
    # Handle variable column names
    if 'Gene' in df.columns and 'Module' in df.columns:
        return pd.Series(df.Module.values, index=df.Gene).to_dict()
    else:
        cols = df.columns
        return pd.Series(df[cols[1]].values, index=df[cols[0]]).to_dict()


def correlation_to_adjacency(expr: pd.DataFrame,
                             method: str,
                             use_abs: bool,
                             threshold: float,
                             weighted: bool) -> pd.DataFrame:
    """
    TODO:
    - compute correlation matrix on expr
    - optionally apply abs()
    - apply threshold to build adjacency
    - remove diagonal
    """
    print(f"Computing {method} correlation on {expr.shape[0]} genes...")
    # Compute correlation
    corr = expr.T.corr(method=method)
    
    # Memory cleanup
    gc.collect()

    if use_abs:
        corr = corr.abs()
    
    # Remove diagonal
    np.fill_diagonal(corr.values, 0)
    
    print(f"Applying threshold > {threshold}...")
    if weighted:
        adj = corr.where(corr > threshold, 0)
    else:
        # Convert to int8 to save memory
        adj = (corr > threshold).astype(np.int8)
    
    # Clean up correlation matrix
    del corr
    gc.collect()
        
    return adj


def graph_from_adjacency(A: pd.DataFrame) -> nx.Graph:
    """
    TODO:
    - convert adjacency DataFrame to NetworkX graph
    - remove isolated nodes
    """
    print("Building NetworkX graph...")
    G = nx.from_pandas_adjacency(A)
    
    isolates = list(nx.isolates(G))
    G.remove_nodes_from(isolates)
    
    print(f"Graph stats: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")
    return G


def color_map_from_modules(nodes: Iterable[str], gene2module: Dict[str, int]) -> Dict[str, str]:
    """
    TODO:
    - assign a color to each node based on its module
    - use matplotlib 'tab10' or another palette
    """
    unique_modules = sorted(list(set(gene2module.values())))
    cmap = plt.get_cmap('tab20')
    
    colors = {}
    for node in nodes:
        mod = gene2module.get(node)
        if mod is not None:
            colors[node] = mcolors.to_hex(cmap(mod % 20))
        else:
            colors[node] = '#cccccc'
    return colors


def compute_hubs(G: nx.Graph, topk: int) -> pd.DataFrame:
    """
    TODO:
    - compute degree for every node
    - (optional) compute betweenness centrality
    - return top-k genes
    """
    degree = dict(G.degree())
    sorted_deg = sorted(degree.items(), key=lambda x: x[1], reverse=True)
    return pd.DataFrame(sorted_deg, columns=['Gene', 'Degree']).head(topk)


# --------------------------
# Main
# --------------------------
if __name__ == "__main__":
    # TODO 1: Verify input files exist
    ensure_exists(MODULES_CSV)
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # TODO 2: Load expression matrix and module mapping
    full_expr = read_expression_matrix(EXPR_CSV)
    gene2module = read_modules_csv(MODULES_CSV)

    # --- MEMORY OPTIMIZATION START ---
    # Filter expression matrix to match only module genes
    valid_genes = [g for g in gene2module.keys() if g in full_expr.index]
    expr = full_expr.loc[valid_genes].copy()
    
    del full_expr # Free memory
    gc.collect()

    # Preprocessing: Log2 and float32 conversion
    expr = np.log2(expr.astype('float32') + 1)

    # Variance filter if too large
    if expr.shape[0] > MAX_GENES:
        print(f"Filtering to top {MAX_GENES} genes by variance to prevent OOM...")
        vars = expr.var(axis=1)
        top_genes = vars.nlargest(MAX_GENES).index
        expr = expr.loc[top_genes]
    # --- MEMORY OPTIMIZATION END ---

    # TODO 3: Load or reconstruct adjacency
    if PRECOMPUTED_ADJ_CSV and PRECOMPUTED_ADJ_CSV.exists():
        pass # implementation skipped as per config
    else:
        adj_matrix = correlation_to_adjacency(
            expr,
            method=CORR_METHOD,
            use_abs=USE_ABS_CORR,
            threshold=ADJ_THRESHOLD,
            weighted=WEIGHTED
        )

    # TODO 4: Build graph
    G = graph_from_adjacency(adj_matrix)

    if G.number_of_nodes() == 0:
        print("Warning: Graph is empty. Try lowering ADJ_THRESHOLD.")
        exit()

    # TODO 5: Compute colors by module
    color_map = color_map_from_modules(G.nodes(), gene2module)
    node_colors = [color_map[n] for n in G.nodes()]

    # TODO 6: Compute hub genes
    hubs_df = compute_hubs(G, TOPK_HUBS)
    print("Top Hubs:")
    print(hubs_df)
    
    hub_genes = set(hubs_df['Gene'])
    node_sizes = [NODE_BASE_SIZE * 4 if n in hub_genes else NODE_BASE_SIZE/2 for n in G.nodes()]

    # TODO 7: Compute layout and draw graph
    plt.figure(figsize=(12, 12))
    print("Computing layout...")
    pos = nx.spring_layout(G, seed=SEED, k=0.2, iterations=30)
    
    # - draw edges
    nx.draw_networkx_edges(G, pos, alpha=EDGE_ALPHA, edge_color='gray')
    # - draw nodes (colored)
    nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=node_sizes, linewidths=0.5, edgecolors='white')
    # - draw labels only for hubs
    hub_labels = {n: n for n in G.nodes() if n in hub_genes}
    nx.draw_networkx_labels(G, pos, labels=hub_labels, font_size=10, font_weight='bold')

    plt.title(f"Network (Top {MAX_GENES} var. genes) | Thresh: {ADJ_THRESHOLD}")
    plt.axis('off')

    # TODO 8: Save network figure
    plt.savefig(OUT_PNG, dpi=300, bbox_inches='tight')
    print(f"Saved: {OUT_PNG}")

    # TODO 9: Save hub genes to CSV
    hubs_df.to_csv(OUT_HUBS, index=False)
    print(f"Saved: {OUT_HUBS}")

    # --- BONUS: EXPORT PENTRU CYTOSCAPE/GEPHI ---
    # 1. Adăugăm atributul 'module' în graf
    nx.set_node_attributes(G, gene2module, name='module')
    
    # 2. Exportăm în format GraphML
    OUT_GRAPHML = OUT_DIR / f"network_{HANDLE}.graphml"
    nx.write_graphml(G, OUT_GRAPHML)
    print(f"Graph exported for Cytoscape: {OUT_GRAPHML}")
    
    # 3. Exportăm și ca listă de muchii (Edge List CSV)
    OUT_EDGES = OUT_DIR / f"edges_{HANDLE}.csv"
    nx.write_edgelist(G, OUT_EDGES, delimiter=',', data=False)
    print(f"Edges exported for manual import: {OUT_EDGES}")
    
    print("Done.")