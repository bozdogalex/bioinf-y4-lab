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
HANDLE = "rauldavid35"  # TODO: Replace with your actual handle

# Input files
EXPR_CSV = Path(f"data/work/{HANDLE}/lab06/expression_matrix.csv")
MODULES_CSV = Path(f"labs/06_wgcna/submissions/{HANDLE}/modules_{HANDLE}.csv")

# Optional: if you saved adjacency in Lab 6, load it here
PRECOMPUTED_ADJ_CSV: Optional[Path] = None

# Parameters for adjacency reconstruction (if PRECOMPUTED_ADJ_CSV is None)
CORR_METHOD = "spearman"
USE_ABS_CORR = True        
ADJ_THRESHOLD = 0.6        
WEIGHTED = False           

# Visualization parameters
SEED = 42
TOPK_HUBS = 10
NODE_BASE_SIZE = 60
EDGE_ALPHA = 0.15

# Outputs
OUT_DIR = Path(f"labs/07_network_viz/submissions/{HANDLE}")
OUT_PNG = OUT_DIR / f"network_{HANDLE}.png"
OUT_HUBS = OUT_DIR / f"hubs_{HANDLE}.csv"


# --------------------------
# Utils
# --------------------------
def ensure_exists(path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Input file not found: {path}")


def read_expression_matrix(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, index_col=0)


def read_modules_csv(path: Path) -> Dict[str, int]:
    df = pd.read_csv(path)
    return pd.Series(df.Module.values, index=df.Gene).to_dict()


def correlation_to_adjacency(expr: pd.DataFrame,
                             method: str,
                             use_abs: bool,
                             threshold: float,
                             weighted: bool) -> pd.DataFrame:

    corr = expr.T.corr(method=method)

    if use_abs:
        corr = corr.abs()

    if weighted:
        adj = corr.where(corr >= threshold, 0.0)
    else:
        adj = (corr >= threshold).astype(int)

    np.fill_diagonal(adj.values, 0)

    return adj


def graph_from_adjacency(A: pd.DataFrame) -> nx.Graph:

    G = nx.from_pandas_adjacency(A)
    
    isolates = list(nx.isolates(G))
    G.remove_nodes_from(isolates)
    
    return G


def color_map_from_modules(nodes: Iterable[str], gene2module: Dict[str, int]) -> list:

    cmap = plt.get_cmap("tab10")
    node_colors = []
    for n in nodes:
        mod_id = gene2module.get(n, 0)
        color = cmap((mod_id - 1) % 10) 
        node_colors.append(color)
    
    return node_colors


def compute_hubs(G: nx.Graph, topk: int) -> pd.DataFrame:

    deg = dict(G.degree())
    
    hubs = sorted(deg.items(), key=lambda x: x[1], reverse=True)[:topk]
    
    return pd.DataFrame(hubs, columns=["Gene", "Degree"])


# --------------------------
# Main
# --------------------------
if __name__ == "__main__":
    print(f"--- Running Network Viz for handle: {HANDLE} ---")

    # TODO 1: Verify input files exist
    if HANDLE == "<handle>":
        print("WARNING: Please set the HANDLE variable in the config section.")
    
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # TODO 2: Load expression matrix and module mapping
    print("Loading data...")
    try:
        expr = read_expression_matrix(EXPR_CSV)
        gene2module = read_modules_csv(MODULES_CSV)
        print(f"Loaded expression matrix: {expr.shape}")
        print(f"Loaded {len(gene2module)} module assignments")
    except FileNotFoundError:
        print("Data files not found. Ensure paths are correct.")
        exit(1)

    # TODO 3: Load or reconstruct adjacency
    print("Constructing adjacency matrix...")
    if PRECOMPUTED_ADJ_CSV and PRECOMPUTED_ADJ_CSV.exists():
        print(f"Loading precomputed adjacency from {PRECOMPUTED_ADJ_CSV}")
        adj = pd.read_csv(PRECOMPUTED_ADJ_CSV, index_col=0)
        common_genes = adj.index.intersection(list(gene2module.keys()))
        adj = adj.loc[common_genes, common_genes]
    else:
        print(f"Computing correlation ({CORR_METHOD}) and applying threshold ({ADJ_THRESHOLD})...")
        valid_genes = [g for g in gene2module.keys() if g in expr.index]
        expr_filtered = expr.loc[valid_genes]
        adj = correlation_to_adjacency(expr_filtered, CORR_METHOD, USE_ABS_CORR, ADJ_THRESHOLD, WEIGHTED)

    # TODO 4: Build graph
    print("Building NetworkX graph...")
    G = graph_from_adjacency(adj)
    print(f"Graph stats: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")

    # TODO 5: Compute colors by module
    print("Assigning node colors...")
    node_colors = color_map_from_modules(G.nodes(), gene2module)

    # TODO 6: Compute hub genes
    print(f"Identifying top {TOPK_HUBS} hub genes...")
    hubs_df = compute_hubs(G, TOPK_HUBS)
    hub_nodes = set(hubs_df["Gene"])
    print(hubs_df)

    node_sizes = [NODE_BASE_SIZE * 3 if n in hub_nodes else NODE_BASE_SIZE for n in G.nodes()]

    # TODO 7: Compute layout and draw graph
    print("Computing layout (Spring Layout)...")
    pos = nx.spring_layout(G, seed=SEED, k=0.15, iterations=50)

    plt.figure(figsize=(10, 10))
    
    nx.draw_networkx_edges(G, pos, alpha=EDGE_ALPHA, edge_color="gray")
    
    nx.draw_networkx_nodes(
        G, pos,
        node_color=node_colors,
        node_size=node_sizes,
        edgecolors="white", # Thin border for nodes
        linewidths=0.5
    )

    hub_labels = {n: n for n in G.nodes() if n in hub_nodes}
    nx.draw_networkx_labels(G, pos, labels=hub_labels, font_size=10, font_weight="bold")

    plt.title(f"Gene Co-Expression Network ({HANDLE})\nThreshold={ADJ_THRESHOLD}, {CORR_METHOD}")
    plt.axis("off")
    plt.tight_layout()

    # TODO 8: Save network figure
    print(f"Saving figure to {OUT_PNG}...")
    plt.savefig(OUT_PNG, dpi=300, bbox_inches='tight')

    # TODO 9: Save hub genes to CSV
    print(f"Saving hub genes to {OUT_HUBS}...")
    hubs_df.to_csv(OUT_HUBS, index=False)

    print("Done!")