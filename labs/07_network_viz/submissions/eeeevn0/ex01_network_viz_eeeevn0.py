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

HANDLE = "eeeevn0"

# Input files
EXPR_CSV = Path(f"data/work/{HANDLE}/lab06/expression_matrix.csv")
MODULES_CSV = Path(f"labs/07_network_viz/submissions/{HANDLE}/modules_{HANDLE}.csv")

# Optional: if you saved adjacency in Lab 6, load it here
PRECOMPUTED_ADJ_CSV: Optional[Path] = None

# Parameters for adjacency reconstruction (if PRECOMPUTED_ADJ_CSV is None)
CORR_METHOD = "spearman"   # TODO: choose "pearson" or "spearman"
USE_ABS_CORR = True        # TODO: use absolute correlations?
ADJ_THRESHOLD = 0.6        # TODO: correlation threshold
WEIGHTED = False           # TODO: use weighted adjacency or binary?

# Visualization parameters
SEED = 42
TOPK_HUBS = 10
NODE_BASE_SIZE = 60
EDGE_ALPHA = 0.15

# Outputs
OUT_DIR = Path(f"labs/07_network_viz/submissions/{HANDLE}")
OUT_PNG = OUT_DIR / f"network_{HANDLE}.png"
OUT_HUBS = OUT_DIR / f"hubs_{HANDLE}.csv"
OUT_EDGES = OUT_DIR / f"edges_{HANDLE}.csv"



def ensure_exists(path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Expected file not found: {path}")


def read_expression_matrix(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, index_col=0)
    return df


def read_modules_csv(path: Path) -> Dict[str, int]:
    df = pd.read_csv(path)
    if not {"Gene", "Module"}.issubset(df.columns):
        raise ValueError("modules CSV must contain columns: Gene, Module")
    return {row["Gene"]: int(row["Module"]) for _, row in df.iterrows()}


def correlation_to_adjacency(expr: pd.DataFrame,
                             method: str,
                             use_abs: bool,
                             threshold: float,
                             weighted: bool) -> pd.DataFrame:
    corr = expr.T.corr(method=method)
    if use_abs:
        corr = corr.abs()
    np.fill_diagonal(corr.values, 0.0)
    
    if weighted:
        A = corr.where(corr >= threshold, other=0.0)
    else:
        A = (corr >= threshold).astype(float)
    return A


def graph_from_adjacency(A: pd.DataFrame) -> nx.Graph:
    G = nx.from_pandas_adjacency(A)
    isolated = list(nx.isolates(G))
    G.remove_nodes_from(isolated)
    return G


def color_map_from_modules(nodes: Iterable[str], gene2module: Dict[str, int]) -> Dict[str, str]:
    nodes = list(nodes)
    modules_present = sorted({gene2module[n] for n in nodes if n in gene2module})
    if not modules_present:
        return {n: "#d3d3d3" for n in nodes}

    cmap = plt.get_cmap("tab20", len(modules_present))
    module_to_color = {m: cmap(i) for i, m in enumerate(modules_present)}

    node_to_color: Dict[str, str] = {}
    for n in nodes:
        mod = gene2module.get(n, None)
        if mod is None:
            node_to_color[n] = "#d3d3d3"  
        else:
            node_to_color[n] = module_to_color[mod]

    return node_to_color


def compute_hubs(G: nx.Graph, topk: int) -> pd.DataFrame:
    if G.number_of_nodes() == 0:
        return pd.DataFrame(columns=["Gene", "Degree", "Betweenness"])
    deg_dict = dict(G.degree())
    betw_dict = nx.betweenness_centrality(G)
    df = pd.DataFrame({
        "Gene": list(deg_dict.keys()),
        "Degree": list(deg_dict.values()),
        "Betweenness": [betw_dict[g] for g in deg_dict.keys()],
    })
    df = df.sort_values(by=["Degree", "Betweenness"], ascending=False)
    return df.head(topk)


if __name__ == "__main__":
    # TODO 1: Verify input files exist
    ensure_exists(EXPR_CSV)
    ensure_exists(MODULES_CSV)
    if PRECOMPUTED_ADJ_CSV is not None:
        ensure_exists(PRECOMPUTED_ADJ_CSV)

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # TODO 2: Load expression matrix and module mapping
    expr = read_expression_matrix(EXPR_CSV)
    gene2module = read_modules_csv(MODULES_CSV)

    # TODO 3: Load or reconstruct adjacency
    if PRECOMPUTED_ADJ_CSV is not None:
        A = pd.read_csv(PRECOMPUTED_ADJ_CSV, index_col=0)
        genes = expr.index.intersection(A.index)
        A = A.loc[genes, genes]
    else:
        A = correlation_to_adjacency(
            expr=expr,
            method=CORR_METHOD,
            use_abs=USE_ABS_CORR,
            threshold=ADJ_THRESHOLD,
            weighted=WEIGHTED,
        )

    genes_in_modules = [g for g in A.index if g in gene2module]
    if genes_in_modules:
        A = A.loc[genes_in_modules, genes_in_modules]

    # TODO 4: Build graph
    G = graph_from_adjacency(A)
    print(f"Graph has {G.number_of_nodes()} nodes and {G.number_of_edges()} edges")

    if G.number_of_nodes() == 0:
        print("Graph is empty after filtering. Nothing to visualize.")
    else:
        # TODO 5: Compute colors by module
        node_color_map = color_map_from_modules(G.nodes, gene2module)
        node_colors = [node_color_map[n] for n in G.nodes]

        # TODO 6: Compute hub genes
        hubs_df = compute_hubs(G, TOPK_HUBS)
        deg_dict = dict(G.degree())
        max_deg = max(deg_dict.values()) if deg_dict else 1
        node_sizes = [
            NODE_BASE_SIZE * (1.0 + 3.0 * (deg_dict[n] / max_deg))
            for n in G.nodes
        ]

        # TODO 7: Compute layout and draw graph
        pos = nx.spring_layout(G, seed=SEED)
        plt.figure(figsize=(10, 8))
        nx.draw_networkx_edges(
            G,
            pos,
            alpha=EDGE_ALPHA,
            width=0.5,
        )

        nx.draw_networkx_nodes(
            G,
            pos,
            node_color=node_colors,
            node_size=node_sizes,
            linewidths=0.5,
            edgecolors="black",
        )

        hub_genes = set(hubs_df["Gene"])
        labels = {n: n for n in G.nodes if n in hub_genes}
        nx.draw_networkx_labels(
            G,
            pos,
            labels=labels,
            font_size=8,
        )
        plt.axis("off")
        plt.tight_layout()

        # TODO 8: Save network figure
        plt.savefig(OUT_PNG, dpi=300, bbox_inches="tight")
        plt.close()

        # TODO 9: Save hub genes to CSV
        hubs_df.to_csv(OUT_HUBS, index=False)

        # print completion message
        print(f"Saved network figure to: {OUT_PNG}")
        print(f"Saved top-{TOPK_HUBS} hub genes to: {OUT_HUBS}")

        #bonus
        nx.write_edgelist(
            G,
            OUT_EDGES,
            delimiter=",",
            data=False 
        )
        print(f"Saved edge list to: {OUT_EDGES}")
