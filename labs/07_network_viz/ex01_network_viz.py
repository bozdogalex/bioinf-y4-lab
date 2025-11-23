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
HANDLE = "<handle>"

# Input files
EXPR_CSV = Path(f"data/work/{HANDLE}/lab06/expression_matrix.csv")
MODULES_CSV = Path(f"labs/06_networks/submissions/{HANDLE}/modules_{HANDLE}.csv")

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
OUT_DIR = Path(f"labs/07_networkviz/submissions/{HANDLE}")
OUT_PNG = OUT_DIR / f"network_{HANDLE}.png"
OUT_HUBS = OUT_DIR / f"hubs_{HANDLE}.csv"


# --------------------------
# Utils
# --------------------------
def ensure_exists(path: Path) -> None:
    """TODO: check that a file exists."""
    pass


def read_expression_matrix(path: Path) -> pd.DataFrame:
    """
    TODO:
    - read CSV
    - set index to gene names
    """
    return pd.DataFrame()  # placeholder


def read_modules_csv(path: Path) -> Dict[str, int]:
    """
    TODO:
    - read CSV with columns: Gene, Module
    - return dict: gene -> module_id
    """
    return {}  # placeholder


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
    return pd.DataFrame()  # placeholder


def graph_from_adjacency(A: pd.DataFrame) -> nx.Graph:
    """
    TODO:
    - convert adjacency DataFrame to NetworkX graph
    - remove isolated nodes
    """
    return nx.Graph()  # placeholder


def color_map_from_modules(nodes: Iterable[str], gene2module: Dict[str, int]) -> Dict[str, str]:
    """
    TODO:
    - assign a color to each node based on its module
    - use matplotlib 'tab10' or another palette
    """
    return {}  # placeholder


def compute_hubs(G: nx.Graph, topk: int) -> pd.DataFrame:
    """
    TODO:
    - compute degree for every node
    - (optional) compute betweenness centrality
    - return top-k genes
    """
    return pd.DataFrame()  # placeholder


# --------------------------
# Main
# --------------------------
if __name__ == "__main__":
    # TODO 1: Verify input files exist
    # ensure_exists(EXPR_CSV)
    # ensure_exists(MODULES_CSV)
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # TODO 2: Load expression matrix and module mapping
    # expr = read_expression_matrix(EXPR_CSV)
    # gene2module = read_modules_csv(MODULES_CSV)

    # TODO 3: Load or reconstruct adjacency
    # if PRECOMPUTED_ADJ_CSV exists:
    #       load adjacency and filter by module genes
    # else:
    #       compute adjacency from correlations

    # TODO 4: Build graph
    # G = graph_from_adjacency(A)
    # print info about number of nodes and edges

    # TODO 5: Compute colors by module
    # node_colors = [...]

    # TODO 6: Compute hub genes
    # hubs_df = compute_hubs(G, TOPK_HUBS)
    # node_sizes = [...]

    # TODO 7: Compute layout and draw graph
    # - draw edges
    # - draw nodes (colored)
    # - draw labels only for hubs

    # TODO 8: Save network figure
    # plt.savefig(OUT_PNG)

    # TODO 9: Save hub genes to CSV
    # hubs_df.to_csv(OUT_HUBS, index=False)

    # print completion message
