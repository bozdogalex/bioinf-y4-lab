from __future__ import annotations
from pathlib import Path
from typing import Dict, Iterable, Optional

import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt

HANDLE = "MadalinaDenisa"

# Input files (cale corectă)
EXPR_CSV = Path(f"../06_wgcna/data/work/{HANDLE}/lab06/expression_matrix.csv")
MODULES_CSV = Path(f"../06_wgcna/labs/06_networks/submissions/{HANDLE}/modules_{HANDLE}.csv")
PRECOMPUTED_ADJ_CSV = None

CORR_METHOD = "spearman"
USE_ABS_CORR = True
ADJ_THRESHOLD = 0.6
WEIGHTED = False

SEED = 42
TOPK_HUBS = 10
NODE_BASE_SIZE = 60
EDGE_ALPHA = 0.15

OUT_DIR = Path(f"labs/07_networkviz/submissions/{HANDLE}")
OUT_PNG = OUT_DIR / f"network_{HANDLE}.png"
OUT_HUBS = OUT_DIR / f"hubs_{HANDLE}.csv"

def ensure_exists(path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Missing file: {path}")

def read_expression_matrix(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, index_col=0)
    return df

def read_modules_csv(path: Path) -> Dict[str, int]:
    df = pd.read_csv(path)
    return dict(zip(df["Gene"], df["Module"]))

def correlation_to_adjacency(expr: pd.DataFrame,
                             method: str,
                             use_abs: bool,
                             threshold: float,
                             weighted: bool) -> pd.DataFrame:
    corr = expr.corr(method=method)
    if use_abs:
        corr = corr.abs()
    A = corr.copy()
    A[A < threshold] = 0
    if not weighted:
        A[A >= threshold] = 1
    np.fill_diagonal(A.values, 0)
    return A

def graph_from_adjacency(A: pd.DataFrame) -> nx.Graph:
    G = nx.from_pandas_adjacency(A)
    G.remove_nodes_from(list(nx.isolates(G)))
    return G

def color_map_from_modules(nodes: Iterable[str], gene2module: Dict[str, int]) -> Dict[str, str]:
    cmap = plt.get_cmap("tab10")
    colors = {}
    for n in nodes:
        mod = gene2module.get(n, -1)
        colors[n] = cmap(mod % 10)
    return colors

def compute_hubs(G: nx.Graph, topk: int) -> pd.DataFrame:
    deg = dict(G.degree())
    df = pd.DataFrame({"Gene": list(deg.keys()), "Degree": list(deg.values())})
    return df.sort_values("Degree", ascending=False).head(topk)

if __name__ == "__main__":
    ensure_exists(EXPR_CSV)
    ensure_exists(MODULES_CSV)
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    expr = read_expression_matrix(EXPR_CSV)
    gene2module = read_modules_csv(MODULES_CSV)

    A = correlation_to_adjacency(expr, CORR_METHOD, USE_ABS_CORR, ADJ_THRESHOLD, WEIGHTED)

    G = graph_from_adjacency(A)

    node_colors_dict = color_map_from_modules(G.nodes(), gene2module)
    node_colors = [node_colors_dict[n] for n in G.nodes()]

    hubs_df = compute_hubs(G, TOPK_HUBS)
    hubs_set = set(hubs_df["Gene"])
    node_sizes = [NODE_BASE_SIZE * 3 if n in hubs_set else NODE_BASE_SIZE for n in G.nodes()]

    pos = nx.spring_layout(G, seed=SEED)
    plt.figure(figsize=(12, 10))

    nx.draw_networkx_edges(G, pos, alpha=EDGE_ALPHA)
    nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=node_sizes)

    labels = {n: n for n in G.nodes() if n in hubs_set}
    nx.draw_networkx_labels(G, pos, labels, font_size=8)

    plt.title(f"Gene Co-Expression Network — {HANDLE}")
    plt.axis("off")
    plt.savefig(OUT_PNG, dpi=300)

    hubs_df.to_csv(OUT_HUBS, index=False)
    print("Network visualization exported.")
