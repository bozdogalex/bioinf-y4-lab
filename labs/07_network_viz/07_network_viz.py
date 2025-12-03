"""
Exercise 8 — Network visualization + hub genes
"""

from __future__ import annotations
from pathlib import Path
from typing import Dict, Iterable, Optional

import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt

# --------------------------
# CONFIG — updated for your handle + Lab 6 params
# --------------------------
HANDLE = "Nistor-Iuliana03"

# Inputs
EXPR_CSV = Path(f"../../data/work/{HANDLE}/lab06/expression_matrix.csv")
MODULES_CSV = Path(f"../06_networks/submissions/{HANDLE}/modules_{HANDLE}.csv")

# No precomputed adjacency
PRECOMPUTED_ADJ_CSV: Optional[Path] = None

# Parameters matching Lab 6
CORR_METHOD = "pearson"
USE_ABS_CORR = True
ADJ_THRESHOLD = 0.2
WEIGHTED = False

# Visualization
SEED = 42
TOPK_HUBS = 3
NODE_BASE_SIZE = 200
EDGE_ALPHA = 0.25

# Outputs
OUT_DIR = Path(f"submissions/{HANDLE}")
OUT_PNG = OUT_DIR / f"network_{HANDLE}.png"
OUT_HUBS = OUT_DIR / f"hubs_{HANDLE}.csv"


# --------------------------
# Utils
# --------------------------
def ensure_exists(path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Missing file: {path}")


def read_expression_matrix(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, index_col=0)
    if df.empty:
        raise ValueError("Expression matrix is empty.")
    return df


def read_modules_csv(path: Path) -> Dict[str, int]:
    df = pd.read_csv(path)
    return dict(zip(df["Gene"].astype(str), df["Module"].astype(int)))


def correlation_to_adjacency(expr: pd.DataFrame,
                             method: str,
                             use_abs: bool,
                             threshold: float,
                             weighted: bool) -> pd.DataFrame:
    corr = expr.T.corr(method=method)
    if use_abs:
        corr = corr.abs()

    if weighted:
        A = corr.copy()
        A[A < threshold] = 0.0
    else:
        A = (corr >= threshold).astype(float)

    np.fill_diagonal(A.values, 0.0)
    return A


def graph_from_adjacency(A: pd.DataFrame) -> nx.Graph:
    G = nx.from_pandas_adjacency(A)
    isolates = list(nx.isolates(G))
    if isolates:
        G.remove_nodes_from(isolates)
    return G


def color_map_from_modules(nodes: Iterable[str], gene2module: Dict[str, int]):
    cmap = plt.get_cmap("tab10")
    colors = {}
    for n in nodes:
        m = gene2module.get(n, 0)
        colors[n] = cmap((m - 1) % 10) if m > 0 else "#CCCCCC"
    return colors


def compute_hubs(G: nx.Graph, topk: int):
    deg = dict(G.degree())
    hubs = sorted(deg.items(), key=lambda x: x[1], reverse=True)[:topk]
    df = pd.DataFrame(hubs, columns=["Gene", "Degree"])
    return df


# --------------------------
# MAIN
# --------------------------
if __name__ == "__main__":
    ensure_exists(EXPR_CSV)
    ensure_exists(MODULES_CSV)
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    expr = read_expression_matrix(EXPR_CSV)
    gene2module = read_modules_csv(MODULES_CSV)

    if PRECOMPUTED_ADJ_CSV is not None:
        A = pd.read_csv(PRECOMPUTED_ADJ_CSV, index_col=0)
    else:
        A = correlation_to_adjacency(expr, CORR_METHOD, USE_ABS_CORR, ADJ_THRESHOLD, WEIGHTED)

    common = sorted(set(A.index) & set(gene2module.keys()))
    A = A.loc[common, common]

    G = graph_from_adjacency(A)
    print(f"Graph: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")

    node_colors_map = color_map_from_modules(G.nodes(), gene2module)
    node_colors = [node_colors_map[n] for n in G.nodes()]

    hubs_df = compute_hubs(G, TOPK_HUBS)
    hubs = set(hubs_df["Gene"])

    node_sizes = [NODE_BASE_SIZE * (1.5 if n in hubs else 1.0) for n in G.nodes()]

    pos = nx.spring_layout(G, seed=SEED)

    plt.figure(figsize=(10, 8))
    nx.draw_networkx_edges(G, pos, alpha=EDGE_ALPHA)
    nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=node_sizes)
    nx.draw_networkx_labels(G, pos, {h: h for h in hubs}, font_size=8)

    plt.axis("off")
    plt.tight_layout()
    plt.savefig(OUT_PNG, dpi=300)
    plt.close()

    hubs_df.to_csv(OUT_HUBS, index=False)

    print(f"Saved figure to: {OUT_PNG}")
    print(f"Saved hubs to: {OUT_HUBS}")
