from __future__ import annotations
from pathlib import Path
from typing import Dict, Iterable, Optional

import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt


HANDLE = "Razvann19"

EXPR_CSV = Path(f"labs/07_network_viz/submissions/Razvann19/expression_matrix.csv")
MODULES_CSV = Path(f"labs/07_network_viz/submissions/Razvann19/modules_Razvann19.csv")
PRECOMPUTED_ADJ_CSV: Optional[Path] = None

CORR_METHOD = "spearman"  
USE_ABS_CORR = True       
ADJ_THRESHOLD = 0.6        
WEIGHTED = False           

SEED = 42
TOPK_HUBS = 10
NODE_BASE_SIZE = 60
EDGE_ALPHA = 0.3

OUT_DIR = Path(f"labs/07_network_viz/submissions/Razvann19")
OUT_PNG = OUT_DIR / f"network_Razvann19.png"
OUT_HUBS = OUT_DIR / f"hubs_Razvann19.csv"


def ensure_exists(path: Path) -> None:
    """Verifică dacă fișierul există."""
    if not path.exists():
        raise FileNotFoundError(f"Nu am găsit fișierul: {path}")


def read_expression_matrix(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, index_col=0)
    if df.empty:
        raise ValueError("Matricea de expresie este goală.")
    return df


def read_modules_csv(path: Path) -> Dict[str, int]:
    """Verifică dacă fișierul cu module este citit corect"""
    df = pd.read_csv(path)
    if not {"Gene", "Module"}.issubset(df.columns):
        raise ValueError("Fișierul modules.csv trebuie să conțină coloanele: Gene, Module")
    return dict(zip(df["Gene"].astype(str), df["Module"].astype(int)))


def correlation_to_adjacency(expr: pd.DataFrame,
                             method: str,
                             use_abs: bool,
                             threshold: float,
                             weighted: bool) -> pd.DataFrame:
    """Reconstruiește matricea de adiacență din corelații."""
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
    """Construiește un graf din matricea de adiacență."""
    G = nx.from_pandas_adjacency(A)
    isolates = list(nx.isolates(G))
    if isolates:
        G.remove_nodes_from(isolates)
    return G


def color_map_from_modules(nodes: Iterable[str], gene2module: Dict[str, int]) -> Dict[str, str]:
    cmap = plt.get_cmap("tab10")
    colors: Dict[str, str] = {}
    for n in nodes:
        m = gene2module.get(n, 0)
        colors[n] = cmap((m - 1) % 10) if m > 0 else "#CCCCCC"
    return colors


def compute_hubs(G: nx.Graph, topk: int) -> pd.DataFrame:
    """Calcularea hub-urilor (noduri cu grad mare)."""
    deg = dict(G.degree())
    btw = nx.betweenness_centrality(G, normalized=True, seed=SEED) if G.number_of_nodes() <= 5000 else {n: np.nan for n in G.nodes()}
    hubs = (
        pd.DataFrame(
            {
                "Gene": list(deg.keys()),
                "Degree": list(deg.values()),
                "Betweenness": [btw.get(n, np.nan) for n in deg.keys()],
            }
        )
        .sort_values(["Degree", "Betweenness"], ascending=False)
        .head(topk)
        .reset_index(drop=True)
    )
    return hubs



if __name__ == "__main__":
    ensure_exists(EXPR_CSV)
    ensure_exists(MODULES_CSV)
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    expr = read_expression_matrix(EXPR_CSV)
    gene2module = read_modules_csv(MODULES_CSV)


    if PRECOMPUTED_ADJ_CSV is not None:
        ensure_exists(PRECOMPUTED_ADJ_CSV)
        A = pd.read_csv(PRECOMPUTED_ADJ_CSV, index_col=0)
        common = sorted(set(A.index) & set(gene2module.keys()))
        A = A.loc[common, common]
    else:
        A = correlation_to_adjacency(expr, CORR_METHOD, USE_ABS_CORR, ADJ_THRESHOLD, WEIGHTED)
        common = sorted(set(A.index) & set(gene2module.keys()))
        A = A.loc[common, common]

    G = graph_from_adjacency(A)
    print(f"Grafic: {G.number_of_nodes()} noduri, {G.number_of_edges()} muchii.")

    node_colors_map = color_map_from_modules(G.nodes(), gene2module)
    node_colors = [node_colors_map[n] for n in G.nodes()]

    hubs_df = compute_hubs(G, TOPK_HUBS)
    hubs_set = set(hubs_df["Gene"])
    node_sizes = [NODE_BASE_SIZE * (1.5 if n in hubs_set else 1.0) for n in G.nodes()]

    pos = nx.spring_layout(G, seed=SEED, k=0.1, iterations=200)  # Ajustat pentru a distribui mai bine nodurile
    plt.figure(figsize=(15, 15))
    nx.draw_networkx_edges(G, pos, alpha=EDGE_ALPHA, width=0.7)
    nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=node_sizes, linewidths=0.0)
    
    hub_labels = {g: g for g in hubs_set if g in G.nodes()}
    nx.draw_networkx_labels(G, pos, labels=hub_labels, font_size=8)

    plt.axis("off")
    plt.tight_layout()
    plt.savefig(OUT_PNG, dpi=300)
    plt.close()
    print(f"Grafic: {OUT_PNG}")

    hubs_df.to_csv(OUT_HUBS, index=False)
    print(f"Gene: {OUT_HUBS}")
