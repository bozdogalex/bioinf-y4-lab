"""
Exercițiu 8 — Vizualizarea rețelelor de co-expresie + gene hub

Obiectiv:
- Incărcați modulele detectate în Lab 6 și reconstruiți o rețea (din corelații)
- Vizualizați graful, colorând nodurile după modul
- Evidențiați genele hub (grad mare) și exportați figura (.png)

Intrări:
- Matricea de expresie folosită în Lab 6: data/work/<handle>/lab06/expression_matrix.csv
- Mapping gene→modul din Lab 6: labs/06_networks/submissions/<handle>/modules_<handle>.csv

Ieșiri:
- labs/07_network_viz/submissions/<handle>/network.png
- labs/07_network_viz/submissions/<handle>/hubs.csv  (opțional, listă gene hub)

Notă:
- Dacă aveți deja o matrice de adiacență salvată din Lab 6, o puteți încărca în loc să o reconstruiți.
- În acest exercițiu ne concentrăm pe VIZUALIZARE (nu refacem detectarea de module).
"""

from __future__ import annotations
from pathlib import Path
from typing import Dict, Iterable, Optional

import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt

HANDLE = "MariusJalba"

EXPR_CSV = Path(f"data/work/{HANDLE}/lab06/expression_matrix.csv")
MODULES_CSV = Path(f"labs/06_networks/submissions/{HANDLE}/modules_{HANDLE}.csv")

PRECOMPUTED_ADJ_CSV: Optional[Path] = None

CORR_METHOD = "spearman"
USE_ABS_CORR = True        
ADJ_THRESHOLD = 0.6        
WEIGHTED = False           

SEED = 42                 
TOPK_HUBS = 10             
NODE_BASE_SIZE = 60        
EDGE_ALPHA = 0.15          

# Ieșiri
OUT_DIR = Path(f"labs/07_network_viz/submissions/{HANDLE}")
OUT_PNG = OUT_DIR / f"network_{HANDLE}.png"
OUT_HUBS = OUT_DIR / f"hubs_{HANDLE}.csv"


def ensure_exists(path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Nu am găsit: {path}")


def read_expression_matrix(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, index_col=0)
    if df.empty:
        raise ValueError("Matricea de expresie nu contine date.")
    return df


def read_modules_csv(path: Path) -> Dict[str, int]:
    df = pd.read_csv(path)
    if not {"Gene", "Module"}.issubset(df.columns):
        raise ValueError("modules.csv nu contine coloanele 'Gene' și 'Module'.")
    return dict(zip(df["Gene"].astype(str), df["Module"].astype(int)))


def correlation_to_adjacency(expr: pd.DataFrame,
                             method: str,
                             use_abs: bool,
                             threshold: float,
                             weighted: bool) -> pd.DataFrame:
    correlation = expr.T.corr(method=method)
    if use_abs:
        correlation = correlation.abs()
    if weighted:
        Adjancency = corr.copy()
        Adjancency[Adjancency < threshold] = 0.0
    else:
        Adjancency = (correlation >= threshold).astype(float)
    np.fill_diagonal(Adjancency.values, 0.0)
    return Adjancency


def graph_from_adjacency(Adjancency: pd.DataFrame) -> nx.Graph:
    Graph = nx.from_pandas_adjacency(Adjancency)
    isolates = list(nx.isolates(Graph))
    if isolates:
        Graph.remove_nodes_from(isolates)
    return Graph


def color_map_from_modules(nodes: Iterable[str], gene2module: Dict[str, int]) -> Dict[str, str]:
    cmap = plt.get_cmap("tab10")
    colors: Dict[str, str] = {}
    for node in nodes:
        m = gene2module.get(node, 0)
        colors[node] = cmap((m - 1) % 10) if m > 0 else "#CCCCCC"
    return colors


def compute_hubs(Graph: nx.Graph, topk: int) -> pd.DataFrame:
    node_degrees = dict(Graph.degree())

    if Graph.number_of_nodes() <= 5000:
        betweenness = nx.betweenness_centrality(Graph)
    else:
        betweenness = {node: np.nan for node in Graph.nodes()}

    hub_table = (
        pd.DataFrame({
            "Gene": list(node_degrees.keys()),
            "Degree": list(node_degrees.values()),
            "Betweenness": [betweenness[node] for node in node_degrees.keys()]
        })
        .sort_values(by="Degree", ascending=False)
        .head(topk)
        .reset_index(drop=True)
    )

    return hub_table

if __name__ == "__main__":
    ensure_exists(EXPR_CSV)
    ensure_exists(MODULES_CSV)
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    expr = read_expression_matrix(EXPR_CSV)
    gene2module = read_modules_csv(MODULES_CSV)

    if PRECOMPUTED_ADJ_CSV is not None:
        ensure_exists(PRECOMPUTED_ADJ_CSV)
        Adjancency = pd.read_csv(PRECOMPUTED_ADJ_CSV, index_col=0)
        common = sorted(set(Adjancency.index) & set(gene2module.keys()))
        Adjancency = Adjancency.loc[common, common]
    else:
        Adjancency = correlation_to_adjacency(expr, CORR_METHOD, USE_ABS_CORR, ADJ_THRESHOLD, WEIGHTED)
        common = sorted(set(Adjancency.index) & set(gene2module.keys()))
        Adjancency = Adjancency.loc[common, common]

    Graph = graph_from_adjacency(Adjancency)
    print(f"Grafic: {Graph.number_of_nodes()} noduri, {Graph.number_of_edges()} muchii.")

    node_colors_map = color_map_from_modules(Graph.nodes(), gene2module)
    node_colors = [node_colors_map[node] for node in Graph.nodes()]

    hubs_df = compute_hubs(Graph, TOPK_HUBS)
    hubs_set = set(hubs_df["Gene"])
    node_sizes = [NODE_BASE_SIZE * (1.5 if node in hubs_set else 1.0) for node in Graph.nodes()]

    pos = nx.spring_layout(Graph, seed=SEED, k=None)
    plt.figure(figsize=(12, 10))
    nx.draw_networkx_edges(Graph, pos, alpha=EDGE_ALPHA, width=0.5)
    nx.draw_networkx_nodes(Graph, pos, node_color=node_colors, node_size=node_sizes, linewidths=0.0)

    hub_labels = {g: g for g in hubs_set if g in Graph.nodes()}
    nx.draw_networkx_labels(Graph, pos, labels=hub_labels, font_size=8)

    plt.axis("off")
    plt.tight_layout()
    plt.savefig(OUT_PNG, dpi=300)
    plt.close()
    print(f"Am salvat figura în: {OUT_PNG}")

    hubs_df.to_csv(OUT_HUBS, index=False)
    print(f"Am salvat hub genes în: {OUT_HUBS}")
