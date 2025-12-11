from __future__ import annotations
from pathlib import Path
from typing import Dict, Iterable

import numpy as np
import pandas as pd
import networkx as nx

INPUT_CSV = Path("data/work/Razvann19/lab06/GSE115469_Data.csv")
OUTPUT_DIR = Path("labs/06_wgcna/submissions/Razvann19")
OUTPUT_CSV = OUTPUT_DIR / "modules_Razvann19.csv"

CORR_METHOD = "pearson"
VARIANCE_THRESHOLD = 0.5
ADJ_THRESHOLD = 0.6
USE_ABS_CORR = True
MAKE_UNDIRECTED = True

MAX_ROWS = 8000       # câte rânduri citesc din CSV
N_TOP_GENES = 4000    # câte gene păstrez pentru rețea


def read_expression_matrix(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Nu am găsit {path}.")
    df = pd.read_csv(path, index_col=0, nrows=MAX_ROWS)
    if df.empty:
        raise ValueError("Matricea de expresie este goală.")
    df = df.astype(np.float32)
    return df


def log_and_filter(df: pd.DataFrame, variance_threshold: float) -> pd.DataFrame:
    df_log = np.log2(df + 1.0).astype(np.float32)
    var = df_log.var(axis=1)
    top_genes = var.sort_values(ascending=False).head(N_TOP_GENES).index
    df_filt = df_log.loc[top_genes]
    return df_filt


def correlation_matrix(df: pd.DataFrame, method: str = "pearson", use_abs: bool = True) -> pd.DataFrame:
    corr = df.T.corr(method=method)
    corr = corr.fillna(0.0)
    if use_abs:
        corr = corr.abs()
    np.fill_diagonal(corr.values, 1.0)
    return corr


def adjacency_from_correlation(corr: pd.DataFrame, threshold: float, weighted: bool = False) -> pd.DataFrame:
    if weighted:
        A = corr.copy()
        A[A < threshold] = 0
    else:
        A = (corr >= threshold).astype(int)
    np.fill_diagonal(A.values, 0.0)
    return A


def graph_from_adjacency(A: pd.DataFrame, undirected: bool = True) -> nx.Graph:
    if undirected:
        G = nx.from_pandas_adjacency(A)
    else:
        G = nx.from_pandas_adjacency(A, create_using=nx.DiGraph)
    isolates = list(nx.isolates(G))
    if isolates:
        G.remove_nodes_from(isolates)
    return G


def detect_modules_louvain_or_greedy(G: nx.Graph) -> Dict[str, int]:
    try:
        from networkx.algorithms.community import louvain_communities
        communities_iterable = louvain_communities(G, seed=42)
        communities = [set(c) for c in communities_iterable]
    except Exception:
        from networkx.algorithms.community import greedy_modularity_communities
        communities_iterable = greedy_modularity_communities(G)
        communities = [set(c) for c in communities_iterable]

    mapping: Dict[str, int] = {}
    for midx, comm in enumerate(communities, start=1):
        for gene in comm:
            mapping[gene] = midx
    return mapping


def save_modules_csv(mapping: Dict[str, int], out_csv: Path) -> None:
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    df_modules = pd.DataFrame({"Gene": list(mapping.keys()), "Module": list(mapping.values())})
    df_modules = df_modules.sort_values(["Module", "Gene"])
    df_modules.to_csv(out_csv, index=False)


if __name__ == "__main__":
    expr = read_expression_matrix(INPUT_CSV)
    expr_pp = log_and_filter(expr, variance_threshold=VARIANCE_THRESHOLD)

    print("Shape după filtrare:", expr_pp.shape)

    corr = correlation_matrix(expr_pp, method=CORR_METHOD, use_abs=USE_ABS_CORR)
    adj = adjacency_from_correlation(corr, threshold=ADJ_THRESHOLD, weighted=False)

    G = graph_from_adjacency(adj, undirected=MAKE_UNDIRECTED)
    print(f"Grafic creat cu {G.number_of_nodes()} noduri și {G.number_of_edges()} muchii.")

    gene_to_module = detect_modules_louvain_or_greedy(G)
    print(f"S-au detectat {len(set(gene_to_module.values()))} module.")

    save_modules_csv(gene_to_module, OUTPUT_CSV)
    print(f"Am salvat mapping-ul gene→modul în: {OUTPUT_CSV}")
