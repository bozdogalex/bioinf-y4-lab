from __future__ import annotations
from pathlib import Path
from typing import Dict, Iterable

import numpy as np
import pandas as pd
import networkx as nx


INPUT_CSV = Path("data/work/eeeevn0/lab06/expression_matrix1.csv")
OUTPUT_DIR = Path("labs/06_wgcna/submissions/eeeevn0")
OUTPUT_CSV = OUTPUT_DIR / "modules_eeeevn0.csv"

CORR_METHOD = "spearman"   
VARIANCE_THRESHOLD = 0.001  
ADJ_THRESHOLD = 0.6        
USE_ABS_CORR = True        
MAKE_UNDIRECTED = True   


def read_expression_matrix(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(
            f"Nu am găsit {path}. E nevoie de matricea de expresie la această locație."
        )
    df = pd.read_csv(path, index_col=0)
    if df.empty:
        raise ValueError("Matricea de expresie este goală.")
    return df


def log_and_filter(df: pd.DataFrame,
                   variance_threshold: float) -> pd.DataFrame:
    df_log = np.log2(df + 1)
    df_filt = df_log.loc[df_log.var(axis=1) > variance_threshold]
    return df_filt


def correlation_matrix(df: pd.DataFrame,
                       method: str = "spearman",
                       use_abs: bool = True) -> pd.DataFrame:
    corr = df.T.corr(method=method)
    if use_abs:
        corr = corr.abs()
    return corr


def adjacency_from_correlation(corr: pd.DataFrame,
                               threshold: float,
                               weighted: bool = False) -> pd.DataFrame:
    if weighted:
        A = corr.copy()
        A[A < threshold] = 0
    else:
        A = (corr >= threshold).astype(int)
    np.fill_diagonal(A.values, 0.0)
    return A


def graph_from_adjacency(A: pd.DataFrame,
                         undirected: bool = True) -> nx.Graph:
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
        communities = louvain_communities(G, seed=42)
    except Exception:
        from networkx.algorithms.community import greedy_modularity_communities
        communities_iterable: Iterable[Iterable[str]] = greedy_modularity_communities(G)
        communities = [set(c) for c in communities_iterable]

    mapping: Dict[str, int] = {}
    for midx, comm in enumerate(communities, start=1):
        for gene in comm:
            mapping[gene] = midx
    return mapping


def save_modules_csv(mapping: Dict[str, int], out_csv: Path) -> None:
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    df_modules = (
        pd.DataFrame({"Gene": list(mapping.keys()), "Module": list(mapping.values())})
        .sort_values(["Module", "Gene"])
    )
    df_modules.to_csv(out_csv, index=False)


if __name__ == "__main__":
    expr = read_expression_matrix(INPUT_CSV)
    expr_pp = log_and_filter(expr, variance_threshold=VARIANCE_THRESHOLD)

    corr = correlation_matrix(expr_pp, method=CORR_METHOD, use_abs=USE_ABS_CORR)
    adj = adjacency_from_correlation(corr, threshold=ADJ_THRESHOLD, weighted=False)

    G = graph_from_adjacency(adj, undirected=MAKE_UNDIRECTED)
    print(f"Grafic creat cu {G.number_of_nodes()} noduri și {G.number_of_edges()} muchii.")

    gene_to_module = detect_modules_louvain_or_greedy(G)
    print(f"S-au detectat {len(set(gene_to_module.values()))} module.")

    save_modules_csv(gene_to_module, OUTPUT_CSV)
    print(f"Am salvat mapping-ul gene→modul în: {OUTPUT_CSV}")
