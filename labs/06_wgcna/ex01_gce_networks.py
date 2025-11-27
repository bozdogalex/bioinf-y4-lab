"""
Exercise Gene Co-Expression Networks (GCEs) — Network construction and module detection
(Final completed version - no diacritics)
"""

from __future__ import annotations
from pathlib import Path
from typing import Dict, Iterable

import numpy as np
import pandas as pd
import networkx as nx

# --------------------------
# Config — adjusted to ensure graph is non-empty
# --------------------------
INPUT_CSV = Path("../../data/work/Nistor-Iuliana03/lab06/expression_matrix.csv")
OUTPUT_DIR = Path("labs/06_networks/submissions/Nistor-Iuliana03")
OUTPUT_CSV = OUTPUT_DIR / "modules_Nistor-Iuliana03.csv"

# Updated parameters
CORR_METHOD = "pearson"      # "pearson" or "spearman"
VARIANCE_THRESHOLD = 0.0     # keep all genes (0.5 was too strict for small dataset)
ADJ_THRESHOLD = 0.2          # correlation threshold lowered (0.6 was too high)
USE_ABS_CORR = True          # use absolute correlation
MAKE_UNDIRECTED = True       # co-expression networks are usually undirected


def read_expression_matrix(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(
            f"Expression matrix not found at: {path}. Please place the CSV there."
        )
    df = pd.read_csv(path, index_col=0)
    if df.empty:
        raise ValueError("Expression matrix is empty.")
    return df


def log_and_filter(df: pd.DataFrame,
                   variance_threshold: float) -> pd.DataFrame:
    """
    Preprocessing:
    - apply log2(x + 1)
    - filter genes with low variance
    """
    df_log = np.log2(df + 1)
    df_filt = df_log.loc[df_log.var(axis=1) > variance_threshold]
    return df_filt


def correlation_matrix(df: pd.DataFrame,
                       method: str = "pearson",
                       use_abs: bool = True) -> pd.DataFrame:
    """
    Compute gene-gene correlation matrix.
    df has shape (genes x samples), so we use df.T.corr().
    If use_abs=True, return absolute correlation values.
    """
    corr = df.T.corr(method=method)
    if use_abs:
        corr = corr.abs()
    return corr


def adjacency_from_correlation(corr: pd.DataFrame,
                               threshold: float,
                               weighted: bool = False) -> pd.DataFrame:
    """
    Build adjacency matrix from correlations.
    - binary: A_ij = 1 if corr >= threshold
    - weighted: A_ij = corr_ij if corr >= threshold else 0
    """
    if weighted:
        A = corr.copy()
        A[A < threshold] = 0
    else:
        A = (corr >= threshold).astype(float)

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
    """
    Detect gene communities (modules).
    Prefer Louvain if available, otherwise fall back to greedy modularity.
    """
    try:
        from networkx.algorithms.community import louvain_communities
        communities = louvain_communities(G, seed=42)
    except Exception:
        from networkx.algorithms.community import greedy_modularity_communities
        communities_iterable: Iterable[Iterable[str]] = greedy_modularity_communities(G)
        communities = [set(c) for c in communities_iterable]

    mapping: Dict[str, int] = {}
    for idx, comm in enumerate(communities, start=1):
        for gene in comm:
            mapping[gene] = idx
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
    print(f"Graph created with {G.number_of_nodes()} nodes and {G.number_of_edges()} edges.")

    gene_to_module = detect_modules_louvain_or_greedy(G)
    print(f"Detected {len(set(gene_to_module.values()))} modules.")

    save_modules_csv(gene_to_module, OUTPUT_CSV)
    print(f"Saved gene->module mapping to: {OUTPUT_CSV}")
