from __future__ import annotations
from pathlib import Path
from typing import Dict, Iterable

import numpy as np
import pandas as pd
import networkx as nx

HANDLE = "Ana-Maria-Bojan" 
INPUT_CSV = Path(f"data/work/{HANDLE}/lab06/expression_matrix.csv")
OUTPUT_DIR = Path(f"labs/06_wgcna/submissions/{HANDLE}")
OUTPUT_CSV = OUTPUT_DIR / f"modules_{HANDLE}.csv"

CORR_METHOD = "spearman"    # Metoda de corelație
VARIANCE_THRESHOLD = 0.5    # Prag pentru filtrare gene
ADJ_THRESHOLD = 0.6         # Prag pentru |cor| 
USE_ABS_CORR = True         # Folosiți |cor| la prag
MAKE_UNDIRECTED = True      # Rețelele de co-expresie sunt neorientate




def read_expression_matrix(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(
            f"Nu am găsit {path}. Puneți matricea de expresie la această locație."
        )
    df = pd.read_csv(path, index_col=0)
    if df.empty:
        raise ValueError("Matricea de expresie este goală.")
    return df


def log_and_filter(df: pd.DataFrame,
                     variance_threshold: float) -> pd.DataFrame:
    """
    Preprocesare: aplică log2(x+1) și filtrează genele cu varianță scăzută
    """
    df_log = np.log2(df + 1)
    df_filt = df_log.loc[df_log.var(axis=1) > variance_threshold]
    return df_filt


def correlation_matrix(df: pd.DataFrame,
                       method: str = "spearman",
                       use_abs: bool = True) -> pd.DataFrame:
    """
    CALCULEAZA MATRICEA DE CORELAȚIE ÎNTRE GENE.
    """
    # df.T.corr() calculează corelația între coloanele lui df.T, adică între rândurile lui df (genele).
    corr = df.T.corr(method=method) 
    
    if use_abs:
        corr = corr.abs()
        
    return corr


def adjacency_from_correlation(corr: pd.DataFrame,
                               threshold: float,
                               weighted: bool = False) -> pd.DataFrame:
    """
    Construiți matricea de adiacență din corelații.
    """
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
        # print(f"Eliminat {len(isolates)} gene izolate.")
        G.remove_nodes_from(isolates)
    return G


def detect_modules_louvain_or_greedy(G: nx.Graph) -> Dict[str, int]:
    """
    DETECTEAZA COMUNITĂȚI (MODULE) ȘI ÎNTOARCE DICT GENE -> MODUL_ID.
    """
    try:
        # Tentativa Louvain (necesita networkx > 2.8)
        from networkx.algorithms.community import louvain_communities
        # Folosim 'weight' in caz ca s-a ales retea ponderata (weighted=True)
        communities_iterable: Iterable[Iterable[str]] = louvain_communities(G, seed=42, weight='weight')
        print("Algoritm: Louvain")
    except Exception:
        # Fallback pe Greedy Modularity (disponibil universal)
        from networkx.algorithms.community import greedy_modularity_communities
        communities_iterable: Iterable[Iterable[str]] = greedy_modularity_communities(G)
        print("Algoritm: Greedy Modularity (Fallback)")
    
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
    # 1) Pregatire date (citire)
    print(f"1. Citire date din: {INPUT_CSV}")
    expr = read_expression_matrix(INPUT_CSV)
    
    # 2) Preprocesare (log2 si filtrare)
    expr_pp = log_and_filter(expr, variance_threshold=VARIANCE_THRESHOLD)
    print(f"2. Preprocesare: {expr_pp.shape[0]} gene rămase după filtrare.")

    # 3) Corelație -> Adiacență
    corr = correlation_matrix(expr_pp, method=CORR_METHOD, use_abs=USE_ABS_CORR)
    adj = adjacency_from_correlation(corr, threshold=ADJ_THRESHOLD, weighted=False)
    
    # 4) Graf + Module
    G = graph_from_adjacency(adj, undirected=MAKE_UNDIRECTED)
    print(f"3. Grafic creat cu {G.number_of_nodes()} noduri și {G.number_of_edges()} muchii.")

    gene_to_module = detect_modules_louvain_or_greedy(G)
    print(f"4. S-au detectat {len(set(gene_to_module.values()))} module.")

    save_modules_csv(gene_to_module, OUTPUT_CSV)
    print(f"5. Am salvat mapping-ul gene→modul în: {OUTPUT_CSV}")