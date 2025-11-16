from __future__ import annotations
from pathlib import Path
from typing import Dict, Iterable

import numpy as np
import pandas as pd
import networkx as nx



INPUT_CSV = Path("data/work/StanaAndrei/lab06/gene_expr_matrix.csv")
OUTPUT_DIR = Path("labs/06_networks/submissions/StanaAndrei")
OUTPUT_CSV = OUTPUT_DIR / "modules_StanaAndrei.csv"

CORR_METHOD = "spearman"   # "pearson" sau "spearman"
VARIANCE_THRESHOLD = 0.5   # prag pentru filtrare gene
ADJ_THRESHOLD = 0.6        # prag pentru |cor| (ex: 0.6)
USE_ABS_CORR = True        # True => folosiți |cor| la prag
MAKE_UNDIRECTED = True     # rețelele de co-expresie sunt de obicei neorientate


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
    Preprocesare:
    - aplică log2(x+1)
    - filtrează genele cu varianță scăzută
    """
    df_log = np.log2(df + 1)
    df_filt = df_log.loc[df_log.var(axis=1) > variance_threshold]
    return df_filt


def correlation_matrix(df: pd.DataFrame,
                       method: str = "spearman",
                       use_abs: bool = True) -> pd.DataFrame:
    """
    TODO: calculați matricea de corelație între gene (rânduri).
    Hint:
      - df este (gene x probe); pentru corelație între gene, folosiți df.T.corr(method=...)
      - dacă use_abs=True, întoarceți |cor|
    """
    # df (gene x probe) -> df.T (probe x gene)
    # .corr() calculează corelația între coloane (gene)
    corr_mat = df.T.corr(method=method)
    
    if use_abs:
        return corr_mat.abs()
    return corr_mat


def adjacency_from_correlation(corr: pd.DataFrame,
                               threshold: float,
                               weighted: bool = False) -> pd.DataFrame:
    """
    Construiți matricea de adiacență din corelații.
    - binară: A_ij = 1 dacă corr_ij >= threshold, altfel 0
    - ponderată: A_ij = corr_ij dacă corr_ij >= threshold, altfel 0
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
        G.remove_nodes_from(isolates)
    return G


def detect_modules_louvain_or_greedy(G: nx.Graph) -> Dict[str, int]:
    """
    TODO: detectați comunități (module) și întoarceți un dict gene -> modul_id.
    Variante:
      - încercați louvain_communities(G, seed=42) dacă e disponibil
      - altfel greedy_modularity_communities(G)
    """
    # Schelet cu fallback pe greedy_modularity_communities:
    try:
        # Varianta 1: Louvain (preferată)
        from networkx.algorithms.community import louvain_communities
        communities_iterable = louvain_communities(G, seed=42)
        # networkx >= 3.0 întoarce seturi, versiunile mai vechi frozensets/liste
        communities = [set(c) for c in communities_iterable]

    except ImportError:
        # Varianta 2: Fallback la Greedy (dacă Louvain lipsește)
        print("Atenție: 'louvain_communities' nu a fost găsit. Se folosește 'greedy_modularity_communities'.")
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
    print(f"Se citesc datele din: {INPUT_CSV}")
    expr = read_expression_matrix(INPUT_CSV)
    print(f"Matrice inițială: {expr.shape[0]} gene, {expr.shape[1]} probe.")

    expr_pp = log_and_filter(expr, variance_threshold=VARIANCE_THRESHOLD)
    print(f"Matrice filtrată: {expr_pp.shape[0]} gene, {expr_pp.shape[1]} probe.")

    print(f"Se calculează corelația ({CORR_METHOD}, |corr|={USE_ABS_CORR})...")
    corr = correlation_matrix(expr_pp, method=CORR_METHOD, use_abs=USE_ABS_CORR)
    
    print(f"Se construiește adiacența (prag={ADJ_THRESHOLD})...")
    adj = adjacency_from_correlation(corr, threshold=ADJ_THRESHOLD, weighted=False)

    G = graph_from_adjacency(adj, undirected=MAKE_UNDIRECTED)
    print(f"Grafic creat cu {G.number_of_nodes()} noduri și {G.number_of_edges()} muchii.")

    print("Se detectează modulele...")
    gene_to_module = detect_modules_louvain_or_greedy(G)
    
    num_modules = len(set(gene_to_module.values()))
    if num_modules > 0:
        print(f"S-au detectat {num_modules} module.")
    else:
        print("Nu s-au detectat module (e posibil ca graful să fie gol).")

    if num_modules > 0:
        save_modules_csv(gene_to_module, OUTPUT_CSV)
        print(f"Am salvat mapping-ul gene→modul în: {OUTPUT_CSV}")