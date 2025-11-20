"""
Exercițiu Gene Co-Expression Networks (GCEs) — Construirea rețelei și detectarea modulelor

Obiectiv:
- Să construiți o rețea de co-expresie dintr-o matrice de expresie RNA-Seq
- Să detectați module (comunități) de gene folosind un algoritm de tip Louvain (sau alternativ)
"""
#python -u labs/06_wgcna/submissions/florina-lucaciu/ex01_gce_networks.py

from __future__ import annotations
from pathlib import Path
from typing import Dict, Iterable

import numpy as np
import pandas as pd
import networkx as nx


# --------------------------
# Config — completați după nevoie
# --------------------------
INPUT_CSV = Path("data/work/florina-lucaciu/lab06/expression_matrix.csv")
OUTPUT_DIR = Path("labs/06_wgcna/submissions/florina-lucaciu")
OUTPUT_CSV = OUTPUT_DIR / "modules_florina-lucaciu.csv"

CORR_METHOD = "spearman"   # "pearson" sau "spearman"
VARIANCE_THRESHOLD = 0.5   # prag pentru filtrare gene
ADJ_THRESHOLD = 0.6        # prag pentru |cor| (ex: 0.6)
USE_ABS_CORR = True        # True => folosiți |cor| la prag
MAKE_UNDIRECTED = True     # rețelele de co-expresie sunt de obicei neorientate


def read_expression_matrix(path: Path) -> pd.DataFrame:
    """Citește matricea de expresie din CSV."""
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
    print(f"Filtrare: {len(df)} gene → {len(df_filt)} gene (varianță > {variance_threshold})")
    return df_filt


def correlation_matrix(df: pd.DataFrame,
                       method: str = "spearman",
                       use_abs: bool = True) -> pd.DataFrame:
    """
    Calculează matricea de corelație între gene (rânduri).
    
    Parameters:
    -----------
    df : pd.DataFrame
        Matricea de expresie (gene × probe)
    method : str
        Metoda de corelație: "pearson" sau "spearman"
    use_abs : bool
        Dacă True, returnează valoarea absolută a corelației
    
    Returns:
    --------
    pd.DataFrame
        Matricea de corelație (gene × gene)
    """
    # Transpunem pentru a calcula corelația între gene (rânduri)
    # df.T.corr() calculează corelația între coloane ale df.T, adică între rânduri ale df
    corr = df.T.corr(method=method)
    
    if use_abs:
        corr = corr.abs()
    
    return corr


def adjacency_from_correlation(corr: pd.DataFrame,
                               threshold: float,
                               weighted: bool = False) -> pd.DataFrame:
    """
    Construiește matricea de adiacență din corelații.
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
    """Construiește graful NetworkX din matricea de adiacență."""
    if undirected:
        G = nx.from_pandas_adjacency(A)
    else:
        G = nx.from_pandas_adjacency(A, create_using=nx.DiGraph)
    
    isolates = list(nx.isolates(G))
    if isolates:
        print(f"Eliminare {len(isolates)} noduri izolate")
        G.remove_nodes_from(isolates)
    
    return G


def detect_modules_louvain_or_greedy(G: nx.Graph) -> Dict[str, int]:
    """
    Detectează comunități (module) și returnează un dict gene -> modul_id.
    
    Încearcă mai întâi algoritmul Louvain, apoi cade pe greedy modularity
    dacă Louvain nu e disponibil.
    """
    try:
        # Încearcă Louvain (disponibil în NetworkX >= 2.8)
        from networkx.algorithms.community import louvain_communities
        communities_sets = louvain_communities(G, seed=42)
        print("Folosesc algoritmul Louvain pentru detectarea modulelor")
    except (ImportError, AttributeError):
        # Fallback pe greedy modularity
        from networkx.algorithms.community import greedy_modularity_communities
        print("Folosesc algoritmul Greedy Modularity pentru detectarea modulelor")
        communities_iterable: Iterable[Iterable[str]] = greedy_modularity_communities(G)
        communities_sets = [set(c) for c in communities_iterable]
    
    # Creează mapping-ul gene -> modul_id
    mapping: Dict[str, int] = {}
    for midx, comm in enumerate(communities_sets, start=1):
        for gene in comm:
            mapping[gene] = midx
    
    return mapping


def save_modules_csv(mapping: Dict[str, int], out_csv: Path) -> None:
    """Salvează mapping-ul gene -> modul într-un CSV."""
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    df_modules = (
        pd.DataFrame({"Gene": list(mapping.keys()), "Module": list(mapping.values())})
        .sort_values(["Module", "Gene"])
    )
    df_modules.to_csv(out_csv, index=False)


def analyze_modules(mapping: Dict[str, int]) -> None:
    """Afișează statistici despre modulele detectate."""
    module_sizes = pd.Series(mapping.values()).value_counts().sort_index()
    print("\n=== Statistici Module ===")
    print(f"Număr total de module: {len(module_sizes)}")
    print(f"Dimensiune medie modul: {module_sizes.mean():.1f} gene")
    print(f"Dimensiune mediană modul: {module_sizes.median():.1f} gene")
    print(f"\nTop 5 cele mai mari module:")
    print(module_sizes.head(5).to_string())


if __name__ == "__main__":
    print("=== Începe analiza rețelei de co-expresie genică ===\n")
    
    # 1. Citire date
    expr = read_expression_matrix(INPUT_CSV)
    print(f"Matrice citită: {expr.shape[0]} gene × {expr.shape[1]} probe\n")
    
    # 2. Preprocesare
    expr_pp = log_and_filter(expr, variance_threshold=VARIANCE_THRESHOLD)
    
    # 3. Calcul corelații și adiacență
    print(f"\nCalcul matricea de corelație ({CORR_METHOD})...")
    corr = correlation_matrix(expr_pp, method=CORR_METHOD, use_abs=USE_ABS_CORR)
    
    print(f"Construire matrice de adiacență (prag = {ADJ_THRESHOLD})...")
    adj = adjacency_from_correlation(corr, threshold=ADJ_THRESHOLD, weighted=False)
    
    # 4. Construire graf
    print(f"\nConstruire graf NetworkX...")
    G = graph_from_adjacency(adj, undirected=MAKE_UNDIRECTED)
    print(f"Graf creat cu {G.number_of_nodes()} noduri și {G.number_of_edges()} muchii.")
    
    # 5. Detectare module
    print(f"\nDetectare module...")
    gene_to_module = detect_modules_louvain_or_greedy(G)
    
    # 6. Analiză și salvare
    analyze_modules(gene_to_module)
    
    save_modules_csv(gene_to_module, OUTPUT_CSV)
    print(f"\n Am salvat mapping-ul gene→modul în: {OUTPUT_CSV}")