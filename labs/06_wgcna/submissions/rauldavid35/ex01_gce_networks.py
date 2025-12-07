"""
Exercițiu Gene Co-Expression Networks (GCEs) — Construirea rețelei și detectarea modulelor

Obiectiv:
- Să construiți o rețea de co-expresie dintr-o matrice de expresie RNA-Seq
- Să detectați module (comunități) de gene folosind un algoritm de tip Louvain (sau alternativ)

Instrucțiuni (în laborator):
1) Pregătire date
   - Descărcați și pregătiți matricea de expresie (ex: GSE115469) într-un CSV cu:
     * rânduri = gene (index), coloane = probe (sample IDs)
   - Salvați fișierul la: data/work/<handle>/lab06/expression_matrix.csv

2) Preprocesare
   - log2(x + 1)
   - filtrare gene cu varianță scăzută

3) Corelație → Adiacență
   - completați funcția `correlation_matrix`
   - funcția `adjacency_from_correlation` este deja implementată

4) Graf + Module
   - construiți graful cu NetworkX
   - detectați modulele (Louvain sau alternativă)
   - exportați mapping-ul gene → modul în submissions/<handle>/modules_<handle>.csv

Notă:
- Documentați în <github_handle>_notes.md: metrica de corelație, pragul, observații scurte.
"""

from __future__ import annotations
from pathlib import Path
from typing import Dict, Iterable

import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt


# --------------------------
# Config — completați după nevoie
# --------------------------
INPUT_CSV = Path("data/work/rauldavid35/lab06/expression_matrix.csv")
OUTPUT_DIR = Path("labs/06_wgcna/submissions/rauldavid35")
OUTPUT_CSV = OUTPUT_DIR / "modules_rauldavid35.csv"

CORR_METHOD = "spearman"   # TODO: "pearson" sau "spearman"
VARIANCE_THRESHOLD = 0.05   # prag pentru filtrare gene (pe date log2)
ADJ_THRESHOLD = 0.6    # prag pentru |cor| (ex: 0.6–0.8)
USE_ABS_CORR = True       # True => folosiți |cor| la prag
MAKE_UNDIRECTED = True     # rețelele de co-expresie sunt de obicei neorientate


def read_expression_matrix(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, index_col=0)
    return df


def log_and_filter(df: pd.DataFrame,
                   variance_threshold: float) -> pd.DataFrame:
    """
    Preprocesare:
    - aplică log2(x+1)
    - filtrează genele cu varianță scăzută
    """
    df_log = np.log2(df + 1)

    variances = df_log.var(axis=1)

    mask = variances >= variance_threshold
    df_filtered = df_log.loc[mask].copy()

    print(
        f"După log2(x+1) + filtrare varianță: "
        f"{df_filtered.shape[0]} gene păstrate din {df.shape[0]}."
    )
    return df_filtered


def correlation_matrix(df: pd.DataFrame,
                       method: str = "spearman",
                       use_abs: bool = True) -> pd.DataFrame:
    """
    TODO: calculați matricea de corelație între gene (rânduri).
    """
    corr = df.T.corr(method=method)

    if use_abs:
        corr = corr.abs()

    return corr


def adjacency_from_correlation(corr: pd.DataFrame,
                               threshold: float,
                               weighted: bool = False) -> pd.DataFrame:
    """
    Construiți matricea de adiacență din corelații.
    - binară: A_ij = 1 dacă corr_ij >= threshold, altfel 0
    - ponderată: A_ij = corr_ij dacă corr_ij >= threshold, altfel 0
    """
    if weighted:
        A = corr.where(corr >= threshold, other=0.0)
    else:
        A = (corr >= threshold).astype(int)

    np.fill_diagonal(A.values, 0)

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
    try:
        from networkx.algorithms.community import louvain_communities

        communities = louvain_communities(G, seed=42)
        method_used = "louvain"
    except Exception:
        from networkx.algorithms.community import greedy_modularity_communities

        communities = greedy_modularity_communities(G)
        method_used = "greedy_modularity"

    print(
        f"Comunități detectate cu metoda: {method_used}. "
        f"Număr de comunități: {len(communities)}"
    )

    mapping: Dict[str, int] = {}
    for module_id, community in enumerate(communities):
        for gene in community:
            mapping[str(gene)] = module_id

    return mapping


def save_modules_csv(mapping: Dict[str, int], out_csv: Path) -> None:
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    df_modules = (
        pd.DataFrame({"Gene": list(mapping.keys()), "Module": list(mapping.values())})
        .sort_values(["Module", "Gene"])
    )
    df_modules.to_csv(out_csv, index=False)

def visualize_network(G: nx.Graph, module_mapping: Dict[str, int], filename: Path):
    plt.figure(figsize=(12, 10))
    
    pos = nx.spring_layout(G, k=0.15, seed=42)
    
    node_colors = [module_mapping.get(str(n)) for n in G.nodes()]
    
    nx.draw_networkx_nodes(G, pos, node_color=node_colors, cmap=plt.cm.tab20, node_size=100)
    nx.draw_networkx_edges(G, pos, alpha=0.3)
    
    if G.number_of_nodes() < 100:
        nx.draw_networkx_labels(G, pos, font_size=8)
    
    plt.title(f"Network: {G.number_of_nodes()} nodes, {len(set(module_mapping.values()))} modules")
    plt.axis("off")
    plt.tight_layout()
    
    plt.savefig(filename, dpi=150) 
    print(f"Plot salvat cu succes la: {filename}")
    plt.close()


if __name__ == "__main__":
    expr_df = read_expression_matrix(INPUT_CSV)
    print(
        f"Am citit matricea de expresie din {INPUT_CSV} "
        f"cu dimensiunea: {expr_df.shape[0]} gene x {expr_df.shape[1]} probe."
    )

    expr_preproc = log_and_filter(expr_df, VARIANCE_THRESHOLD)

    corr = correlation_matrix(expr_preproc, method=CORR_METHOD, use_abs=USE_ABS_CORR)

    A = adjacency_from_correlation(corr, threshold=ADJ_THRESHOLD, weighted=False)

    G = graph_from_adjacency(A, undirected=MAKE_UNDIRECTED)

    print(f"Grafic creat cu {G.number_of_nodes()} noduri și {G.number_of_edges()} muchii.")

    gene_to_module = detect_modules_louvain_or_greedy(G)
    print(f"S-au detectat {len(set(gene_to_module.values()))} module.")
    
    if G.number_of_nodes() > 0:
        visualize_network(G, gene_to_module, OUTPUT_DIR) # <--- Trimitem calea fișierului
    else:
        print("Nu se poate genera plotul: Graful nu are noduri.")

    save_modules_csv(gene_to_module, OUTPUT_CSV)
    print(f"Am salvat mapping-ul gene→modul în: {OUTPUT_CSV}")