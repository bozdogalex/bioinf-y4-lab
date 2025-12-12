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
from networkx.algorithms import community


# --------------------------
# Config — completați după nevoie
# --------------------------
INPUT_CSV =Path("data/work/ioanamhl/lab06/expression_matrix2.csv")
OUTPUT_DIR =Path("labs/06_wgcna/submissions/ioanamhl")
OUTPUT_CSV =OUTPUT_DIR / "modules_ioanamhl.csv"

CORR_METHOD = "spearman"   # TODO: "pearson" sau "spearman"
VARIANCE_THRESHOLD =0.5  # prag pentru filtrare gene
ADJ_THRESHOLD =0.6     # prag pentru |cor| (ex: 0.6)
USE_ABS_CORR =True        # True => folosiți |cor| la prag
MAKE_UNDIRECTED =True      # rețelele de co-expresie sunt de obicei neorientate


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
    keep = variances >= variance_threshold
    df_filt = df_log.loc[keep].copy()

    return df_filt

def correlation_matrix(df: pd.DataFrame,
                       method: str = "spearman",
                       use_abs: bool = True) -> pd.DataFrame:
    """
    TODO: calculați matricea de corelație între gene (rânduri).
    """
    # TODO: înlocuiți acest placeholder cu implementarea voastră
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
        A = corr.copy()
        A[A < threshold] = 0.0
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
    TODO: detectați comunități (module) și întoarceți un dict gene -> modul_id.
    Variante:
      - încercați louvain_communities(G, seed=42) dacă e disponibil
      - altfel greedy_modularity_communities(G)
    """
    try:
        # networkx >= 3.0
        comms = community.louvain_communities(G, seed=42)
    except AttributeError:
        # fallback: greedy modularity
        comms = community.greedy_modularity_communities(G)

    gene_to_module: Dict[str, int] = {}
    for module_id, nodes in enumerate(comms):
        for node in nodes:
            gene_to_module[str(node)] = module_id

    return gene_to_module


def save_modules_csv(mapping: Dict[str, int], out_csv: Path) -> None:
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    df_modules = (
        pd.DataFrame({"Gene": list(mapping.keys()), "Module": list(mapping.values())})
        .sort_values(["Module", "Gene"])
    )
    df_modules.to_csv(out_csv, index=False)


if __name__ == "__main__":
    expr = read_expression_matrix(INPUT_CSV)
    expr_proc = log_and_filter(expr, VARIANCE_THRESHOLD)

    corr = correlation_matrix(expr_proc, method=CORR_METHOD, use_abs=USE_ABS_CORR)
    A = adjacency_from_correlation(corr, threshold=ADJ_THRESHOLD, weighted=False)
    G = graph_from_adjacency(A, undirected=MAKE_UNDIRECTED)

    print(f"Grafic creat cu {G.number_of_nodes()} noduri și {G.number_of_edges()} muchii.")

    gene_to_module = detect_modules_louvain_or_greedy(G)
    print(f"S-au detectat {len(set(gene_to_module.values()))} module.")

    save_modules_csv(gene_to_module, OUTPUT_CSV)
    print(f"Am salvat mapping-ul gene→modul în: {OUTPUT_CSV}")
