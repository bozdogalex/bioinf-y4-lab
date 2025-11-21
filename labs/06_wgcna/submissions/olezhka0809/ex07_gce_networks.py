#!/usr/bin/env python3

from __future__ import annotations
from pathlib import Path
from typing import Dict

import numpy as np
import pandas as pd
import networkx as nx


# --------------------------
# Config — adaptează doar HANDLE / praguri dacă vrei
# --------------------------
HANDLE = "olezhka0809"

INPUT_CSV = Path(f"/workspaces/bioinf-y4-lab/labs/06_wgcna/submissions/olezhka0809/expression_matrix.csv")
OUTPUT_DIR = Path(f"/workspaces/bioinf-y4-lab/labs/06_wgcna/submissions/{HANDLE}")
OUTPUT_CSV = OUTPUT_DIR / f"modules_{HANDLE}.csv"

CORR_METHOD = "spearman"     # "pearson" sau "spearman"
VARIANCE_THRESHOLD = 0.5     # prag pentru filtrare gene (varianță minimă)
ADJ_THRESHOLD = 0.6          # prag pentru |cor|
USE_ABS_CORR = True          # True => folosim |cor| la prag
MAKE_UNDIRECTED = True       # GCE = rețea neorientată


# --------------------------
# 1) Citire și preprocesare
# --------------------------

def read_expression_matrix(path: Path) -> pd.DataFrame:
    """
    Citește matricea de expresie.
    Presupune CSV cu:
      - index: gene
      - coloane: sample IDs
    """
    if not path.exists():
        raise FileNotFoundError(f"Nu am găsit fișierul de expresie: {path}")
    df = pd.read_csv(path, index_col=0)
    print(f"[INFO] Matrice expresie: {df.shape[0]} gene × {df.shape[1]} probe")
    return df


def log_and_filter(df: pd.DataFrame,
                   variance_threshold: float) -> pd.DataFrame:
    """
    Preprocesare:
    - aplică log2(x+1)
    - filtrează genele cu varianță scăzută
    """
    # log2(x + 1)
    df_log = np.log2(df + 1.0)

    # varianță pe gene (pe rând, across samples)
    gene_var = df_log.var(axis=1)
    mask = gene_var >= variance_threshold
    df_filt = df_log.loc[mask].copy()

    print(f"[INFO] După log2(x+1) și filtrare varianță >= {variance_threshold}: "
          f"{df_filt.shape[0]} gene păstrate din {df.shape[0]}")
    return df_filt


# --------------------------
# 2) Corelație & adiacență
# --------------------------

def correlation_matrix(df: pd.DataFrame,
                       method: str = "spearman",
                       use_abs: bool = True) -> pd.DataFrame:
    """
    Calculează matricea de corelație între gene.
    - df: rânduri = gene, coloane = probe
    - method: "pearson" sau "spearman"
    - use_abs: dacă True, întoarce |cor|
    """
    # Corr între gene => corelăm pe rânduri => df.T.corr(...)
    corr = df.T.corr(method=method)
    if use_abs:
        corr = corr.abs()

    # Asigură diagonală bine definită (nu contează în rețea)
    np.fill_diagonal(corr.values, 1.0)

    print(f"[INFO] Matrice de corelație ({method}, abs={use_abs}) calculată: "
          f"{corr.shape[0]} × {corr.shape[1]}")
    return corr


def adjacency_from_correlation(corr: pd.DataFrame,
                               threshold: float,
                               weighted: bool = False) -> pd.DataFrame:
    """
    Construiți matricea de adiacență din corelații.
    - binară: A_ij = 1 dacă corr_ij >= threshold, altfel 0
    - ponderată: A_ij = corr_ij dacă corr_ij >= threshold, altfel 0
    """
    A = corr.copy()

    if weighted:
        A.values[A.values < threshold] = 0.0
    else:
        A = (A >= threshold).astype(float)

    # scoatem self-loops explicit
    np.fill_diagonal(A.values, 0.0)

    # mic rezumat
    num_edges = (A.values > 0).sum() // 2  # fiind simetric
    print(f"[INFO] Matrice adiacență (threshold={threshold}, weighted={weighted}): "
          f"{A.shape[0]} noduri, ~{num_edges} muchii")
    return A


# --------------------------
# 3) Graph & módulos
# --------------------------

def graph_from_adjacency(A: pd.DataFrame,
                         undirected: bool = True) -> nx.Graph:
    if undirected:
        G = nx.from_pandas_adjacency(A)
    else:
        G = nx.from_pandas_adjacency(A, create_using=nx.DiGraph)

    isolates = list(nx.isolates(G))
    if isolates:
        print(f"[INFO] Elimin {len(isolates)} noduri izolate (fără muchii).")
        G.remove_nodes_from(isolates)

    print(f"[INFO] Graf final: {G.number_of_nodes()} noduri, {G.number_of_edges()} muchii.")
    return G


def detect_modules_louvain_or_greedy(G: nx.Graph) -> Dict[str, int]:
    """
    Detectează comunități (module) și întoarce un dict gene -> modul_id.
    Variante:
      - încearcă louvain_communities(G, seed=42) dacă e disponibil
      - altfel greedy_modularity_communities(G)
    """
    try:
        # NetworkX >= 2.8 + extra community
        from networkx.algorithms.community import louvain_communities
        print("[INFO] Folosesc Louvain pentru detectarea modulelor.")
        communities = louvain_communities(G, seed=42)
    except Exception:
        from networkx.algorithms.community import greedy_modularity_communities
        print("[INFO] Louvain nu e disponibil, folosesc greedy_modularity_communities.")
        communities = list(greedy_modularity_communities(G))

    # Construim mapping gene -> modul_id
    gene_to_module: Dict[str, int] = {}
    for mod_id, comm in enumerate(communities):
        for node in comm:
            gene_to_module[str(node)] = mod_id

    print(f"[INFO] Detectate {len(communities)} module.")
    return gene_to_module


def save_modules_csv(mapping: Dict[str, int], out_csv: Path) -> None:
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    df_modules = (
        pd.DataFrame({"Gene": list(mapping.keys()), "Module": list(mapping.values())})
        .sort_values(["Module", "Gene"])
    )
    df_modules.to_csv(out_csv, index=False)
    print(f"[OK] Am salvat mapping-ul gene→modul în: {out_csv}")


# --------------------------
# MAIN
# --------------------------

if __name__ == "__main__":
    print("\n==============================================")
    print(" GENE CO-EXPRESSION NETWORKS (GCE) - LAB 06  ")
    print("==============================================\n")

    # 1) Citește matricea de expresie
    expr = read_expression_matrix(INPUT_CSV)

    # 2) Preprocesare: log2(x+1) + filtrare varianță mică
    expr_proc = log_and_filter(expr, VARIANCE_THRESHOLD)

    # 3) Corelație între gene
    corr = correlation_matrix(
        expr_proc,
        method=CORR_METHOD,
        use_abs=USE_ABS_CORR,
    )

    # 4) Adiacență din corelație
    #    (eu am pus weighted=True, dar poți schimba în False dacă vrei rețea binară)
    A = adjacency_from_correlation(
        corr,
        threshold=ADJ_THRESHOLD,
        weighted=True,
    )

    # 5) Graf de co-expresie
    G = graph_from_adjacency(A, undirected=MAKE_UNDIRECTED)

    # 6) Module (comunități)
    gene_to_module = detect_modules_louvain_or_greedy(G)

    # 7) Salvare CSV
    save_modules_csv(gene_to_module, OUTPUT_CSV)

    print("\n[INFO] Rezumat specific pentru notes.md:")
    print(f" - metrica: {CORR_METHOD}")
    print(f" - use_abs_corr: {USE_ABS_CORR}")
    print(f" - VARIANCE_THRESHOLD: {VARIANCE_THRESHOLD}")
    print(f" - ADJ_THRESHOLD: {ADJ_THRESHOLD}")
    print(f" - număr module: {len(set(gene_to_module.values()))}")
