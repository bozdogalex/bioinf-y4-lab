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


INPUT_CSV = Path("../../../../data/sample/expression_data.csv")
OUTPUT_DIR = Path("")
OUTPUT_CSV = OUTPUT_DIR / "modules_AlexTGoCreative.csv"

CORR_METHOD = "spearman"   # "pearson" sau "spearman"
VARIANCE_THRESHOLD = 0.1   # prag pentru filtrare gene
ADJ_THRESHOLD = 0.6     # prag pentru |cor| (ex: 0.6)
USE_ABS_CORR = True        # True => folosiți |cor| la prag
MAKE_UNDIRECTED = True      # rețelele de co-expresie sunt de obicei neorientate


def read_expression_matrix(path: Path) -> pd.DataFrame:
    """
    Citește matricea de expresie din CSV.
    Presupunem că prima coloană conține numele genelor (index).
    """
    df = pd.read_csv(path, index_col=0)
    print(f"Matricea citită: {df.shape[0]} gene × {df.shape[1]} probe")
    return df


def log_and_filter(df: pd.DataFrame,
                   variance_threshold: float) -> pd.DataFrame:
    """
    Preprocesare:
    - aplică log2(x+1)
    - filtrează genele cu varianță scăzută
    """
    # Log-transformare
    df_log = np.log2(df + 1)
    
    # Calculează varianța pe rânduri (gene)
    variances = df_log.var(axis=1)
    
    # Filtrează genele cu varianță scăzută
    mask = variances >= variance_threshold
    df_filtered = df_log[mask]
    
    print(f"După filtrare: {df_filtered.shape[0]} gene rămase (din {df_log.shape[0]})")
    return df_filtered

def correlation_matrix(df: pd.DataFrame,
                       method: str = "spearman",
                       use_abs: bool = True) -> pd.DataFrame:
    """
    Calculați matricea de corelație între gene (rânduri).
    method: 'pearson' sau 'spearman'
    use_abs: dacă True, returnează valorile absolute ale corelațiilor
    """
    # Calculăm corelația între rânduri (gene)
    # Transpunem pentru a avea genele pe coloane pentru .corr()
    corr = df.T.corr(method=method)
    
    if use_abs:
        corr = corr.abs()
    
    print(f"Matricea de corelație calculată: {corr.shape[0]}×{corr.shape[1]} (metoda: {method})")
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
        # Matricea ponderată: păstrăm valorile corelației peste prag
        A[A < threshold] = 0
    else:
        # Matricea binară: 1 dacă >= threshold, 0 altfel
        A = (A >= threshold).astype(int)
    
    # Eliminăm diagonal (auto-corelațiile)
    np.fill_diagonal(A.values, 0)
    
    num_edges = (A > 0).sum().sum() // 2  # împărțim la 2 pentru că e simetrică
    print(f"Matricea de adiacență: {num_edges} muchii peste pragul {threshold}")
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
    Detectați comunități (module) și întoarceți un dict gene -> modul_id.
    Variante:
      - încercați louvain_communities(G, seed=42) dacă e disponibil
      - altfel greedy_modularity_communities(G)
    """
    try:
        # Încercăm să folosim algoritmul Louvain
        from networkx.algorithms.community import louvain_communities
        communities = louvain_communities(G, seed=42)
        print("Folosim algoritmul Louvain pentru detectarea modulelor.")
    except ImportError:
        # Dacă Louvain nu e disponibil, folosim greedy modularity
        from networkx.algorithms.community import greedy_modularity_communities
        communities = greedy_modularity_communities(G)
        print("Folosim algoritmul Greedy Modularity pentru detectarea modulelor.")
    
    # Construim mapping gene -> modul_id
    gene_to_module = {}
    for module_id, community in enumerate(communities):
        for gene in community:
            gene_to_module[gene] = module_id
    
    return gene_to_module


def save_modules_csv(mapping: Dict[str, int], out_csv: Path) -> None:
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    df_modules = (
        pd.DataFrame({"Gene": list(mapping.keys()), "Module": list(mapping.values())})
        .sort_values(["Module", "Gene"])
    )
    df_modules.to_csv(out_csv, index=False)


if __name__ == "__main__":
    print("=" * 60)
    print("Gene Co-Expression Networks — Construirea rețelei și modulelor")
    print("=" * 60)
    
    # 1. Citim matricea de expresie
    print("\n[1] Citim matricea de expresie...")
    df_expr = read_expression_matrix(INPUT_CSV)
    
    # 2. Preprocesare: log2(x+1) și filtrare
    print("\n[2] Preprocesare: log-transformare și filtrare varianță...")
    df_filtered = log_and_filter(df_expr, VARIANCE_THRESHOLD)
    
    # 3. Calculăm matricea de corelație
    print("\n[3] Calculăm matricea de corelație...")
    corr = correlation_matrix(df_filtered, method=CORR_METHOD, use_abs=USE_ABS_CORR)
    
    # 4. Construim matricea de adiacență
    print("\n[4] Construim matricea de adiacență...")
    A = adjacency_from_correlation(corr, threshold=ADJ_THRESHOLD, weighted=False)
    
    # 5. Construim graful
    print("\n[5] Construim graful NetworkX...")
    G = graph_from_adjacency(A, undirected=MAKE_UNDIRECTED)
    print(f"Grafic creat cu {G.number_of_nodes()} noduri și {G.number_of_edges()} muchii.")
    
    # 6. Detectăm modulele
    print("\n[6] Detectăm modulele (community detection)...")
    gene_to_module = detect_modules_louvain_or_greedy(G)
    print(f"S-au detectat {len(set(gene_to_module.values()))} module.")
    
    # 7. Salvăm rezultatele
    print("\n[7] Salvăm mapping-ul gene→modul...")
    save_modules_csv(gene_to_module, OUTPUT_CSV)
    print(f"Am salvat mapping-ul gene→modul în: {OUTPUT_CSV}")
    
    print("\n" + "=" * 60)
    print("✓ Procesare completă!")
    print("=" * 60)
