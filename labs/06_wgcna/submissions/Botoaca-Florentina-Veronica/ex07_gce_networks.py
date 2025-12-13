from __future__ import annotations
from pathlib import Path
from typing import Dict, Iterable, List, Set

import numpy as np
import pandas as pd
import networkx as nx


HANDLE = "Botoaca-Florentina-Veronica"
INPUT_CSV = Path(f"/workspaces/bioinf-y4-lab/labs/01_intro&databases/data/work/Botoaca-Florentina-Veronica/lab01/expression_matrix.csv")
OUTPUT_DIR = Path(f"/workspaces/bioinf-y4-lab/labs/06_wgcna/submissions/Botoaca-Florentina-Veronica")
OUTPUT_CSV = OUTPUT_DIR / f"modules_{HANDLE}.csv"

CORR_METHOD = "spearman"
VARIANCE_PERCENTILE = 50  # Păstrează top 50% gene după varianță
ADJ_THRESHOLD = 0.6
USE_ABS_CORR = True
MAKE_UNDIRECTED = True


def read_expression_matrix(path: Path) -> pd.DataFrame:
    if not path.exists():
        # Încearcă path alternativ
        alt_path = Path(f"data/work/{HANDLE}/lab06/expression_matrix.csv")
        if alt_path.exists():
            print(f"Folosesc calea alternativă: {alt_path}")
            path = alt_path
        else:
            raise FileNotFoundError(
                f"Nu am găsit {path}. Generează mai întâi datele!"
            )
    
    df = pd.read_csv(path, index_col=0)
    print(f"Matrice încărcată: {df.shape[0]} gene, {df.shape[1]} probe")
    print(f"Valori min/max: {df.values.min():.2f} / {df.values.max():.2f}")
    
    # Verifică dacă datele arată OK
    if df.isnull().any().any():
        print("⚠️ Avertisment: Există valori NaN în date!")
        df = df.fillna(df.mean())
    
    return df


def log_and_filter(df: pd.DataFrame,
                   variance_percentile: float = 50) -> pd.DataFrame:
    """
    Preprocesare cu prag dinamic bazat pe percentile
    """
    print("\n=== PREPROCESARE DATES ===")
    print(f"Date intrare: {df.shape[0]} gene")
    
    # Log transform
    df_log = np.log2(df + 1)
    print("✓ Transformare log2(x+1) aplicată")
    
    # Calculează varianțe
    variances = df_log.var(axis=1)
    
    # Determină pragul
    threshold = np.percentile(variances, variance_percentile)
    
    print(f"\nStatistici varianță:")
    print(f"  Min:    {variances.min():.8f}")
    print(f"  Medie:  {variances.mean():.8f}")
    print(f"  Max:    {variances.max():.8f}")
    print(f"  Prag ({variance_percentile}%): {threshold:.8f}")
    
    # Afișează distribuția
    percentiles = [0, 25, 50, 75, 95, 100]
    print(f"\nDistribuție varianțe:")
    for p in percentiles:
        val = np.percentile(variances, p)
        print(f"  {p:3d}%: {val:.8f}")
    
    # Filtrare
    mask = variances > threshold
    n_passing = mask.sum()
    
    print(f"\n✓ Gene care trec pragul: {n_passing} / {len(variances)}")
    
    if n_passing == 0:
        print("\n❌ CRITIC: Nicio genă nu trece pragul!")
        print("Folosesc prag de urgență (percentila 0)...")
        threshold = np.percentile(variances, 0)  # prag minim
        mask = variances > threshold
        n_passing = mask.sum()
        print(f"Gene cu prag de urgență: {n_passing}")
    
    df_filt = df_log.loc[mask]
    
    print(f"✓ Gene rămase după filtrare: {df_filt.shape[0]}")
    return df_filt


def correlation_matrix(df: pd.DataFrame,
                       method: str = "spearman",
                       use_abs: bool = True) -> pd.DataFrame:
    """
    Calculează matricea de corelație între gene
    """
    print("\n=== CALCUL CORELAȚII ===")
    print(f"Metodă: {method}")
    
    # Transpunem pentru corelație între gene
    corr = df.T.corr(method=method)
    
    if use_abs:
        corr = corr.abs()
        print("✓ Folosesc valorile absolute ale corelației")
    
    print(f"✓ Matrice corelație: {corr.shape[0]}×{corr.shape[1]}")
    print(f"  Valori: [{corr.values.min():.3f}, {corr.values.max():.3f}]")
    
    return corr


def adjacency_from_correlation(corr: pd.DataFrame,
                               threshold: float,
                               weighted: bool = False) -> pd.DataFrame:
    """
    Construiește matricea de adiacență
    """
    print(f"\n=== CONSTRUIRE ADIACENȚĂ ===")
    print(f"Prag pentru muchie: {threshold}")
    
    if weighted:
        A = corr.copy()
        A[A < threshold] = 0
        print("✓ Matrice ponderată")
    else:
        A = (corr >= threshold).astype(int)
        print("✓ Matrice binară")
    
    # Diagonala la 0
    np.fill_diagonal(A.values, 0.0)
    
    # Numără muchiile
    if not weighted:
        n_edges = (A > 0).sum().sum() / 2
        density = n_edges / (len(A) * (len(A) - 1) / 2)
        print(f"✓ Muchii: {int(n_edges)}")
        print(f"✓ Densitate graf: {density:.6f}")
    
    return A


def graph_from_adjacency(A: pd.DataFrame,
                         undirected: bool = True) -> nx.Graph:
    print(f"\n=== CONSTRUIRE GRAF ===")
    
    if undirected:
        G = nx.from_pandas_adjacency(A)
        print("✓ Graf neorientat creat")
    else:
        G = nx.from_pandas_adjacency(A, create_using=nx.DiGraph)
        print("✓ Graf orientat creat")
    
    # Elimină noduri izolate
    isolates = list(nx.isolates(G))
    if isolates:
        print(f"✂️ Elimin {len(isolates)} noduri izolate")
        G.remove_nodes_from(isolates)
    
    print(f"✓ Noduri finale: {G.number_of_nodes()}")
    print(f"✓ Muchii finale: {G.number_of_edges()}")
    
    if G.number_of_nodes() == 0:
        print("⚠️ ATENȚIE: Graful este gol!")
        print("  Motiv posibil: prag ADJ_THRESHOLD prea mare")
        print(f"  Sugestie: scade ADJ_THRESHOLD la 0.4 sau 0.3")
    
    return G


def detect_modules_louvain_or_greedy(G: nx.Graph) -> Dict[str, int]:
    """
    Detectează module (comunități)
    """
    print(f"\n=== DETECTARE MODULE ===")
    
    if G.number_of_nodes() == 0:
        print("❌ Graful este gol - nu pot detecta module")
        return {}
    
    try:
        from networkx.algorithms.community import louvain_communities
        communities_list = list(louvain_communities(G, seed=42, resolution=1.0))
        communities = [set(c) for c in communities_list]
        print("✓ Algoritm: Louvain")
    except Exception as e:
        print(f"⚠️ Louvain indisponibil: {e}")
        from networkx.algorithms.community import greedy_modularity_communities
        communities_iterable = greedy_modularity_communities(G)
        communities = [set(c) for c in communities_iterable]
        print("✓ Algoritm: Greedy Modularity (fallback)")
    
    mapping: Dict[str, int] = {}
    for midx, comm in enumerate(communities, start=1):
        for gene in comm:
            mapping[gene] = midx
    
    print(f"✓ Module detectate: {len(communities)}")
    
    if len(communities) > 0:
        sizes = [len(c) for c in communities]
        print(f"✓ Dimensiuni module:")
        print(f"  Min: {min(sizes)} gene")
        print(f"  Max: {max(sizes)} gene")
        print(f"  Medie: {np.mean(sizes):.1f} gene")
    
    return mapping


def save_modules_csv(mapping: Dict[str, int], out_csv: Path) -> None:
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    
    if not mapping:
        print("\n❌ Niciun modul detectat - creăm fișier gol")
        pd.DataFrame(columns=["Gene", "Module"]).to_csv(out_csv, index=False)
    else:
        df_modules = (
            pd.DataFrame({"Gene": list(mapping.keys()), 
                         "Module": list(mapping.values())})
            .sort_values(["Module", "Gene"])
        )
        df_modules.to_csv(out_csv, index=False)
        print(f"✓ Fișier salvat: {out_csv}")
        print(f"✓ Total gene în module: {len(df_modules)}")


def main():
    print("=" * 60)
    print("LABORATOR 6: GENE CO-EXPRESSION NETWORKS")
    print(f"Handle: {HANDLE}")
    print("=" * 60)
    
    try:
        # 1. Încarcă date
        expr = read_expression_matrix(INPUT_CSV)
        
        # 2. Preprocesare
        expr_pp = log_and_filter(expr, variance_percentile=VARIANCE_PERCENTILE)
        
        # Verifică dacă avem suficiente gene
        if expr_pp.shape[0] < 10:
            print("\n⚠️ PREA PUȚINE GENE! Ajustez parametrii...")
            # Reîncepe cu prag mai mic
            expr_pp = log_and_filter(expr, variance_percentile=20)
        
        if expr_pp.shape[0] == 0:
            print("❌ EROARE: Nicio genă rămasă după filtrare")
            print("   Verifică datele de intrare sau scade VARIANCE_PERCENTILE")
            return
        
        # 3. Corelații
        corr = correlation_matrix(expr_pp, method=CORR_METHOD, use_abs=USE_ABS_CORR)
        
        # 4. Adiacență
        adj = adjacency_from_correlation(corr, threshold=ADJ_THRESHOLD, weighted=False)
        
        # 5. Graf
        G = graph_from_adjacency(adj, undirected=MAKE_UNDIRECTED)
        
        if G.number_of_nodes() == 0:
            print("\n❌ Graful este gol - încerc prag mai mic pentru adiacență")
            # Reîncepe cu prag mai mic pentru muchii
            adj = adjacency_from_correlation(corr, threshold=0.4, weighted=False)
            G = graph_from_adjacency(adj, undirected=MAKE_UNDIRECTED)
        
        # 6. Module
        gene_to_module = detect_modules_louvain_or_greedy(G)
        
        # 7. Salvează
        save_modules_csv(gene_to_module, OUTPUT_CSV)
        
        print("\n" + "=" * 60)
        print("REZUMAT FINAL:")
        print(f"- Gene inițiale: {expr.shape[0]}")
        print(f"- Gene după filtrare: {expr_pp.shape[0]}")
        print(f"- Noduri în graf: {G.number_of_nodes()}")
        print(f"- Muchii în graf: {G.number_of_edges()}")
        print(f"- Module detectate: {len(set(gene_to_module.values()))}")
        print("=" * 60)
        
    except Exception as e:
        print(f"\n❌ EROARE: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    main()