"""
Assignment: Gene Co-Expression Networks, Visualization & Diseasome
Handle: AlexTGoCreative
"""

from __future__ import annotations
from pathlib import Path
from typing import Dict
import gzip
import urllib.request
import io

import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

HANDLE = "AlexTGoCreative"
BASE_DIR = Path(__file__).resolve().parents[4]

# Folosim datele create pentru lab06
EXPR_CSV = BASE_DIR / "data" / "work" / HANDLE / "lab06" / "expression_matrix.csv"

# Output direct în directorul curent (suntem deja în submission/AlexTGoCreative)
OUT_DIR = Path(__file__).resolve().parent
MODULES_CSV = OUT_DIR / f"modules_tp53_{HANDLE}.csv"
NETWORK_PNG = OUT_DIR / f"network_tp53_{HANDLE}.png"
HUBS_CSV = OUT_DIR / f"hubs_tp53_{HANDLE}.csv"
ENRICHMENT_CSV = OUT_DIR / f"enrichment_{HANDLE}.csv"

CORR_METHOD = "spearman"
VARIANCE_THRESHOLD = 0.01  # Prag foarte mic pentru a păstra toate genele
TOP_GENES = 500
ADJ_THRESHOLD = 0.6  # Prag mai mic pentru a păstra mai multe conexiuni
USE_ABS_CORR = True
TOPK_HUBS = 10
SEED = 42


# =============================================================================
# TASK 1: DATE SI PREPROCESARE 
# =============================================================================
def read_expression_matrix(path: Path) -> pd.DataFrame:
    """Citeste matricea de expresie din CSV."""
    print(f"   Citire: {path}")
    df = pd.read_csv(path, index_col=0)
    return df


def log_and_filter(df: pd.DataFrame, 
                   variance_threshold: float,
                   top_n: int = None) -> pd.DataFrame:
    """
    Preprocesare:
    - log2(x+1)
    - filtrare gene cu varianta scazuta
    - selectare top N gene dupa varianta
    """
    # log2(x + 1)
    df_log = np.log2(df + 1)
    
    # Calculeaza varianta pentru fiecare gena
    gene_variance = df_log.var(axis=1)
    
    # Filtreaza genele cu varianta sub prag
    df_filtered = df_log[gene_variance >= variance_threshold]
    print(f"   Gene dupa filtrare varianta: {len(df_filtered)} din {len(df_log)}")
    
    # Selecteaza top N gene dupa varianta
    if top_n is not None and len(df_filtered) > top_n:
        gene_variance_filtered = df_filtered.var(axis=1)
        top_genes = gene_variance_filtered.nlargest(top_n).index
        df_filtered = df_filtered.loc[top_genes]
        print(f"   Selectate top {top_n} gene dupa varianta")
    
    return df_filtered


# =============================================================================
# TASK 2: RETEA SI MODULE 
# =============================================================================
def correlation_matrix(df: pd.DataFrame, 
                       method: str = "spearman",
                       use_abs: bool = True) -> pd.DataFrame:
    """Calculeaza matricea de corelatie intre gene."""
    corr = df.T.corr(method=method)
    
    if use_abs:
        corr = corr.abs()
    
    # Seteaza diagonala la 0 (fara self-loops)
    np.fill_diagonal(corr.values, 0)
    
    return corr


def adjacency_from_correlation(corr: pd.DataFrame,
                               threshold: float) -> pd.DataFrame:
    """Construieste matricea de adiacenta binara din corelatii."""
    adj = (corr >= threshold).astype(int)
    return adj


def graph_from_adjacency(A: pd.DataFrame) -> nx.Graph:
    """Converteste matricea de adiacenta in graf NetworkX."""
    G = nx.from_pandas_adjacency(A)
    # Elimina nodurile izolate
    isolates = list(nx.isolates(G))
    if isolates:
        G.remove_nodes_from(isolates)
    return G


def detect_modules_louvain(G: nx.Graph) -> Dict[str, int]:
    """Detecteaza comunitati (module) folosind algoritmul Louvain."""
    gene_to_module: Dict[str, int] = {}
    
    try:
        communities = nx.community.louvain_communities(G, seed=SEED)
        print("   Algoritm: Louvain")
    except AttributeError:
        communities = list(nx.community.greedy_modularity_communities(G))
        print("   Algoritm: Greedy Modularity (fallback)")
    
    for module_id, community in enumerate(communities):
        for gene in community:
            gene_to_module[gene] = module_id
    
    return gene_to_module


def save_modules_csv(mapping: Dict[str, int], out_csv: Path) -> None:
    """Salveaza mapping-ul gene -> modul in CSV."""
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame({
        "Gene": list(mapping.keys()),
        "Module": list(mapping.values())
    }).sort_values(["Module", "Gene"])
    df.to_csv(out_csv, index=False)
    print(f"   Salvat: {out_csv.name}")


# =============================================================================
# TASK 3: VIZUALIZARE SI HUB GENES
# =============================================================================
def compute_hubs(G: nx.Graph, topk: int) -> pd.DataFrame:
    """Calculeaza hub genes bazat pe degree si betweenness centrality."""
    degrees = dict(G.degree())
    betweenness = nx.betweenness_centrality(G)
    
    df = pd.DataFrame({
        'Gene': list(degrees.keys()),
        'Degree': list(degrees.values()),
        'Betweenness': [betweenness[node] for node in degrees.keys()]
    })
    
    df = df.sort_values('Degree', ascending=False).head(topk)
    df = df.reset_index(drop=True)
    
    return df


def visualize_network(G: nx.Graph, 
                      gene2module: Dict[str, int],
                      hubs_df: pd.DataFrame,
                      output_path: Path) -> None:
    """Vizualizeaza reteaua cu noduri colorate dupa modul."""
    
    plt.figure(figsize=(16, 12))
    
    # Layout
    pos = nx.spring_layout(G, seed=SEED, k=2/np.sqrt(G.number_of_nodes()))
    
    # Culori pentru module
    n_modules = len(set(gene2module.values()))
    cmap = plt.colormaps.get_cmap('tab10')
    
    node_colors = []
    for node in G.nodes():
        module_id = gene2module.get(node, -1)
        if module_id >= 0:
            node_colors.append(cmap(module_id % 10))
        else:
            node_colors.append('lightgray')
    
    # Dimensiuni noduri bazate pe degree
    degrees = dict(G.degree())
    max_degree = max(degrees.values()) if degrees else 1
    node_sizes = [100 + (degrees[n] / max_degree) * 400 for n in G.nodes()]
    
    # Hub genes
    hub_genes = set(hubs_df['Gene'].tolist())
    
    # Desenare muchii
    nx.draw_networkx_edges(G, pos, alpha=0.15, edge_color='gray', width=0.5)
    
    # Desenare noduri
    nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=node_sizes, alpha=0.8)
    
    # Etichete doar pentru hub-uri
    hub_labels = {n: n for n in G.nodes() if n in hub_genes}
    nx.draw_networkx_labels(G, pos, labels=hub_labels, font_size=7, font_weight='bold')
    
    # Legenda module
    handles = []
    for i in range(n_modules):
        handles.append(plt.Line2D([0], [0], marker='o', color='w', 
                                   markerfacecolor=cmap(i % 10), markersize=10,
                                   label=f'Module {i}'))
    plt.legend(handles=handles, loc='upper left', fontsize=8)
    
    plt.title(f"Gene Co-Expression Network - TP53 ({HANDLE})\n"
              f"{G.number_of_nodes()} nodes, {G.number_of_edges()} edges, {n_modules} modules",
              fontsize=12)
    plt.axis('off')
    plt.tight_layout()
    
    output_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path, dpi=150, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"   Salvat: {output_path.name}")

def get_module_genes(gene2module: Dict[str, int], module_id: int) -> list:
    """Returneaza genele dintr-un modul specific."""
    return [gene for gene, mod in gene2module.items() if mod == module_id]


def enrichment_analysis_gprofiler(gene_list: list, output_csv: Path) -> pd.DataFrame:
    """
    Analiza de imbogatire folosind g:Profiler API.
    Returneaza rezultatele GO/KEGG.
    """
    try:
        from gprofiler import GProfiler
        
        gp = GProfiler(return_dataframe=True)
        results = gp.profile(organism='hsapiens', query=gene_list)
        
        if not results.empty:
            # Filtram rezultatele relevante
            results_filtered = results[results['source'].isin(['GO:BP', 'GO:MF', 'GO:CC', 'KEGG'])]
            results_filtered = results_filtered.sort_values('p_value').head(20)
            
            output_csv.parent.mkdir(parents=True, exist_ok=True)
            results_filtered.to_csv(output_csv, index=False)
            print(f"   Salvat: {output_csv.name}")
            return results_filtered
        else:
            print("   Nu s-au gasit rezultate de imbogatire.")
            return pd.DataFrame()
            
    except ImportError:
        print("   NOTA: gprofiler-official nu este instalat.")
        print("   Pentru analiza de imbogatire, folositi manual:")
        print("   https://biit.cs.ut.ee/gprofiler/gost")
        print(f"   Gene pentru analiza: {', '.join(gene_list[:20])}...")
        
        # Salvam lista de gene pentru analiza manuala
        genes_file = output_csv.parent / f"genes_for_enrichment_{HANDLE}.txt"
        with open(genes_file, 'w') as f:
            f.write('\n'.join(gene_list))
        print(f"   Lista gene salvata in: {genes_file.name}")
        
        return pd.DataFrame()

if __name__ == "__main__":
    print("=" * 70)
    print("ASSIGNMENT: Gene Co-Expression Networks, Visualization & Diseasome")
    print(f"Handle: {HANDLE}")
    print("=" * 70)
    
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    
    # =========================================================================
    # TASK 1: Date si preprocesare 
    # =========================================================================
    print("\n[TASK 1] Date si preprocesare")
    print("-" * 40)
    
    expr_df = read_expression_matrix(EXPR_CSV)
    print(f"   Dimensiuni initiale: {expr_df.shape[0]} gene x {expr_df.shape[1]} probe")
    
    expr_filtered = log_and_filter(expr_df, VARIANCE_THRESHOLD, TOP_GENES)
    print(f"   Dimensiuni dupa preprocesare: {expr_filtered.shape}")
    
    # Verificare: daca nu au ramas gene dupa filtrare, reducem criteriile
    if len(expr_filtered) == 0:
        print(f"\n   ⚠ ATENȚIE: Nicio genă nu a trecut filtrul de varianță!")
        print(f"   Datele sample sunt prea mici. Reduc pragul la 0.0...")
        expr_filtered = log_and_filter(expr_df, 0.0, TOP_GENES)
        print(f"   Dimensiuni după refiltrare: {expr_filtered.shape}")
    
    # =========================================================================
    # TASK 2: Retea si module 
    # =========================================================================
    print("\n[TASK 2] Retea si module")
    print("-" * 40)
    
    print("   Calculare matrice corelatie...")
    corr = correlation_matrix(expr_filtered, CORR_METHOD, USE_ABS_CORR)
    
    print(f"   Construire adiacenta (prag={ADJ_THRESHOLD})...")
    adj = adjacency_from_correlation(corr, ADJ_THRESHOLD)
    
    print("   Construire graf...")
    G = graph_from_adjacency(adj)
    print(f"   Graf: {G.number_of_nodes()} noduri, {G.number_of_edges()} muchii")
    
    # Adaugam atribut degree pentru fiecare nod
    for node in G.nodes():
        G.nodes[node]['degree'] = G.degree(node)

    # Exportam reteaua in format GraphML pentru Cytoscape (pentru Bonus)
    nx.write_graphml(G, OUT_DIR / "tp53_network.graphml")
    edges = nx.to_pandas_edgelist(G)
    edges.to_csv(OUT_DIR / "tp53_edges.csv", index=False)
    print(f"   Exportat pentru Cytoscape: tp53_network.graphml, tp53_edges.csv")

    
    print("   Detectare module (Louvain)...")
    gene2module = detect_modules_louvain(G)
    n_modules = len(set(gene2module.values()))
    print(f"   Module detectate: {n_modules}")
    
    save_modules_csv(gene2module, MODULES_CSV)
    
    # =========================================================================
    # TASK 3: Vizualizare si hub genes
    # =========================================================================
    print("\n[TASK 3] Vizualizare si hub genes")
    print("-" * 40)
    
    print(f"   Calculare hub genes (top {TOPK_HUBS})...")
    hubs_df = compute_hubs(G, TOPK_HUBS)
    
    print("   Hub genes:")
    for _, row in hubs_df.iterrows():
        module_id = gene2module.get(row['Gene'], -1)
        print(f"      {row['Gene']}: degree={row['Degree']}, "
              f"betweenness={row['Betweenness']:.4f}, module={module_id}")
    
    hubs_df.to_csv(HUBS_CSV, index=False)
    print(f"   Salvat: {HUBS_CSV.name}")
    
    print("   Vizualizare retea...")
    visualize_network(G, gene2module, hubs_df, NETWORK_PNG)
    
    # =========================================================================
    # TASK 4: Interpretare biologica si diseasome
    # =========================================================================
    print("\n[TASK 4] Interpretare biologica si diseasome")
    print("-" * 40)
    
    # Alegem modulul cu cele mai multe hub-uri sau cel mai mare
    module_sizes = {}
    for gene, mod in gene2module.items():
        module_sizes[mod] = module_sizes.get(mod, 0) + 1
    
    largest_module = max(module_sizes, key=module_sizes.get)
    print(f"   Analizam modulul {largest_module} ({module_sizes[largest_module]} gene)")
    
    module_genes = get_module_genes(gene2module, largest_module)
    
    print("   Analiza de imbogatire GO/KEGG...")
    enrichment_df = enrichment_analysis_gprofiler(module_genes, ENRICHMENT_CSV)
    
    # =========================================================================
    # Rezumat final
    # =========================================================================
    print("\n" + "=" * 70)
    print("✓ ASSIGNMENT COMPLET - Toate task-urile rezolvate!")
    print("=" * 70)
    print("\nFișiere generate (livrabile):")
    print(f"  1. {MODULES_CSV.name} (Task 2: gene → module)")
    print(f"  2. {NETWORK_PNG.name} (Task 3: vizualizare rețea)")
    print(f"  3. {HUBS_CSV.name} (Task 3: hub genes)")
    print(f"  4. assigment.py (codul folosit)")
    if ENRICHMENT_CSV.exists():
        print(f"  5. {ENRICHMENT_CSV.name} (Task 4: analiză îmbogățire)")
    print(f"\nFișiere bonus (pentru Cytoscape/Gephi):")
    print(f"  - tp53_network.graphml")
    print(f"  - tp53_edges.csv")
    print(f"\nRămâne de făcut:")
    print(f"  - report_{HANDLE}.pdf (Task 4: raport max 2 pagini)")
    print(f"  - Pentru bonus: screenshot Cytoscape/Gephi + observații")
    print("=" * 70)


