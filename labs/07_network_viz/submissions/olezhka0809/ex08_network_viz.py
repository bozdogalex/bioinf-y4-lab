#!/usr/bin/env python3
"""
Exercise 8 — Visualization of Co-Expression Networks + Hub Genes

Autor: olezhka0809

Pași:
- Încarcă matricea de expresie și modulele (din Lab 6)
- Reconstruiește matricea de adiacență din corelații
- Construiește graful de co-expresie
- Colorează nodurile după modul
- Calculează genele hub (top degree)
- Desenează și salvează rețeaua (.png)
- Salvează hub genes în CSV
"""

from __future__ import annotations
from pathlib import Path
from typing import Dict, Iterable, Optional

import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt


# --------------------------
# Config — complete with your values
# --------------------------
HANDLE = "olezhka0809"

# Input files – le folosim pe cele din labul 6
EXPR_CSV = Path(f"/workspaces/bioinf-y4-lab/labs/06_wgcna/submissions/olezhka0809/expression_matrix.csv")
MODULES_CSV = Path(f"/workspaces/bioinf-y4-lab/labs/06_wgcna/submissions/olezhka0809/modules_olezhka0809.csv")

# Optional: dacă ai salva adiacența în Lab 6 ai putea pune aici fișierul
PRECOMPUTED_ADJ_CSV: Optional[Path] = None

# Parameters for adjacency reconstruction (if PRECOMPUTED_ADJ_CSV is None)
CORR_METHOD = "spearman"   # "pearson" sau "spearman"
USE_ABS_CORR = True        # folosim |cor| ?
ADJ_THRESHOLD = 0.6        # prag pe corelație
WEIGHTED = True            # weighted adjacency (edge weight = corr) sau binar

# Visualization parameters
SEED = 42
TOPK_HUBS = 10
NODE_BASE_SIZE = 60
EDGE_ALPHA = 0.15

# Outputs (atenție la numele folderului: 07_network_viz)
OUT_DIR = Path(f"/workspaces/bioinf-y4-lab/labs/07_network_viz/submissions/{HANDLE}")
OUT_PNG = OUT_DIR / f"network_{HANDLE}.png"
OUT_HUBS = OUT_DIR / f"hubs_{HANDLE}.csv"


# --------------------------
# Utils
# --------------------------
def ensure_exists(path: Path) -> None:
    """Verifică dacă fișierul există; altfel, oprește scriptul cu un mesaj clar."""
    if not path.exists():
        raise FileNotFoundError(f"[ERROR] Nu am găsit fișierul: {path.resolve()}")


def read_expression_matrix(path: Path) -> pd.DataFrame:
    """
    Citește matricea de expresie.
    Asumăm:
      - rânduri = gene
      - coloane = probe
      - prima coloană e index (gene)
    """
    df = pd.read_csv(path, index_col=0)
    print(f"[INFO] Matrice expresie încărcată: {df.shape[0]} gene × {df.shape[1]} probe")
    return df


def read_modules_csv(path: Path) -> Dict[str, int]:
    """
    Citește CSV cu coloane: Gene, Module
    Returnează: dict gene -> module_id
    """
    df = pd.read_csv(path)
    if not {"Gene", "Module"}.issubset(df.columns):
        raise ValueError(f"[ERROR] {path} trebuie să conțină coloanele 'Gene' și 'Module'")
    mapping = dict(zip(df["Gene"], df["Module"]))
    print(f"[INFO] Mapping module încărcat: {len(mapping)} gene cu module.")
    return mapping


def correlation_to_adjacency(expr: pd.DataFrame,
                             method: str,
                             use_abs: bool,
                             threshold: float,
                             weighted: bool) -> pd.DataFrame:
    """
    - calculează matricea de corelație între gene (rânduri)
    - opțional aplică abs()
    - aplică threshold pentru a construi adiacența
    - setează diagonala la 0
    """
    print(f"[INFO] Calculez corelația ({method}) pentru {expr.shape[0]} gene...")
    # corelație între gene => corr pe transpose (samples ca rânduri)
    corr = expr.T.corr(method=method)
    if use_abs:
        corr = corr.abs()

    # diagonala = 0 (nu vrem self-loops)
    np.fill_diagonal(corr.values, 0.0)

    # threshold
    if weighted:
        A = corr.where(corr >= threshold, 0.0)
    else:
        A = (corr >= threshold).astype(float)

    # asigurăm simetria
    A = (A + A.T) / 2.0

    # mici mesaje
    n = A.shape[0]
    num_edges = np.count_nonzero(np.triu(A.values, k=1))
    print(f"[INFO] Matrice adiacență (threshold={threshold}, weighted={weighted}): "
          f"{n} noduri, ~{num_edges} muchii (triunghiul superior).")
    return A


def graph_from_adjacency(A: pd.DataFrame) -> nx.Graph:
    """
    Transformă matricea de adiacență în graf NetworkX și elimină nodurile izolate.
    """
    # Pentru weighted adjacency, NetworkX ia automat valorile ca 'weight'
    G = nx.from_pandas_adjacency(A)

    isolates = list(nx.isolates(G))
    if isolates:
        print(f"[INFO] Elimin {len(isolates)} noduri izolate (fără muchii).")
        G.remove_nodes_from(isolates)

    print(f"[INFO] Graf final: {G.number_of_nodes()} noduri, {G.number_of_edges()} muchii.")
    return G


def color_map_from_modules(nodes: Iterable[str],
                           gene2module: Dict[str, int]) -> Dict[str, str]:
    """
    Atribuie o culoare fiecărui nod în funcție de modul.
    Folosim paleta 'tab10' (sau mai multe cicluri dacă sunt >10 module).
    """
    nodes = list(nodes)
    modules = {gene2module.get(n, -1) for n in nodes}
    if -1 in modules:
        print("[WARN] Unele gene nu au modul asociat (vor fi gri).")

    modules_sorted = sorted(m for m in modules if m != -1)
    cmap = plt.cm.get_cmap("tab10", max(len(modules_sorted), 1))

    module2color: Dict[int, str] = {}
    for i, m in enumerate(modules_sorted):
        module2color[m] = cmap(i)

    default_color = "#cccccc"  # pentru gene fără modul

    node2color: Dict[str, str] = {}
    for n in nodes:
        m = gene2module.get(n, -1)
        node2color[n] = module2color.get(m, default_color)

    return node2color


def compute_hubs(G: nx.Graph, topk: int) -> pd.DataFrame:
    """
    - calculează degree pentru fiecare nod
    - calculează betweenness centrality (opțional dar util pentru raport)
    - întoarce top-k gene ordonate după degree desc.
    """
    print("[INFO] Calculez nodurile hub (degree + betweenness)...")
    # Dacă ai weighted edges și vrei degree ponderat, pune weight="weight"
    degree_dict = dict(G.degree(weight="weight" if WEIGHTED else None))
    betw_dict = nx.betweenness_centrality(G, weight="weight" if WEIGHTED else None)

    df = pd.DataFrame({
        "Gene": list(degree_dict.keys()),
        "Degree": list(degree_dict.values()),
        "Betweenness": [betw_dict[g] for g in degree_dict.keys()],
    })

    df = df.sort_values(["Degree", "Betweenness"], ascending=False).head(topk).reset_index(drop=True)
    print("[INFO] Top hub genes:")
    for _, row in df.iterrows():
        print(f"  {row['Gene']}: degree={row['Degree']}, betw={row['Betweenness']:.3f}")
    return df


# --------------------------
# Main
# --------------------------
if __name__ == "__main__":
    print("\n" + "=" * 46)
    print(" GENE CO-EXPRESSION NETWORKS (GCE) - LAB 07 ")
    print("=" * 46 + "\n")

    # 1) Verifică input files
    ensure_exists(EXPR_CSV)
    ensure_exists(MODULES_CSV)
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # 2) Load expression matrix și module mapping
    expr_all = read_expression_matrix(EXPR_CSV)
    gene2module = read_modules_csv(MODULES_CSV)

    # limităm la gene care apar și în module, și în expresie
    common_genes = sorted(set(expr_all.index) & set(gene2module.keys()))
    if not common_genes:
        raise ValueError("[ERROR] Nu există gene comune între expression_matrix și modules CSV!")
    expr = expr_all.loc[common_genes]
    print(f"[INFO] Lucrez cu {expr.shape[0]} gene (cu module) × {expr.shape[1]} probe")

    # 3) Load or reconstruct adjacency
    if PRECOMPUTED_ADJ_CSV is not None:
        ensure_exists(PRECOMPUTED_ADJ_CSV)
        A_full = pd.read_csv(PRECOMPUTED_ADJ_CSV, index_col=0)
        # sub-matrice doar pe genele comune
        A = A_full.loc[common_genes, common_genes]
        print(f"[INFO] Am încărcat adiacența precomputată: {A.shape[0]}×{A.shape[1]}")
    else:
        A = correlation_to_adjacency(
            expr=expr,
            method=CORR_METHOD,
            use_abs=USE_ABS_CORR,
            threshold=ADJ_THRESHOLD,
            weighted=WEIGHTED,
        )

    # 4) Build graph
    G = graph_from_adjacency(A)

    # 5) Compute colors by module
    node_color_map = color_map_from_modules(G.nodes(), gene2module)
    node_colors = [node_color_map[n] for n in G.nodes()]

    # 6) Compute hub genes
    hubs_df = compute_hubs(G, TOPK_HUBS)
    hub_genes = set(hubs_df["Gene"])

    # Adaugă module ca atribut pe noduri (dacă ai gene2module)
    nx.set_node_attributes(G, gene2module, "module")

    # Adaugă degree și betweenness ca atribute (le ai deja în compute_hubs)
    degree_dict = dict(G.degree(weight="weight"))
    betw_dict = nx.betweenness_centrality(G, weight="weight", normalized=True)

    nx.set_node_attributes(G, degree_dict, "degree")
    nx.set_node_attributes(G, betw_dict, "betweenness")

    # Marchează hub-urile
    hub_genes = set(hubs_df["Gene"])
    hub_flag = {g: (g in hub_genes) for g in G.nodes()}
    nx.set_node_attributes(G, hub_flag, "is_hub")


    # node sizes proporțional cu degree
    degree_dict = dict(G.degree(weight="weight" if WEIGHTED else None))
    node_sizes = [NODE_BASE_SIZE + 40 * degree_dict[n] for n in G.nodes()]
    

    # 1) GraphML – recomandat pentru Cytoscape și Gephi
    graphml_path = OUT_DIR / f"network_{HANDLE}.graphml"
    nx.write_graphml(G, graphml_path)
    print(f"[OK] GraphML salvat în: {graphml_path}")

    # 2) (Opțional) Edge list CSV – dacă preferi CSV
    edges_df = nx.to_pandas_edgelist(G)
    edges_csv = OUT_DIR / f"network_edges_{HANDLE}.csv"
    edges_df.to_csv(edges_csv, index=False)
    print(f"[OK] Edge list CSV salvat în: {edges_csv}")

    # 3) (Opțional) Node attributes CSV
    nodes_df = pd.DataFrame.from_dict(dict(G.nodes(data=True)), orient="index")
    nodes_df.index.name = "Gene"
    nodes_csv = OUT_DIR / f"network_nodes_{HANDLE}.csv"
    nodes_df.to_csv(nodes_csv)
    print(f"[OK] Node attributes CSV salvat în: {nodes_csv}")

    # 7) Layout + draw graph
    print("[INFO] Desenez rețeaua...")
    plt.figure(figsize=(10, 8))
    pos = nx.spring_layout(G, seed=SEED)

    # edges
    nx.draw_networkx_edges(
        G, pos,
        alpha=EDGE_ALPHA,
        width=1.0,
        edge_color="gray"
    )

    # nodes
    nx.draw_networkx_nodes(
        G, pos,
        node_color=node_colors,
        node_size=node_sizes,
        edgecolors="black",
        linewidths=0.5
    )

    # labels doar pentru hub genes
    labels = {n: n for n in G.nodes() if n in hub_genes}
    nx.draw_networkx_labels(
        G, pos,
        labels=labels,
        font_size=8
    )

    plt.title(f"Gene Co-Expression Network (TP53-related) — {HANDLE}")
    plt.axis("off")
    plt.tight_layout()
    plt.savefig(OUT_PNG, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"[OK] Rețea salvată în: {OUT_PNG}")

    # 8) Save hub genes
    hubs_df.to_csv(OUT_HUBS, index=False)
    print(f"[OK] Hub genes salvate în: {OUT_HUBS}\n")

    print("=" * 46)
    print(" GATA — poți folosi PNG + CSV în raport ")
    print("=" * 46)
