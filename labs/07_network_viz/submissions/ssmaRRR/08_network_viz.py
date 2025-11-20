from __future__ import annotations
from pathlib import Path
from typing import Dict, Iterable, Optional

import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt


# --------------------------
# Config
# --------------------------
HANDLE = "ssmaRRR"

# Fișier expresie (același ca în Lab 6)
EXPR_CSV = Path(f"data/work/{HANDLE}/lab06/expression_matrix.csv")

# Mapping gene -> modul din Lab 6 (ai folosit 06_wgcna)
MODULES_CSV = Path(f"labs/06_wgcna/submissions/{HANDLE}/modules_{HANDLE}.csv")

# Nu folosim adiacență pre-calculată (reconstruim din expresie)
PRECOMPUTED_ADJ_CSV: Optional[Path] = None

# Parametri preprocesare (identici cu Lab 6)
CORR_METHOD = "spearman"      # "pearson" sau "spearman"
VARIANCE_THRESHOLD = 0.5      # prag pentru filtrare gene (var > 0.5)
USE_ABS_CORR = True           # True => folosim |cor|
ADJ_THRESHOLD = 0.6           # prag pentru |cor|
WEIGHTED = False              # False => 0/1; True => păstrează valorile corr peste prag

# Parametri de vizualizare
SEED = 42                     # pentru layout determinist
TOPK_HUBS = 10                # câte gene hub etichetăm
NODE_BASE_SIZE = 60           # mărimea de bază a nodurilor
EDGE_ALPHA = 0.15             # transparența muchiilor

# Ieșiri
OUT_DIR = Path(f"labs/07_network_viz/submissions/{HANDLE}")
OUT_PNG = OUT_DIR / f"network_{HANDLE}.png"
#OUT_HUBS = OUT_DIR / f"hubs_{HANDLE}.csv"


# --------------------------
# Utils
# --------------------------
def ensure_exists(path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Nu am găsit: {path}")


def read_expression_matrix(path: Path) -> pd.DataFrame:
    ensure_exists(path)
    df = pd.read_csv(path, index_col=0)
    if df.empty:
        raise ValueError("Matricea de expresie este goală.")
    return df


def log_and_filter(df: pd.DataFrame,
                   variance_threshold: float) -> pd.DataFrame:
    """
    Preprocesare identică cu Lab 6:
    - log2(x+1)
    - păstrăm doar genele cu varianță > variance_threshold
    """
    df_log = np.log2(df + 1)
    df_filt = df_log.loc[df_log.var(axis=1) > variance_threshold]
    if df_filt.empty:
        raise ValueError("După filtrare nu a mai rămas nicio genă. Verificați pragul de varianță.")
    return df_filt


def correlation_matrix(df: pd.DataFrame,
                       method: str = "spearman",
                       use_abs: bool = True) -> pd.DataFrame:
    """
    Corelație între gene (rânduri):
    - df: (gene x probe)
    - întoarce matrice (gene x gene)
    """
    corr = df.T.corr(method=method)
    if use_abs:
        corr = corr.abs()
    return corr


def adjacency_from_correlation(corr: pd.DataFrame,
                               threshold: float,
                               weighted: bool = False) -> pd.DataFrame:
    """
    Matrice de adiacență din corelații:
    - binară: 1 dacă corr_ij >= threshold, altfel 0
    - ponderată: corr_ij dacă corr_ij >= threshold, altfel 0
    """
    if weighted:
        A = corr.copy()
        A[A < threshold] = 0.0
    else:
        A = (corr >= threshold).astype(float)
    np.fill_diagonal(A.values, 0.0)
    return A


def graph_from_adjacency(A: pd.DataFrame) -> nx.Graph:
    """
    Creează graf neorientat din adiacență și elimină nodurile izolate.
    """
    G = nx.from_pandas_adjacency(A)
    isolates = list(nx.isolates(G))
    if isolates:
        G.remove_nodes_from(isolates)
    return G


def read_modules_csv(path: Path) -> Dict[str, int]:
    """
    Așteaptă un CSV cu coloanele: Gene, Module
    (exact cum ai salvat în Lab 6).
    """
    ensure_exists(path)
    df = pd.read_csv(path)
    if not {"Gene", "Module"}.issubset(df.columns):
        raise ValueError("modules.csv trebuie să conțină coloanele: Gene, Module")
    return dict(zip(df["Gene"].astype(str), df["Module"].astype(int)))


def color_map_from_modules(nodes: Iterable[str],
                           gene2module: Dict[str, int]) -> Dict[str, str]:
    """
    Nod -> culoare pe baza modulului.
    Gene fără modul -> gri.
    """
    cmap = plt.get_cmap("tab10")
    colors: Dict[str, str] = {}
    for n in nodes:
        m = gene2module.get(n, 0)
        colors[n] = cmap((m - 1) % 10) if m > 0 else "#CCCCCC"
    return colors


def compute_hubs(G: nx.Graph, topk: int) -> pd.DataFrame:
    """
    Gene hub = noduri cu grad mare.
    Returnează DataFrame cu Gene, Degree, Betweenness.
    """
    deg_dict = dict(G.degree())

    if G.number_of_nodes() <= 5000:
        btw_dict = nx.betweenness_centrality(G, normalized=True, seed=SEED)
    else:
        btw_dict = {n: np.nan for n in G.nodes()}

    hubs = (
        pd.DataFrame(
            {
                "Gene": list(deg_dict.keys()),
                "Degree": list(deg_dict.values()),
                "Betweenness": [btw_dict.get(n, np.nan) for n in deg_dict.keys()],
            }
        )
        .sort_values(["Degree", "Betweenness"], ascending=False)
        .head(topk)
        .reset_index(drop=True)
    )
    return hubs


# --------------------------
# Main
# --------------------------
if __name__ == "__main__":
    # 1) Verificări & output dir
    ensure_exists(EXPR_CSV)
    ensure_exists(MODULES_CSV)
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # 2) Expresie + preprocesare ca în Lab 6
    expr_raw = read_expression_matrix(EXPR_CSV)
    expr = log_and_filter(expr_raw, variance_threshold=VARIANCE_THRESHOLD)

    # 3) Mapping gene -> modul din Lab 6
    gene2module = read_modules_csv(MODULES_CSV)

    # 4) Adiacență: dacă ai vreodată PRECOMPUTED_ADJ_CSV, o poți folosi aici
    if PRECOMPUTED_ADJ_CSV is not None:
        ensure_exists(PRECOMPUTED_ADJ_CSV)
        A = pd.read_csv(PRECOMPUTED_ADJ_CSV, index_col=0)
    else:
        corr = correlation_matrix(expr, method=CORR_METHOD, use_abs=USE_ABS_CORR)
        A = adjacency_from_correlation(corr, threshold=ADJ_THRESHOLD, weighted=WEIGHTED)

    # Păstrăm doar genele care apar și în modules.csv
    common_genes = sorted(set(A.index) & set(gene2module.keys()))
    if not common_genes:
        raise ValueError("Nu există gene comune între adiacență și modules.csv.")
    A = A.loc[common_genes, common_genes]

    # 5) Graf
    G = graph_from_adjacency(A)
    print(f"Grafic: {G.number_of_nodes()} noduri, {G.number_of_edges()} muchii.")

    if G.number_of_nodes() == 0:
        raise ValueError("Graful rezultat nu are noduri (toate genele sunt izolate). Verifică pragul de corelație.")

    # 6) Culori după modul
    node_colors_map = color_map_from_modules(G.nodes(), gene2module)
    node_colors = [node_colors_map[n] for n in G.nodes()]

    # 7) Gene hub
    hubs_df = compute_hubs(G, TOPK_HUBS)
    hubs_set = set(hubs_df["Gene"])
    node_sizes = [
        NODE_BASE_SIZE * (1.5 if n in hubs_set else 1.0)
        for n in G.nodes()
    ]

    # 8) Layout + desen
    pos = nx.spring_layout(G, seed=SEED)

    plt.figure(figsize=(12, 10))
    nx.draw_networkx_edges(G, pos, alpha=EDGE_ALPHA, width=0.5)
    nx.draw_networkx_nodes(
        G,
        pos,
        node_color=node_colors,
        node_size=node_sizes,
        linewidths=0.0,
    )

    # etichete doar pentru hub-uri
    hub_labels = {g: g for g in hubs_set if g in G.nodes()}
    nx.draw_networkx_labels(G, pos, labels=hub_labels, font_size=8)

    plt.axis("off")
    plt.tight_layout()
    plt.savefig(OUT_PNG, dpi=300)
    plt.close()
    print(f"Am salvat figura în: {OUT_PNG}")

    # 9) Export hub-uri
    #hubs_df.to_csv(OUT_HUBS, index=False)
    #print(f"Am salvat hub genes în: {OUT_HUBS}")
