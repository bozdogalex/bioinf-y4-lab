"""
Exercise 8 — Visualization of Co-Expression Networks + Hub Genes

TODO:
- Load the expression matrix and module mapping from Lab 6
- Rebuild (or load) the adjacency matrix
- Construct the graph from adjacency
- Color nodes by module
- Compute hub genes (top degree)
- Visualize and export the network figure (.png)
- Export hub genes to CSV (optional)
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
HANDLE = "blejdeaxoxo"

# Input files
EXPR_CSV = Path(f"/workspaces/bioinf-y4-lab/data/work/blejdeaxoxo/lab06/expression_matrix.csv")
MODULES_CSV = Path(f"/workspaces/bioinf-y4-lab/labs/06_wgcna/submissions/blejdeaxoxo/modules_blejdeaxoxo.csv")

# Optional: if you saved adjacency in Lab 6, load it here
PRECOMPUTED_ADJ_CSV: Optional[Path] = None

# Parameters for adjacency reconstruction (if PRECOMPUTED_ADJ_CSV is None)
CORR_METHOD = "spearman"   # TODO: choose "pearson" or "spearman"
USE_ABS_CORR = True        # TODO: use absolute correlations?
ADJ_THRESHOLD = 0.6        # TODO: correlation threshold
WEIGHTED = True           # TODO: use weighted adjacency or binary?
VARIANCE_THRESHOLD =  0.7

MAX_GENES = 5000

# Visualization parameters
SEED = 42
TOPK_HUBS = 10
NODE_BASE_SIZE = 60
EDGE_ALPHA = 0.15

# Outputs
OUT_DIR = Path(f"/workspaces/bioinf-y4-lab/labs/07_network_viz/submissions/blejdeaxoxo")
OUT_PNG = OUT_DIR / f"network_{HANDLE}.png"
OUT_HUBS = OUT_DIR / f"hubs_{HANDLE}.csv"


# --------------------------
# Utils
# --------------------------
def ensure_exists(path: Path) -> None:
    if path is None:
        raise ValueError("path must be provided")
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"Required file not found: {p}")
    if not p.is_file():
        raise FileNotFoundError(f"Required path exists but is not a file: {p}")


def read_expression_matrix(path: Path, max_genes: Optional[int] = MAX_GENES) -> pd.DataFrame:
    if path is None:
        raise ValueError("path must be provided")

    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"Expression matrix not found: {p}")

    if max_genes is None or max_genes <= 0:
        df = pd.read_csv(p, index_col=0)
        truncated = False
    else:
        nrows_to_read = max_genes + 1
        df = pd.read_csv(p, index_col=0, nrows=nrows_to_read)
        if df.shape[0] > max_genes:
            truncated = True
            df = df.iloc[:max_genes].copy()
        else:
            truncated = False

    if df.shape[0] == 0 or df.shape[1] == 0:
        raise ValueError(f"Empty expression matrix: {p}")

    if truncated:
        print(
            f"Warning: expression matrix has more than {max_genes} genes; "
            "loading first {0} rows only to avoid OOM".format(max_genes)
        )

    df.index = df.index.map(str)

    return df


def read_modules_csv(path: Path) -> Dict[str, int]:
    if path is None:
        raise ValueError("path must be provided")
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"Modules CSV not found: {p}")

    df = pd.read_csv(p)
    if "Gene" not in df.columns or "Module" not in df.columns:
        raise ValueError(f"Expected columns 'Gene' and 'Module' in {p}; got {list(df.columns)}")

    df = df[["Gene", "Module"]].dropna(how="any")

    mapping: Dict[str, int] = {}
    duplicates = []
    for _, row in df.iterrows():
        gene = str(row["Gene"])
        module_val = row["Module"]
        try:
            module_id = int(module_val)
        except Exception:
            try:
                module_id = int(float(module_val))
            except Exception:
                raise ValueError(f"Cannot convert Module value to int for gene {gene}: {module_val!r}")

        if gene in mapping:
            duplicates.append(gene)
        mapping[gene] = module_id

    if duplicates:
        print(f"Warning: duplicate gene entries in {p}; last occurrence used for {len(duplicates)} genes")

    return mapping


def log_and_filter(df: pd.DataFrame,
                   variance_threshold: float) -> pd.DataFrame:
    if df is None:
        raise ValueError("df must be a pandas DataFrame, got None")
    
    df_numeric = df.apply(pd.to_numeric, errors="coerce")
    if df_numeric.isna().any(axis=None):
        df_numeric = df_numeric.dropna(axis=0, how="any")

    df_log = np.log2(df_numeric + 1.0)

    variances = df_log.var(axis=1, ddof=0)

    keep_mask = variances >= variance_threshold
    filtered = df_log.loc[keep_mask].copy()

    return filtered

def correlation_matrix(df: pd.DataFrame,
                       method: str = "spearman",
                       use_abs: bool = True) -> pd.DataFrame:
    if df is None:
        raise ValueError("df must be a pandas DataFrame, got None")
    if df.shape[0] == 0:
        return pd.DataFrame(index=df.index, columns=df.index, dtype=float)

    method = method.lower()
    if method not in {"pearson", "spearman"}:
        raise ValueError("method must be one of: 'pearson', 'spearman'")

    corr = df.T.corr(method=method)
    corr = corr.fillna(0.0)
    np.fill_diagonal(corr.values, 1.0)

    if use_abs:
        corr = corr.abs()

    corr = corr.astype(float)
    return corr


def adjacency_from_correlation(corr: pd.DataFrame,
                               threshold: float,
                               weighted: bool = False) -> pd.DataFrame:
    """
    Construiți matricea de adiacență din corelații.
    - binară: A_ij = 1 dacă corr_ij >= threshold, altfel 0
    - ponderată: A_ij = corr_ij dacă corr_ij >= threshold, altfel 0
    """
    if corr is None:
        raise ValueError("corr must be a pandas DataFrame, got None")
    if corr.shape[0] == 0:
        return corr.copy()

    mask = corr >= threshold
    if weighted:
        adjacency = corr.where(mask, other=0.0)
    else:
        adjacency = mask.astype(float)

    np.fill_diagonal(adjacency.values, 0.0)

    return adjacency.astype(float)


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


def color_map_from_modules(nodes: Iterable[str], gene2module: Dict[str, int]) -> Dict[str, str]:
    from matplotlib import colors as mcolors

    if gene2module is None:
        gene2module = {}

    node_list = [str(n) for n in nodes]
    modules = sorted({gene2module[n] for n in node_list if n in gene2module})

    if not modules:
        return {n: "#808080" for n in node_list}  

    cmap_name = "tab10" if len(modules) <= 10 else "tab20"
    cmap = plt.get_cmap(cmap_name)
    cmap_size = cmap.N

    module2color: Dict[int, str] = {}
    for i, mod in enumerate(modules):
        rgba = cmap(i % cmap_size)
        module2color[mod] = mcolors.to_hex(rgba)

    node2color: Dict[str, str] = {}
    for n in node_list:
        if n in gene2module:
            node2color[n] = module2color[gene2module[n]]
        else:
            node2color[n] = "#808080"

    return node2color


def compute_hubs(G: nx.Graph, topk: int) -> pd.DataFrame:
    if G is None:
        raise ValueError("G must be provided")
    if G.number_of_nodes() == 0:
        return pd.DataFrame(columns=["Gene", "Degree", "WeightedDegree", "Betweenness"])

    deg = dict(G.degree())
    wdeg = dict(G.degree(weight="weight"))

    bet = nx.betweenness_centrality(G, weight="weight") if G.number_of_edges() > 0 else {n: 0.0 for n in G.nodes()}

    rows = []
    for n in G.nodes():
        rows.append({
            "Gene": str(n),
            "Degree": int(deg.get(n, 0)),
            "WeightedDegree": float(wdeg.get(n, 0.0)),
            "Betweenness": float(bet.get(n, 0.0)),
        })

    df = pd.DataFrame(rows)
    df = df.sort_values(by=["WeightedDegree", "Degree", "Betweenness"], ascending=[False, False, False]).reset_index(drop=True)

    if topk is None or topk <= 0:
        topk = len(df)

    return df.head(min(topk, len(df))).copy()


# --------------------------
# Main
# --------------------------
if __name__ == "__main__":
    # 1. Verify input files exist
    ensure_exists(EXPR_CSV)
    ensure_exists(MODULES_CSV)
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # 2. Load expression matrix and module mapping
    try:
        expr = read_expression_matrix(EXPR_CSV, max_genes=5000)
    except Exception as e:
        raise RuntimeError(f"Failed to load expression matrix: {e}")

    try:
        gene2module = read_modules_csv(MODULES_CSV)
    except Exception as e:
        raise RuntimeError(f"Failed to load modules CSV: {e}")

    print(f"Expression matrix: genes={expr.shape[0]}, samples={expr.shape[1]}")
    print(f"Module assignments: genes={len(gene2module)}")

    # 3. Load or reconstruct adjacency
    if PRECOMPUTED_ADJ_CSV is not None and PRECOMPUTED_ADJ_CSV.exists():
        A = pd.read_csv(PRECOMPUTED_ADJ_CSV, index_col=0)
        print(f"Loaded precomputed adjacency: {A.shape}")
    else:
        expr1 = log_and_filter(expr, VARIANCE_THRESHOLD)
        corr = correlation_matrix(expr1, method=CORR_METHOD, use_abs=USE_ABS_CORR)
        A = adjacency_from_correlation(corr, ADJ_THRESHOLD, weighted=True)
        print(f"Computed adjacency: {A.shape}")

    # Restrict adjacency to genes present in both expr and module mapping (if any)
    genes_in_data = set(map(str, expr.index))
    genes_in_adj = set(map(str, A.index))
    genes_in_modules = set(map(str, gene2module.keys()))
    keep_genes = sorted(list(genes_in_adj & genes_in_data))
    if genes_in_modules:
        keep_genes = sorted(list(set(keep_genes) & genes_in_modules))

    if not keep_genes:
        raise RuntimeError("No overlapping genes between expression/adjacency/modules")

    A = A.loc[keep_genes, keep_genes]
    print(f"Filtered adjacency to {len(keep_genes)} genes")

    # 4. Build graph
    G = graph_from_adjacency(A, undirected=True)
    print(f"Graph: nodes={G.number_of_nodes()}, edges={G.number_of_edges()}")

    if G.number_of_nodes() == 0:
        raise RuntimeError("Empty graph — nothing to visualize")

    # 5. Compute colors by module
    node_color_map = color_map_from_modules(G.nodes(), gene2module)
    node_colors = [node_color_map.get(str(n), "#808080") for n in G.nodes()]

    # 6. Compute hub genes
    hubs_df = compute_hubs(G, TOPK_HUBS)
    hubs_df.to_csv(OUT_HUBS, index=False)
    print(f"Saved top hubs to: {OUT_HUBS}")

    deg_dict = dict(G.degree())
    node_sizes = [NODE_BASE_SIZE * (1 + deg_dict.get(n, 0)) for n in G.nodes()]

    # 7. Compute layout and draw graph
    plt.figure(figsize=(12, 10))
    pos = nx.spring_layout(G, seed=SEED)

    if WEIGHTED:
        weights = [d.get("weight", 0.0) for _, _, d in G.edges(data=True)]

        widths = [max(0.5, float(w) * 3.0) for w in weights]
        nx.draw_networkx_edges(G, pos, alpha=EDGE_ALPHA, width=widths, edge_color="#888888")
    else:
        nx.draw_networkx_edges(G, pos, alpha=EDGE_ALPHA, width=0.5, edge_color="#888888")

    nx.draw_networkx_nodes(G, pos,
                           node_color=node_colors,
                           node_size=node_sizes,
                           linewidths=0.5,
                           edgecolors="#222222")

    hub_genes = list(hubs_df["Gene"].astype(str).values)
    hub_labels = {n: str(n) for n in G.nodes() if str(n) in hub_genes}
    if hub_labels:
        nx.draw_networkx_labels(G, pos, labels=hub_labels, font_size=8)

    plt.axis("off")
    plt.title(f"Co-expression network ({HANDLE}) — nodes={G.number_of_nodes()}, edges={G.number_of_edges()}", fontsize=12)

    # 8. Save network figure
    plt.tight_layout()
    plt.savefig(OUT_PNG, dpi=200)
    plt.close()
    print(f"Saved network figure to: {OUT_PNG}")

    # 9. Completion
    print("Done.")