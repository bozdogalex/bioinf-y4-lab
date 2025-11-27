from __future__ import annotations
from pathlib import Path
from typing import Dict
import numpy as np
import pandas as pd
import networkx as nx
from sklearn.feature_selection import VarianceThreshold
try:
    from community import community_louvain
except ImportError:
    community_louvain = None

# --------------------------
# Config
# --------------------------
INPUT_CSV = "data/work/MadalinaDenisa/lab06/expression_matrix.csv"
OUTPUT_DIR = Path("labs/06_networks/submissions/MadalinaDenisa")
OUTPUT_CSV = OUTPUT_DIR / "modules_MadalinaDenisa.csv"

CORR_METHOD = "spearman"   
VARIANCE_THRESHOLD = 0.01  # modificare pentru a nu elimina toate genele
ADJ_THRESHOLD = 0.6     
USE_ABS_CORR = True       
MAKE_UNDIRECTED = True    

def read_expression_matrix(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, index_col=0)
    return df

def log_and_filter(df: pd.DataFrame, variance_threshold: float) -> pd.DataFrame:
    df_log = np.log2(df + 1)
    selector = VarianceThreshold(threshold=variance_threshold)
    filtered = selector.fit_transform(df_log.T).T
    genes_kept = df_log.index[selector.get_support()]
    return pd.DataFrame(filtered, index=genes_kept, columns=df.columns)

def correlation_matrix(df: pd.DataFrame, method: str = "spearman", use_abs: bool = True) -> pd.DataFrame:
    corr = df.T.corr(method=method)
    if use_abs:
        corr = corr.abs()
    return corr

def adjacency_from_correlation(corr: pd.DataFrame, threshold: float, weighted: bool = False) -> pd.DataFrame:
    if weighted:
        A = corr.where(corr >= threshold, 0)
    else:
        A = (corr >= threshold).astype(int)
    return A

def graph_from_adjacency(A: pd.DataFrame, undirected: bool = True) -> nx.Graph:
    if undirected:
        G = nx.from_pandas_adjacency(A)
    else:
        G = nx.from_pandas_adjacency(A, create_using=nx.DiGraph)
    isolates = list(nx.isolates(G))
    if isolates:
        G.remove_nodes_from(isolates)
    return G

def detect_modules_louvain_or_greedy(G: nx.Graph) -> Dict[str, int]:
    if community_louvain:
        partition = community_louvain.best_partition(G)
        return partition
    else:
        from networkx.algorithms.community import greedy_modularity_communities
        comms = greedy_modularity_communities(G)
        mapping = {}
        for i, comm in enumerate(comms):
            for gene in comm:
                mapping[gene] = i
        return mapping

def save_modules_csv(mapping: Dict[str, int], out_csv: Path) -> None:
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    df_modules = pd.DataFrame({"Gene": list(mapping.keys()), "Module": list(mapping.values())})
    df_modules.sort_values(["Module", "Gene"], inplace=True)
    df_modules.to_csv(out_csv, index=False)

if __name__ == "__main__":
    df = read_expression_matrix(INPUT_CSV)
    df_filtered = log_and_filter(df, VARIANCE_THRESHOLD)
    corr = correlation_matrix(df_filtered, CORR_METHOD, USE_ABS_CORR)
    A = adjacency_from_correlation(corr, ADJ_THRESHOLD)
    G = graph_from_adjacency(A, MAKE_UNDIRECTED)
    print(f"Graph created with {G.number_of_nodes()} nodes and {G.number_of_edges()} edges.")
    gene_to_module = detect_modules_louvain_or_greedy(G)
    print(f"Detected {len(set(gene_to_module.values()))} modules.")
    save_modules_csv(gene_to_module, OUTPUT_CSV)
    print(f"Saved modules → {OUTPUT_CSV}")
