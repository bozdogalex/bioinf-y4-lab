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
import os
import requests
import gzip
from io import StringIO

def download_and_prepare():
    handle = "blejdeaxoxo"
    out_path = f"data/work/{handle}/lab06"
    os.makedirs(out_path, exist_ok=True)
    out_file = f"{out_path}/expression_matrix.csv"
    
    url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE115nnn/GSE115469/suppl/GSE115469_Data.csv.gz"
    gz_file = os.path.join(out_path, "GSE115469_Data.csv.gz")
    
    print("Downloading GSE115469_Data.csv.gz …")
    r = requests.get(url, stream=True)
    r.raise_for_status()
    with open(gz_file, "wb") as f:
        for chunk in r.iter_content(chunk_size=1024*1024):
            f.write(chunk)
    print("Downloaded supplementary data.")
    
    print("Reading CSV …")
    with gzip.open(gz_file, "rt") as f:
        df = pd.read_csv(f, index_col=0)
    
    print("Dataframe loaded: shape =", df.shape)
    
    df.to_csv(out_file)
    print(f"Expression matrix saved to: {out_file}")



# --------------------------
# Config — completați după nevoie
# --------------------------
INPUT_CSV = "/workspaces/bioinf-y4-lab/data/work/blejdeaxoxo/lab06/expression_matrix.csv"
OUTPUT_DIR = "/workspaces/bioinf-y4-lab/labs/06_wgcna/submissions/blejdeaxoxo"
OUTPUT_CSV = "/workspaces/bioinf-y4-lab/labs/06_wgcna/submissions/blejdeaxoxo/modules_blejdeaxoxo.csv"

CORR_METHOD = "spearman"   # TODO: "pearson" sau "spearman"
VARIANCE_THRESHOLD =  0.7 # prag pentru filtrare gene
ADJ_THRESHOLD =   0.6  # prag pentru |cor| (ex: 0.6)
USE_ABS_CORR =   True     # True => folosiți |cor| la prag
MAKE_UNDIRECTED =  True    # rețelele de co-expresie sunt de obicei neorientate
MAX_GENES = 5000

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


def detect_modules_louvain_or_greedy(G: nx.Graph) -> Dict[str, int]:
    """
    TODO: detectați comunități (module) și întoarceți un dict gene -> modul_id.
    Variante:
      - încercați louvain_communities(G, seed=42) dacă e disponibil
      - altfel greedy_modularity_communities(G)
    """
    if G is None:
        raise ValueError("G must be a networkx Graph, got None")
    if G.number_of_nodes() == 0:
        return {}
    try:
        from networkx.algorithms.community import louvain_communities
        communities = list(louvain_communities(G, weight="weight", seed=42))
    except Exception:
        from networkx.algorithms.community import greedy_modularity_communities
        communities = list(greedy_modularity_communities(G, weight="weight"))

    communities_sorted = sorted(
        communities,
        key=lambda c: (-len(c), min(map(str, c)))
    )

    mapping: Dict[str, int] = {}
    for module_id, community in enumerate(communities_sorted):
        for node in community:
            mapping[str(node)] = int(module_id)

    print(f"Detected {len(communities_sorted)} modules; nodes assigned={len(mapping)}")

    return mapping


def save_modules_csv(mapping: Dict[str, int], out_csv: Path) -> None:
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    df_modules = (
        pd.DataFrame({"Gene": list(mapping.keys()), "Module": list(mapping.values())})
        .sort_values(["Module", "Gene"])
    )
    df_modules.to_csv(out_csv, index=False)


if __name__ == "__main__":
    df = read_expression_matrix(Path(INPUT_CSV))
    df = log_and_filter(df, VARIANCE_THRESHOLD)
    corr = correlation_matrix(df, method=CORR_METHOD, use_abs=USE_ABS_CORR)
    adjacency = adjacency_from_correlation(corr, ADJ_THRESHOLD, weighted=True)
    G = graph_from_adjacency(adjacency, undirected=MAKE_UNDIRECTED)
    print(f"Grafic creat cu {G.number_of_nodes()} noduri și {G.number_of_edges()} muchii.")

    gene_to_module = detect_modules_louvain_or_greedy(G)
    print(f"S-au detectat {len(set(gene_to_module.values()))} module.")

    save_modules_csv(gene_to_module, Path(OUTPUT_CSV))
    print(f"Am salvat mapping-ul gene→modul în: {OUTPUT_CSV}")
