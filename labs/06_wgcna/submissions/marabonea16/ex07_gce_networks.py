"""
Exercițiu Gene Co-Expression Networks (GCEs) — Construirea rețelei și detectarea modulelor

Obiectiv:
- Să construiți o rețea de co-expresie dintr-o matrice de expresie RNA-Seq
- Să detectați module (comunități) de gene folosind un algoritm de tip Louvain (sau alternativ)

Instrucțiuni (în laborator):
1) Pregătire date
   - Descărcați și pregătiți matricea de expresie (ex: GSE115469) într-un CSV cu:
     * rânduri = gene (index), coloane = probe (sample IDs)
   - Salvați fișierul la: data/work/marabonea16/lab06/expression_matrix.csv

2) Preprocesare
   - log2(x + 1)
   - filtrare gene cu varianță scăzută

3) Corelație → Adiacență
   - completați funcția `correlation_matrix`
   - funcția `adjacency_from_correlation` este deja implementată

4) Graf + Module
   - construiți graful cu NetworkX
   - detectați modulele (Louvain sau alternativă)
   - exportați mapping-ul gene → modul în submissions/marabonea16/modules_marabonea16.csv

Notă:
- Documentați în marabonea16_notes.md: metrica de corelație, pragul, observații scurte.
"""

from __future__ import annotations
from pathlib import Path
from typing import Dict, Iterable

import numpy as np
import pandas as pd
import networkx as nx
import urllib.request
import gzip
import io
import re


INPUT_CSV = Path("data/work/marabonea16/lab06/expression_matrix.csv")
OUTPUT_DIR = Path("labs/06_wgcna/submissions/marabonea16")
OUTPUT_CSV = OUTPUT_DIR / "modules_marabonea16.csv"

CORR_METHOD = "spearman"   # "pearson" sau "spearman"
VARIANCE_THRESHOLD = 0.6   # prag pentru filtrare gene
ADJ_THRESHOLD = 0.6        # prag pentru |cor| (ex: 0.6)
USE_ABS_CORR = True        # True => folosiți |cor| la prag
MAKE_UNDIRECTED = True     # rețelele de co-expresie sunt de obicei neorientate

TOP_N_DEFAULT = 200
SAMPLE_COLS_DEFAULT = 2000

def download_matrix(gse_id: str, out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    series_num = re.match(r"GSE(\d+)", gse_id).group(1)
    series_prefix = "GSE" + series_num[:-3] + "nnn"
    ftp_url = (
        f"ftp://ftp.ncbi.nlm.nih.gov/geo/series/{series_prefix}/{gse_id}/matrix/"
        f"{gse_id}_series_matrix.txt.gz"
    )

    with urllib.request.urlopen(ftp_url) as response:
        with gzip.GzipFile(fileobj=io.BytesIO(response.read())) as gz_file:
            df = pd.read_csv(gz_file, sep="\t", comment="!", index_col=0)

    if df.empty:
        suppl_url = (
            f"ftp://ftp.ncbi.nlm.nih.gov/geo/series/{series_prefix}/{gse_id}/suppl/"
            f"{gse_id}_Data.csv.gz"
        )
        with urllib.request.urlopen(suppl_url) as response:
            data = response.read()
            txt = gzip.decompress(data).decode("utf-8", errors="replace")
            df = pd.read_csv(io.StringIO(txt), index_col=0)

    if SAMPLE_COLS_DEFAULT and df.shape[1] > SAMPLE_COLS_DEFAULT:
        rng = np.random.RandomState(42)
        chosen = rng.choice(df.columns.to_numpy(), size=SAMPLE_COLS_DEFAULT, replace=False)
        df = df[chosen]

    if TOP_N_DEFAULT and df.shape[0] > TOP_N_DEFAULT:
        vars = df.var(axis=1)
        top_genes = vars.sort_values(ascending=False).head(TOP_N_DEFAULT).index
        df = df.loc[top_genes]

    df.to_csv(out_path)


def read_expression_matrix(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(
            f"Nu am găsit {path}. Puneți matricea de expresie la această locație."
        )
    df = pd.read_csv(path, index_col=0)
    if df.empty:
        raise ValueError("Matricea de expresie este goală.")
    return df


def log_and_filter(df: pd.DataFrame,
                   variance_threshold: float) -> pd.DataFrame:
    """
    Preprocesare:
    - aplică log2(x+1)
    - filtrează genele cu varianță scăzută
    """
    df_log = np.log2(df + 1)
    df_filt = df_log.loc[df_log.var(axis=1) > variance_threshold]
    return df_filt


def correlation_matrix(df: pd.DataFrame,
                       method: str = "spearman",
                       use_abs: bool = True) -> pd.DataFrame:
    """
    Hint:
      - df este (gene x probe); pentru corelație între gene, folosiți df.T.corr(method=...)
      - dacă use_abs=True, întoarceți |cor|
    """
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
        A[A < threshold] = 0
    else:
        A = (corr >= threshold).astype(int)
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
    Variante:
      - încercați louvain_communities(G, seed=42) dacă e disponibil
      - altfel greedy_modularity_communities(G)
    """
    try:
        from networkx.algorithms.community import louvain_communities
        communities = louvain_communities(G, seed=42)
    except Exception:
        from networkx.algorithms.community import greedy_modularity_communities
        communities_iterable: Iterable[Iterable[str]] = greedy_modularity_communities(G)
        communities = [set(c) for c in communities_iterable]

    mapping: Dict[str, int] = {}
    for midx, comm in enumerate(communities, start=1):
        for gene in comm:
            mapping[gene] = midx
    return mapping


def save_modules_csv(mapping: Dict[str, int], out_csv: Path) -> None:
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    df_modules = (
        pd.DataFrame({"Gene": list(mapping.keys()), "Module": list(mapping.values())})
        .sort_values(["Module", "Gene"])
    )
    df_modules.to_csv(out_csv, index=False)


if __name__ == "__main__":
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    download_matrix("GSE115469", INPUT_CSV)

    expr = read_expression_matrix(INPUT_CSV)
    expr_pp = log_and_filter(expr, variance_threshold=VARIANCE_THRESHOLD)

    corr = correlation_matrix(expr_pp, method=CORR_METHOD, use_abs=USE_ABS_CORR)
    adj = adjacency_from_correlation(corr, threshold=ADJ_THRESHOLD, weighted=False)

    G = graph_from_adjacency(adj, undirected=MAKE_UNDIRECTED)
    print(f"Grafic creat cu {G.number_of_nodes()} noduri și {G.number_of_edges()} muchii.")

    gene_to_module = detect_modules_louvain_or_greedy(G)
    print(f"S-au detectat {len(set(gene_to_module.values()))} module.")

    save_modules_csv(gene_to_module, OUTPUT_CSV)
    print(f"Am salvat mapping-ul gene→modul în: {OUTPUT_CSV}")
