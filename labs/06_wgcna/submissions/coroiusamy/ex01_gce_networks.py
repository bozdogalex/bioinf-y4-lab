#!/usr/bin/env python
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
import urllib.request
import gzip
import shutil
import os
import numpy as np
import pandas as pd
import networkx as nx


# === CONFIGURARE UTILIZATOR === #
HANDLE = "coroiusamy"

INPUT_CSV = Path(f"data/work/{HANDLE}/lab06/expression_matrix.csv")
OUTPUT_DIR = Path(f"labs/06_wgcna/submissions/{HANDLE}")
OUTPUT_CSV = OUTPUT_DIR / f"modules_{HANDLE}.csv"

# Parametrii pentru construcția rețelei
CORR_METHOD = "spearman"     # "pearson" sau "spearman"
VARIANCE_THRESHOLD = 0.1    # prag pentru filtrare gene
ADJ_THRESHOLD = 0.6          # prag pentru |cor| (ex: 0.6)
USE_ABS_CORR = True          # True => folosiți |cor| la prag
MAKE_UNDIRECTED = True       # rețelele de co-expresie sunt neorientate


def setup_data_file(target_csv: Path) -> None:
    """
    Verifică dacă fișierul CSV există.
    Dacă nu, descarcă și decomprimă arhiva .gz de pe GEO.
    """
    if target_csv.exists():
        print(f"Fișierul de date a fost găsit: {target_csv}")
        return

    print(f"Fișierul {target_csv} nu a fost găsit. Se încearcă descărcarea...")

    data_url = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE115469&format=file&file=GSE115469_Data.csv.gz"
    gz_path = target_csv.parent / "GSE115469_Data.csv.gz"
    target_csv.parent.mkdir(parents=True, exist_ok=True)

    try:
        # 1. Descarcare fisier .gz
        print(f"Informație: Se descarcă {data_url}...")
        with urllib.request.urlopen(data_url) as response, open(gz_path, 'wb') as out_file:
            shutil.copyfileobj(response, out_file)
        print(f"Descărcare finalizată: {gz_path}")

        # 2. Decomprimare fișier .gz
        print(f"Informație: Se decomprimă {gz_path} în {target_csv}...")
        with gzip.open(gz_path, 'rb') as f_in:
            with open(target_csv, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        print("Decomprimare finalizată.")

        # 3. stergere arhiva temporara
        os.remove(gz_path)
        print(f"Arhiva {gz_path} a fost ștearsă.")
        print(f"Datele sunt gata la: {target_csv}")

    except Exception as e:
        print(f"Eroare la descărcarea sau decomprimarea datelor: {e}")
        print("Vă rugăm să descărcați și să plasați manual fișierul.")
        raise


def read_expression_matrix(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Nu am găsit {path}. Rulați din nou scriptul pentru a încerca descărcarea.")

    # Definim un numar maxim de rânduri (gene) de încarcat
    MAX_GENES_TO_LOAD = 10000
    
    print(f"Se încarcă primele {MAX_GENES_TO_LOAD} gene din {path}...")
    
    df = pd.read_csv(path, index_col=0, nrows=MAX_GENES_TO_LOAD)

    if df.empty:
        raise ValueError("Matricea de expresie este goală.")
    
    print(f"Datele au fost încărcate: {df.shape[0]} gene, {df.shape[1]} probe.")
    return df


def log_and_filter(df: pd.DataFrame, variance_threshold: float) -> pd.DataFrame:
    """Preprocesare: log2(x+1) + filtrare gene cu varianță scăzută"""
    df_log = np.log2(df + 1)
    df_filt = df_log.loc[df_log.var(axis=1) > variance_threshold]
    return df_filt


def correlation_matrix(df: pd.DataFrame, method: str = "spearman", use_abs: bool = True) -> pd.DataFrame:
    """Calculează matricea de corelație între gene (rânduri)."""
    corr = df.T.corr(method=method)
    return corr.abs() if use_abs else corr


def adjacency_from_correlation(corr: pd.DataFrame, threshold: float, weighted: bool = False) -> pd.DataFrame:
    """Construiește matricea de adiacență din corelații."""
    if weighted:
        A = corr.copy()
        A[A < threshold] = 0
    else:
        A = (corr >= threshold).astype(int)

    np.fill_diagonal(A.values, 0.0)
    return A


def graph_from_adjacency(A: pd.DataFrame, undirected: bool = True) -> nx.Graph:
    """Creează un graf NetworkX din matricea de adiacență."""
    if undirected:
        G = nx.from_pandas_adjacency(A)
    else:
        G = nx.from_pandas_adjacency(A, create_using=nx.DiGraph)

    isolates = list(nx.isolates(G))
    if isolates:
        G.remove_nodes_from(isolates)

    return G


def detect_modules_louvain_or_greedy(G: nx.Graph) -> Dict[str, int]:
    """Detectează comunități (module) și întoarce un dicționar gene → modul_id."""
    communities = []
    try:
        from networkx.algorithms.community import louvain_communities
        communities_iterable = louvain_communities(G, seed=42)
        communities = [set(c) for c in communities_iterable]
        print("Am folosit algoritmul Louvain pentru detecția modulelor.")

    except ImportError:
        print("AVERTISMENT: 'community' nu e instalat. Se folosește fallback: greedy_modularity_communities.")
        from networkx.algorithms.community import greedy_modularity_communities
        communities_iterable = greedy_modularity_communities(G)
        communities = [set(c) for c in communities_iterable]

    except Exception as e:
        print(f"Eroare la detecția comunităților ({e}). Se folosește fallback pe greedy.")
        from networkx.algorithms.community import greedy_modularity_communities
        communities_iterable = greedy_modularity_communities(G)
        communities = [set(c) for c in communities_iterable]

    mapping = {gene: midx for midx, comm in enumerate(communities, start=1) for gene in comm}
    return mapping


def save_modules_csv(mapping: Dict[str, int], out_csv: Path) -> None:
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    df_modules = (
        pd.DataFrame({"Gene": list(mapping.keys()), "Module": list(mapping.values())})
        .sort_values(["Module", "Gene"])
    )
    df_modules.to_csv(out_csv, index=False)
    print(f"Am salvat mapping-ul gene→modul în: {out_csv}")


if __name__ == "__main__":
    setup_data_file(INPUT_CSV)

    expr = read_expression_matrix(INPUT_CSV)
    expr_pp = log_and_filter(expr, variance_threshold=VARIANCE_THRESHOLD)
    #reducere nr gene
    expr_pp = expr_pp.sample(n=500, random_state=42)
    print(f"Folosim {expr_pp.shape[0]} gene pentru construcția rețelei.")

    corr = correlation_matrix(expr_pp, method=CORR_METHOD, use_abs=USE_ABS_CORR)
    adj = adjacency_from_correlation(corr, threshold=ADJ_THRESHOLD, weighted=False)

    G = graph_from_adjacency(adj, undirected=MAKE_UNDIRECTED)
    print(f"Grafic creat cu {G.number_of_nodes()} noduri și {G.number_of_edges()} muchii.")

    gene_to_module = detect_modules_louvain_or_greedy(G)
    print(f"S-au detectat {len(set(gene_to_module.values()))} module.")

    save_modules_csv(gene_to_module, OUTPUT_CSV)
    print(f"Am salvat mapping-ul gene→modul în: {OUTPUT_CSV}")
