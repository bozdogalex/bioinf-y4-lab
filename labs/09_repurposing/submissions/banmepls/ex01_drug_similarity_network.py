"""
Exercise 9.1 — Drug–Gene Bipartite Network & Drug Similarity Network
"""

from __future__ import annotations
from pathlib import Path
from typing import Dict, Set, Tuple, List

import itertools

import networkx as nx
import pandas as pd

# --------------------------
# Config — adaptați pentru handle-ul vostru
# --------------------------
HANDLE = "banmepls"

# Input: fișier cu coloane cel puțin: drug, gene
DRUG_GENE_CSV = Path(f"data/work/{HANDLE}/lab09/drug_gene_{HANDLE}.csv")

# Output directory & files
OUT_DIR = Path(f"labs/09_repurposing/submissions/{HANDLE}")
OUT_DIR.mkdir(parents=True, exist_ok=True)

OUT_DRUG_SUMMARY = OUT_DIR / f"drug_summary_{HANDLE}.csv"
OUT_DRUG_SIMILARITY = OUT_DIR / f"drug_similarity_{HANDLE}.csv"
OUT_GRAPH_DRUG_GENE = OUT_DIR / f"network_drug_gene_{HANDLE}.gpickle"


def ensure_exists(path: Path) -> None:
    if not path.is_file():
        raise FileNotFoundError(
            f"Fișierul de date lipsește: {path}\n"
            f"Verifică dacă ai descărcat și salvat tabelul drug-gene corect."
        )


def load_drug_gene_table(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    
    # Validare minimă
    required_cols = {"drug", "gene"}
    if not required_cols.issubset(df.columns):
        raise ValueError(f"Fișierul CSV trebuie să conțină cel puțin coloanele: {required_cols}")
    
    return df


def build_drug2genes(df: pd.DataFrame) -> Dict[str, Set[str]]:
    # Grupăm după drug și transformăm genele în set unic
    drug2genes = df.groupby("drug")["gene"].apply(set).to_dict()
    return drug2genes


def build_bipartite_graph(drug2genes: Dict[str, Set[str]]) -> nx.Graph:
    G = nx.Graph()
    
    # Noduri cu bipartiție
    drugs = list(drug2genes.keys())
    genes = set()
    for geneset in drug2genes.values():
        genes.update(geneset)
    
    G.add_nodes_from(drugs, bipartite="drug")
    G.add_nodes_from(genes, bipartite="gene")
    
    # Muchii drug → gene
    edges = []
    for drug, geneset in drug2genes.items():
        for gene in geneset:
            edges.append((drug, gene))
    G.add_edges_from(edges)
    
    print(f"   → {G.number_of_edges()} interacțiuni")
    
    return G


def summarize_drugs(drug2genes: Dict[str, Set[str]]) -> pd.DataFrame:
    data = []
    for drug, geneset in drug2genes.items():
        data.append({
            "drug": drug,
            "num_targets": len(geneset)
        })
    
    df_summary = pd.DataFrame(data).sort_values("num_targets", ascending=False)
    return df_summary


def jaccard_similarity(s1: Set[str], s2: Set[str]) -> float:
    if not s1 and not s2:
        return 0.0
    inter = len(s1 & s2)
    union = len(s1 | s2)
    return inter / union if union > 0 else 0.0


def compute_drug_similarity_edges(
    drug2genes: Dict[str, Set[str]],
    min_sim: float = 0.0,
) -> List[Tuple[str, str, float]]:
    edges = []
    
    for d1, d2 in itertools.combinations(drug2genes.keys(), 2):
        sim = jaccard_similarity(drug2genes[d1], drug2genes[d2])
        if sim >= min_sim:
            edges.append((d1, d2, sim))
    
    return edges


def edges_to_dataframe(edges: List[Tuple[str, str, float]]) -> pd.DataFrame:
    df = pd.DataFrame(edges, columns=["drug1", "drug2", "similarity"])
    df = df.sort_values("similarity", ascending=False)
    return df


# --------------------------
# Main
# --------------------------
if __name__ == "__main__":
    ensure_exists(DRUG_GENE_CSV)
    
    df_drug_gene = load_drug_gene_table(DRUG_GENE_CSV)
    
    drug2genes = build_drug2genes(df_drug_gene)
    
    G_bipartite = build_bipartite_graph(drug2genes)
    
    df_summary = summarize_drugs(drug2genes)
    df_summary.to_csv(OUT_DRUG_SUMMARY, index=False)
    print(f"   → {len(df_summary)} medicamente salvate în {OUT_DRUG_SUMMARY}")
    
    sim_edges = compute_drug_similarity_edges(drug2genes, min_sim=0.1)
    
    df_sim = edges_to_dataframe(sim_edges)
    df_sim.to_csv(OUT_DRUG_SIMILARITY, index=False)
    print(f"   → {len(df_sim)} muchii salvate în {OUT_DRUG_SIMILARITY}")
