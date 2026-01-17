"""
Exercise 9.2 — Disease Proximity and Drug Ranking
"""

from __future__ import annotations
from pathlib import Path
from typing import Dict, Set, List, Tuple

import networkx as nx
import pandas as pd

# --------------------------
# Config
# --------------------------
HANDLE = "banmepls"

# Input: graful bipartit (salvat anterior) SAU tabelul drug-gene
GRAPH_DRUG_GENE = Path(f"labs/09_repurposing/submissions/{HANDLE}/network_drug_gene_{HANDLE}.gpickle")
DRUG_GENE_CSV = Path(f"data/work/{HANDLE}/lab09/drug_gene_{HANDLE}.csv")

# Input: lista genelor bolii
DISEASE_GENES_TXT = Path(f"data/work/{HANDLE}/lab09/disease_genes_{HANDLE}.txt")

# Output directory & file
OUT_DIR = Path(f"labs/09_repurposing/submissions/{HANDLE}")
OUT_DIR.mkdir(parents=True, exist_ok=True)

OUT_DRUG_PRIORITY = OUT_DIR / f"drug_priority_{HANDLE}.csv"


# --------------------------
# Utils
# --------------------------
def ensure_exists(path: Path) -> None:
    if not path.is_file():
        raise FileNotFoundError(f"Fișierul lipsește: {path}")


def load_bipartite_graph_or_build() -> nx.Graph:
    if GRAPH_DRUG_GENE.exists():
        return nx.read_graphml(GRAPH_DRUG_GENE)
    
    print("[INFO] Reconstruire graf bipartit din tabel...")
    df = pd.read_csv(DRUG_GENE_CSV)
    required = {"drug", "gene"}
    if not required.issubset(df.columns):
        raise ValueError(f"CSV-ul trebuie să aibă cel puțin coloanele: {required}")
    
    G = nx.Graph()
    drugs = df["drug"].unique()
    genes = df["gene"].unique()
    
    G.add_nodes_from(drugs, bipartite="drug")
    G.add_nodes_from(genes, bipartite="gene")
    
    for _, row in df.iterrows():
        G.add_edge(row["drug"], row["gene"])
    
    # Salvăm ca GraphML (format standard, ușor de deschis)
    nx.write_graphml(G, GRAPH_DRUG_GENE)
    return G


def load_disease_genes(path: Path) -> Set[str]:
    if not path.is_file():
        return {"TP53", "MDM2", "CDKN1A", "BAX", "PUMA", "NOXA"}
    
    with open(path, "r") as f:
        genes = {line.strip() for line in f if line.strip()}
    return genes


def get_drug_nodes(B: nx.Graph) -> List[str]:
    drugs = [n for n, d in B.nodes(data=True) if d.get("bipartite") == "drug"]
    return drugs


def compute_drug_disease_distance(
    B: nx.Graph,
    drug: str,
    disease_genes: Set[str],
    mode: str = "mean",
    max_dist: int = 5,
) -> float:
    distances = []
    for gene in disease_genes:
        if gene not in B:
            continue
        try:
            dist = nx.shortest_path_length(B, drug, gene)
            distances.append(dist)
        except nx.NetworkXNoPath:
            distances.append(max_dist + 1)  # penalizare
    
    if not distances:
        return max_dist + 1
    
    if mode == "mean":
        return np.mean(distances)
    elif mode == "min":
        return np.min(distances)
    else:
        return np.mean(distances)


def rank_drugs_by_proximity(
    B: nx.Graph,
    disease_genes: Set[str],
    mode: str = "mean",
) -> pd.DataFrame:
    drugs = get_drug_nodes(B)
    
    data = []
    for drug in drugs:
        dist = compute_drug_disease_distance(B, drug, disease_genes, mode=mode)
        data.append({"drug": drug, "distance": dist})
    
    df_rank = pd.DataFrame(data).sort_values("distance")
    return df_rank


# --------------------------
# Main
# --------------------------
if __name__ == "__main__":
    ensure_exists(DRUG_GENE_CSV)
    
    B = load_bipartite_graph_or_build()
    
    disease_genes = load_disease_genes(DISEASE_GENES_TXT)
    
    df_priority = rank_drugs_by_proximity(B, disease_genes, mode="mean")
    
    df_priority.to_csv(OUT_DRUG_PRIORITY, index=False)
    print("Top 10 medicamente cele mai apropiate de genele bolii:")
    print(df_priority.head(10))
