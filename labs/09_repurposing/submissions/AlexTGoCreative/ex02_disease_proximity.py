"""
Exercise 9.2 — Disease Proximity and Drug Ranking

Scop:
- să calculați distanța medie dintre fiecare medicament și un set de gene asociate unei boli
- să ordonați medicamentele în funcție de proximitate (network-based prioritization)

TODO-uri principale:
- încărcați graful bipartit drug–gene (din exercițiul 9.1) sau reconstruiți-l
- încărcați lista de disease genes
- pentru fiecare medicament, calculați distanța minimă / medie până la genele bolii
- exportați un tabel cu medicamente și scorul lor de proximitate
"""

from __future__ import annotations
from pathlib import Path
from typing import Dict, Set, List, Tuple

import networkx as nx
import pandas as pd

# --------------------------
# Config
# --------------------------
HANDLE = "AlexTGoCreative"

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
    """
    TODO:
    - verificați că fișierul există
    - dacă nu, ridicați FileNotFoundError
    """
    if not path.exists():
        raise FileNotFoundError(f"File not found: {path}")


def load_bipartite_graph_or_build() -> nx.Graph:
    """
    TODO:
    - dacă GRAPH_DRUG_GENE există, încărcați-l direct
    - altfel, reconstruiți graful plecând de la DRUG_GENE_CSV
      (puteți reutiliza logica din ex09_drug_similarity_network.py)
    """
    if GRAPH_DRUG_GENE.exists():
        print(f"[INFO] Loading saved bipartite graph from {GRAPH_DRUG_GENE}")
        import pickle
        with open(GRAPH_DRUG_GENE, 'rb') as f:
            return pickle.load(f)
    
    print(f"[INFO] Building bipartite graph from {DRUG_GENE_CSV}")
    ensure_exists(DRUG_GENE_CSV)
    
    # Reconstruim graful
    df = pd.read_csv(DRUG_GENE_CSV, comment='#')
    
    # Construim mapping-ul drug -> gene set
    drug2genes = {}
    for drug, group in df.groupby('drug'):
        drug2genes[drug] = set(group['gene'].values)
    
    # Construim graful bipartit
    B = nx.Graph()
    for drug in drug2genes.keys():
        B.add_node(drug, bipartite="drug")
    
    for drug, genes in drug2genes.items():
        for gene in genes:
            if not B.has_node(gene):
                B.add_node(gene, bipartite="gene")
            B.add_edge(drug, gene)
    
    return B


def load_disease_genes(path: Path) -> Set[str]:
    """
    TODO:
    - încărcați fișierul text cu gene (una pe linie)
    - returnați un set de gene (string)
    """
    genes = set()
    with open(path, 'r') as f:
        for line in f:
            gene = line.strip()
            if gene and not gene.startswith('#'):
                genes.add(gene)
    return genes


def get_drug_nodes(B: nx.Graph) -> List[str]:
    """
    TODO:
    - extrageți lista nodurilor de tip 'drug'
    - presupunem atributul bipartite="drug"
    """
    drugs = [n for n, d in B.nodes(data=True) if d.get("bipartite") == "drug"]
    return drugs


def compute_drug_disease_distance(
    B: nx.Graph,
    drug: str,
    disease_genes: Set[str],
    mode: str = "mean",
    max_dist: int = 5,
) -> float:
    """
    TODO:
    - pentru un medicament:
      - calculați distanța (de ex. shortest_path_length) până la fiecare genă din disease_genes
      - ignorați genele care nu sunt în graf
      - dacă nu există niciun drum, puteți seta o distanță penalizantă (ex. max_dist + 1)
    - returnați media (sau minimul) distanțelor; controlați cu parametrul 'mode'
    """
    distances = []
    
    for gene in disease_genes:
        if gene not in B:
            # Gena nu e în graf, o ignorăm
            continue
        
        try:
            dist = nx.shortest_path_length(B, source=drug, target=gene)
            distances.append(dist)
        except nx.NetworkXNoPath:
            # Nu există drum între drug și gene
            distances.append(max_dist + 1)
    
    if not distances:
        return float('inf')
    
    if mode == "mean":
        return sum(distances) / len(distances)
    elif mode == "min":
        return min(distances)
    elif mode == "median":
        distances_sorted = sorted(distances)
        n = len(distances_sorted)
        if n % 2 == 0:
            return (distances_sorted[n//2 - 1] + distances_sorted[n//2]) / 2
        else:
            return distances_sorted[n//2]
    else:
        return sum(distances) / len(distances)


def rank_drugs_by_proximity(
    B: nx.Graph,
    disease_genes: Set[str],
    mode: str = "mean",
) -> pd.DataFrame:
    """
    TODO:
    - pentru fiecare medicament din graf:
      - calculați scorul de distanță (ex. media distanțelor către genele bolii)
    - construiți un DataFrame cu:
      drug, distance
    - sortați crescător după distance (distanță mai mică = proximitate mai mare)
    """
    drugs = get_drug_nodes(B)
    
    results = []
    for drug in drugs:
        dist = compute_drug_disease_distance(B, drug, disease_genes, mode=mode)
        results.append({'drug': drug, 'distance': dist})
    
    df = pd.DataFrame(results)
    df = df.sort_values('distance')
    
    return df


# --------------------------
# Main
# --------------------------
if __name__ == "__main__":
    print(f"[INFO] Starting Exercise 9.2 for {HANDLE}")
    
    # TODO 1: verificați input-urile
    print(f"[INFO] Checking disease genes file: {DISEASE_GENES_TXT}")
    ensure_exists(DISEASE_GENES_TXT)
    
    # TODO 2: încărcați / construiți graful bipartit
    print("[INFO] Loading bipartite graph...")
    B = load_bipartite_graph_or_build()
    print(f"[INFO] Graph has {B.number_of_nodes()} nodes and {B.number_of_edges()} edges")
    
    # TODO 3: încărcați setul de disease genes
    print("[INFO] Loading disease genes...")
    disease_genes = load_disease_genes(DISEASE_GENES_TXT)
    print(f"[INFO] Loaded {len(disease_genes)} disease genes: {disease_genes}")
    
    # TODO 4: calculați ranking-ul medicamentelor după proximitate
    print("[INFO] Computing drug-disease proximity...")
    df_priority = rank_drugs_by_proximity(B, disease_genes, mode="mean")
    print(f"[INFO] Ranked {len(df_priority)} drugs")
    
    # TODO 5: salvați rezultatele
    df_priority.to_csv(OUT_DRUG_PRIORITY, index=False)
    print(f"[INFO] Saved drug priority ranking to {OUT_DRUG_PRIORITY}")
    print(f"\nDrug priority ranking (top candidates have lower distance):\n{df_priority}")
    
    print("\n[SUCCESS] Exercise 9.2 completed successfully!")
