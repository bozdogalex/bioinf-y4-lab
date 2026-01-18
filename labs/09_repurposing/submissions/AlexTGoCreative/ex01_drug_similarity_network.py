"""
Exercise 9.1 — Drug–Gene Bipartite Network & Drug Similarity Network

Scop:
- să construiți o rețea bipartită drug–gene plecând de la un CSV
- să proiectați layer-ul de medicamente folosind similaritatea dintre seturile de gene
- să exportați un fișier cu muchiile de similaritate între medicamente

TODO:
- încărcați datele drug–gene
- construiți dict-ul drug -> set de gene țintă
- construiți graful bipartit drug–gene (NetworkX)
- calculați similaritatea dintre medicamente (ex. Jaccard)
- construiți graful drug similarity
- exportați tabelul cu muchii: drug1, drug2, weight
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
HANDLE = "AlexTGoCreative"

# Input: fișier cu coloane cel puțin: drug, gene
DRUG_GENE_CSV = Path(f"data/work/{HANDLE}/lab09/drug_gene_{HANDLE}.csv")

# Output directory & files
OUT_DIR = Path(f"labs/09_repurposing/submissions/{HANDLE}")
OUT_DIR.mkdir(parents=True, exist_ok=True)

OUT_DRUG_SUMMARY = OUT_DIR / f"drug_summary_{HANDLE}.csv"
OUT_DRUG_SIMILARITY = OUT_DIR / f"drug_similarity_{HANDLE}.csv"
OUT_GRAPH_DRUG_GENE = OUT_DIR / f"network_drug_gene_{HANDLE}.gpickle"


def ensure_exists(path: Path) -> None:
    """
    TODO:
    - verificați că fișierul există
    - dacă nu, ridicați FileNotFoundError cu un mesaj clar
    """
    if not path.exists():
        raise FileNotFoundError(f"File not found: {path}")


def load_drug_gene_table(path: Path) -> pd.DataFrame:
    """
    TODO:
    - citiți CSV-ul cu pandas
    - validați că există cel puțin coloanele: 'drug', 'gene'
    - returnați DataFrame-ul
    """
    df = pd.read_csv(path, comment='#')
    
    if 'drug' not in df.columns or 'gene' not in df.columns:
        raise ValueError("CSV must contain 'drug' and 'gene' columns")
    
    return df


def build_drug2genes(df: pd.DataFrame) -> Dict[str, Set[str]]:
    """
    TODO:
    - construiți un dict: drug -> set de gene țintă
    - sugestie: folosiți groupby("drug") și aplicați set() pe coloana gene
    """
    drug2genes = {}
    for drug, group in df.groupby('drug'):
        drug2genes[drug] = set(group['gene'].values)
    return drug2genes


def build_bipartite_graph(drug2genes: Dict[str, Set[str]]) -> nx.Graph:
    """
    TODO:
    - construiți graful bipartit:
      - nodurile 'drug' cu atribut bipartite="drug"
      - nodurile 'gene' cu atribut bipartite="gene"
      - muchii drug-gene
    """
    B = nx.Graph()
    
    # Adăugăm nodurile de tip drug
    for drug in drug2genes.keys():
        B.add_node(drug, bipartite="drug")
    
    # Adăugăm nodurile de tip gene și muchiile
    for drug, genes in drug2genes.items():
        for gene in genes:
            if not B.has_node(gene):
                B.add_node(gene, bipartite="gene")
            B.add_edge(drug, gene)
    
    return B


def summarize_drugs(drug2genes: Dict[str, Set[str]]) -> pd.DataFrame:
    """
    TODO:
    - construiți un DataFrame cu:
        drug, num_targets (numărul de gene țintă)
    - returnați DataFrame-ul
    """
    data = []
    for drug, genes in drug2genes.items():
        data.append({'drug': drug, 'num_targets': len(genes)})
    
    df = pd.DataFrame(data)
    df = df.sort_values('num_targets', ascending=False)
    return df


def jaccard_similarity(s1: Set[str], s2: Set[str]) -> float:
    """
    Calculați similaritatea Jaccard între două seturi de gene:
    J(A, B) = |A ∩ B| / |A ∪ B|
    """
    if not s1 and not s2:
        return 0.0
    inter = len(s1 & s2)
    union = len(s1 | s2)
    return inter / union if union > 0 else 0.0


def compute_drug_similarity_edges(
    drug2genes: Dict[str, Set[str]],
    min_sim: float = 0.0,
) -> List[Tuple[str, str, float]]:
    """
    TODO:
    - pentru toate perechile de medicamente (combinații de câte 2),
      calculați similaritatea Jaccard între seturile de gene
    - rețineți doar muchiile cu similaritate >= min_sim
    - returnați o listă de tuple (drug1, drug2, weight)
    """
    edges = []
    drugs = list(drug2genes.keys())
    
    for drug1, drug2 in itertools.combinations(drugs, 2):
        sim = jaccard_similarity(drug2genes[drug1], drug2genes[drug2])
        if sim >= min_sim:
            edges.append((drug1, drug2, sim))
    
    return edges


def edges_to_dataframe(edges: List[Tuple[str, str, float]]) -> pd.DataFrame:
    """
    TODO:
    - transformați lista de muchii (drug1, drug2, weight) într-un DataFrame
      cu coloanele: drug1, drug2, similarity
    """
    if not edges:
        return pd.DataFrame(columns=['drug1', 'drug2', 'similarity'])
    
    df = pd.DataFrame(edges, columns=['drug1', 'drug2', 'similarity'])
    df = df.sort_values('similarity', ascending=False)
    return df


# --------------------------
# Main
# --------------------------
if __name__ == "__main__":
    print(f"[INFO] Starting Exercise 9.1 for {HANDLE}")
    
    # TODO 1: verificați că fișierul de input există
    print(f"[INFO] Checking input file: {DRUG_GENE_CSV}")
    ensure_exists(DRUG_GENE_CSV)
    
    # TODO 2: încărcați tabelul drug-gene
    print("[INFO] Loading drug-gene table...")
    df_drug_gene = load_drug_gene_table(DRUG_GENE_CSV)
    print(f"[INFO] Loaded {len(df_drug_gene)} drug-gene interactions")
    
    # TODO 3: construiți mapping-ul drug -> set de gene
    print("[INFO] Building drug -> genes mapping...")
    drug2genes = build_drug2genes(df_drug_gene)
    print(f"[INFO] Found {len(drug2genes)} unique drugs")
    
    # TODO 4: construiți graful bipartit și salvați-l (opțional)
    print("[INFO] Building bipartite graph...")
    B = build_bipartite_graph(drug2genes)
    print(f"[INFO] Graph has {B.number_of_nodes()} nodes and {B.number_of_edges()} edges")
    
    # Salvăm graful
    import pickle
    with open(OUT_GRAPH_DRUG_GENE, 'wb') as f:
        pickle.dump(B, f)
    print(f"[INFO] Saved bipartite graph to {OUT_GRAPH_DRUG_GENE}")
    
    # TODO 5: generați și salvați sumarul pe medicamente
    print("[INFO] Generating drug summary...")
    df_summary = summarize_drugs(drug2genes)
    df_summary.to_csv(OUT_DRUG_SUMMARY, index=False)
    print(f"[INFO] Saved drug summary to {OUT_DRUG_SUMMARY}")
    print(f"\nDrug summary:\n{df_summary}")
    
    # TODO 6: calculați similaritatea între medicamente
    print("\n[INFO] Computing drug similarity...")
    edges = compute_drug_similarity_edges(drug2genes, min_sim=0.0)
    print(f"[INFO] Found {len(edges)} drug-drug similarity edges")
    
    df_similarity = edges_to_dataframe(edges)
    df_similarity.to_csv(OUT_DRUG_SIMILARITY, index=False)
    print(f"[INFO] Saved drug similarity to {OUT_DRUG_SIMILARITY}")
    print(f"\nTop drug similarities:\n{df_similarity.head(10)}")
    
    print("\n[SUCCESS] Exercise 9.1 completed successfully!")
