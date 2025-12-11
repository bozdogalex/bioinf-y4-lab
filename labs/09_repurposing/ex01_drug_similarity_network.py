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
HANDLE = "<handle>"

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
    pass


def load_drug_gene_table(path: Path) -> pd.DataFrame:
    """
    TODO:
    - citiți CSV-ul cu pandas
    - validați că există cel puțin coloanele: 'drug', 'gene'
    - returnați DataFrame-ul
    """

    return pd.DataFrame()  # placeholder


def build_drug2genes(df: pd.DataFrame) -> Dict[str, Set[str]]:
    """
    TODO:
    - construiți un dict: drug -> set de gene țintă
    - sugestie: folosiți groupby("drug") și aplicați set() pe coloana gene
    """

    return {}  # placeholder


def build_bipartite_graph(drug2genes: Dict[str, Set[str]]) -> nx.Graph:
    """
    TODO:
    - construiți graful bipartit:
      - nodurile 'drug' cu atribut bipartite="drug"
      - nodurile 'gene' cu atribut bipartite="gene"
      - muchii drug-gene
    """

    return nx.Graph()  # placeholder


def summarize_drugs(drug2genes: Dict[str, Set[str]]) -> pd.DataFrame:
    """
    TODO:
    - construiți un DataFrame cu:
        drug, num_targets (numărul de gene țintă)
    - returnați DataFrame-ul
    """

    return pd.DataFrame()  # placeholder


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

    return []  # placeholder


def edges_to_dataframe(edges: List[Tuple[str, str, float]]) -> pd.DataFrame:
    """
    TODO:
    - transformați lista de muchii (drug1, drug2, weight) într-un DataFrame
      cu coloanele: drug1, drug2, similarity
    """
    return pd.DataFrame()  # placeholder


# --------------------------
# Main
# --------------------------
if __name__ == "__main__":
    # TODO 1: verificați că fișierul de input există

    # TODO 2: încărcați tabelul drug-gene

    # TODO 3: construiți mapping-ul drug -> set de gene

    # TODO 4: construiți graful bipartit și salvați-l (opțional)

    # TODO 5: generați și salvați sumarul pe medicamente

    # TODO 6: calculați similaritatea între medicamente

    print("[INFO] Exercise 9.1 skeleton — completați TODO-urile și rulați scriptul.")
