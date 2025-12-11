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
HANDLE = "<handle>"

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

    pass


def load_bipartite_graph_or_build() -> nx.Graph:
    """
    TODO:
    - dacă GRAPH_DRUG_GENE există, încărcați-l direct
    - altfel, reconstruiți graful plecând de la DRUG_GENE_CSV
      (puteți reutiliza logica din ex09_drug_similarity_network.py)
    """

    return nx.Graph()  # placeholder


def load_disease_genes(path: Path) -> Set[str]:
    """
    TODO:
    - încărcați fișierul text cu gene (una pe linie)
    - returnați un set de gene (string)
    """

    return set()  # placeholder


def get_drug_nodes(B: nx.Graph) -> List[str]:
    """
    TODO:
    - extrageți lista nodurilor de tip 'drug'
    - presupunem atributul bipartite="drug"
    """
    # drugs = [n for n, d in B.nodes(data=True) if d.get("bipartite") == "drug"]
    # return drugs
    return []  # placeholder


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

    return 0.0  # placeholder


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

    return pd.DataFrame()  # placeholder


# --------------------------
# Main
# --------------------------
if __name__ == "__main__":
    # TODO 1: verificați input-urile

    # TODO 2: încărcați / construiți graful bipartit

    # TODO 3: încărcați setul de disease genes

    # TODO 4: calculați ranking-ul medicamentelor după proximitate

    # TODO 5: salvați rezultatele

    print("[INFO] Exercise 9.2 skeleton — completați TODO-urile și rulați scriptul.")
