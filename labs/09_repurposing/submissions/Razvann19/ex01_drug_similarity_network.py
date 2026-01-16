"""
Exercise 9.1 — Drug–Gene Bipartite Network & Drug Similarity Network

Rezolvă:
- Task 1: bipartite graph + drug_summary_<handle>.csv
- Task 2: drug similarity (Jaccard) + drug_similarity_<handle>.csv
"""

from __future__ import annotations
from pathlib import Path
from typing import Dict, Set, Tuple, List


import itertools

import networkx as nx
import pandas as pd


HANDLE = "Razvann19"

DRUG_GENE_CSV = Path(__file__).parent / f"drug_gene_{HANDLE}.csv"

OUT_DIR = Path(f"labs/09_repurposing/submissions/{HANDLE}")
OUT_DIR.mkdir(parents=True, exist_ok=True)

OUT_DRUG_SUMMARY = OUT_DIR / f"drug_summary_{HANDLE}.csv"
OUT_DRUG_SIMILARITY = OUT_DIR / f"drug_similarity_{HANDLE}.csv"


def ensure_exists(path: Path) -> None:
    """Verifică existența fișierului și dă o eroare clară dacă lipsește."""
    if not path.exists():
        raise FileNotFoundError(
            f"[ERROR] Nu găsesc fișierul de input: {path}\n"
            f"Verifică HANDLE-ul și path-ul din DRUG_GENE_CSV."
        )
    if not path.is_file():
        raise FileNotFoundError(f"[ERROR] Path-ul există dar nu e fișier: {path}")


def load_drug_gene_table(path: Path) -> pd.DataFrame:
    """
    - citește CSV-ul
    - validează coloanele minime: 'drug', 'gene'
    - curăță rânduri invalide (NaN / string gol)
    """
    df = pd.read_csv(path)

    required = {"drug", "gene"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(
            f"[ERROR] CSV-ul nu conține coloanele necesare {required}. "
            f"Lipsesc: {missing}. Coloane găsite: {list(df.columns)}"
        )

    df = df[["drug", "gene"]].copy()
    df["drug"] = df["drug"].astype(str).str.strip()
    df["gene"] = df["gene"].astype(str).str.strip()
    df = df[(df["drug"] != "") & (df["gene"] != "")]
    df = df.dropna(subset=["drug", "gene"])

    return df


def build_drug2genes(df: pd.DataFrame) -> Dict[str, Set[str]]:
    """Construiește mapping drug -> set(gene)."""
    drug2genes = (
        df.groupby("drug")["gene"]
        .apply(lambda s: set(s.dropna().astype(str)))
        .to_dict()
    )
    drug2genes = {d: gs for d, gs in drug2genes.items() if gs}
    return drug2genes


def build_bipartite_graph(drug2genes: Dict[str, Set[str]]) -> nx.Graph:
    B = nx.Graph()

    for drug, genes in drug2genes.items():
        B.add_node(drug, bipartite="drug", type="drug")
        for gene in genes:
            B.add_node(gene, bipartite="gene", type="gene")
            B.add_edge(drug, gene)

    return B


def summarize_drugs(drug2genes: Dict[str, Set[str]]) -> pd.DataFrame:
    """Returnează DataFrame cu: drug, num_targets."""
    rows = [{"drug": drug, "num_targets": len(genes)} for drug, genes in drug2genes.items()]
    out = pd.DataFrame(rows).sort_values(["num_targets", "drug"], ascending=[False, True])
    return out


def jaccard_similarity(s1: Set[str], s2: Set[str]) -> float:
    """J(A, B) = |A ∩ B| / |A ∪ B|"""
    if not s1 and not s2:
        return 0.0
    inter = len(s1 & s2)
    union = len(s1 | s2)
    return inter / union if union > 0 else 0.0


def compute_drug_similarity_edges(
    drug2genes: Dict[str, Set[str]],
    min_sim: float = 0.0,
) -> List[Tuple[str, str, float]]:
    drugs = sorted(drug2genes.keys())
    edges: List[Tuple[str, str, float]] = []

    for d1, d2 in itertools.combinations(drugs, 2):
        sim = jaccard_similarity(drug2genes[d1], drug2genes[d2])
        if sim >= min_sim:
            edges.append((d1, d2, float(sim)))

    return edges


def edges_to_dataframe(edges: List[Tuple[str, str, float]]) -> pd.DataFrame:
    """Transformă muchiile în DataFrame: drug1, drug2, similarity."""
    df = pd.DataFrame(edges, columns=["drug1", "drug2", "similarity"])
    df = df.sort_values(["similarity", "drug1", "drug2"], ascending=[False, True, True])
    return df


if __name__ == "__main__":
    ensure_exists(DRUG_GENE_CSV)
    df = load_drug_gene_table(DRUG_GENE_CSV)

    drug2genes = build_drug2genes(df)
    if not drug2genes:
        raise ValueError("[ERROR] Nu am putut construi drug2genes (dataset gol sau invalid).")

    B = build_bipartite_graph(drug2genes)

    drug_summary = summarize_drugs(drug2genes)
    drug_summary.to_csv(OUT_DRUG_SUMMARY, index=False)


    edges = compute_drug_similarity_edges(drug2genes, min_sim=0.0)
    sim_df = edges_to_dataframe(edges)
    sim_df.to_csv(OUT_DRUG_SIMILARITY, index=False)

    print("[INFO] Gata ")
    print(f"  - {OUT_DRUG_SUMMARY}")
    print(f"  - {OUT_DRUG_SIMILARITY}")
