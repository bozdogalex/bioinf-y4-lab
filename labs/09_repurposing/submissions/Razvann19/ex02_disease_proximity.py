"""
Exercise 9.2 — Disease Proximity and Drug Ranking (+ Vizualizare Task 4)

Rezolvă:
- Task 3: drug_priority_<handle>.csv (ranking după distanța medie drug ↔ disease genes)
- Task 4: network_drug_gene_<handle>.png (bipartite graph, color + size)

Input așteptat (în același folder cu scriptul):
- drug_gene_<handle>.csv   (coloane: drug, gene, ... optional)
- disease_genes_<handle>.txt (o genă pe linie)
"""

from __future__ import annotations
from pathlib import Path
from typing import Dict, Set, List, Tuple

import networkx as nx
import pandas as pd
import matplotlib.pyplot as plt


HANDLE = "Razvann19"

HERE = Path(__file__).parent
DRUG_GENE_CSV = HERE / f"drug_gene_{HANDLE}.csv"
DISEASE_GENES_TXT = HERE / f"disease_genes_{HANDLE}.txt"

OUT_DIR = Path(f"labs/09_repurposing/submissions/{HANDLE}")
OUT_DIR.mkdir(parents=True, exist_ok=True)

OUT_DRUG_PRIORITY = OUT_DIR / f"drug_priority_{HANDLE}.csv"
OUT_NETWORK_PNG = OUT_DIR / f"network_drug_gene_{HANDLE}.png"


def ensure_exists(path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(f"[ERROR] Nu găsesc fișierul: {path}")
    if not path.is_file():
        raise FileNotFoundError(f"[ERROR] Path-ul există dar nu e fișier: {path}")


def load_drug_gene_table(path: Path) -> pd.DataFrame:
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
    drug2genes = (
        df.groupby("drug")["gene"]
        .apply(lambda s: set(s.dropna().astype(str)))
        .to_dict()
    )
    return {d: gs for d, gs in drug2genes.items() if gs}


def build_bipartite_graph(drug2genes: Dict[str, Set[str]]) -> nx.Graph:
    B = nx.Graph()
    for drug, genes in drug2genes.items():
        B.add_node(drug, bipartite="drug", type="drug")
        for gene in genes:
            B.add_node(gene, bipartite="gene", type="gene")
            B.add_edge(drug, gene)
    return B


def load_disease_genes(path: Path) -> Set[str]:
    genes: Set[str] = set()
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            g = line.strip()
            if g:
                genes.add(g)
    return genes


def get_drug_nodes(B: nx.Graph) -> List[str]:
    return [n for n, d in B.nodes(data=True) if d.get("bipartite") == "drug"]


def compute_drug_disease_distance(
    B: nx.Graph,
    drug: str,
    disease_genes: Set[str],
    mode: str = "mean",
    penalize_if_missing: float = 6.0,
) -> float:
    dists: List[int] = []

    for gene in disease_genes:
        if gene not in B:
            continue
        try:
            d = nx.shortest_path_length(B, source=drug, target=gene)
            dists.append(d)
        except nx.NetworkXNoPath:
            continue

    if not dists:
        return float(penalize_if_missing)

    if mode == "min":
        return float(min(dists))
    return float(sum(dists) / len(dists))


def rank_drugs_by_proximity(
    B: nx.Graph,
    disease_genes: Set[str],
    mode: str = "mean",
) -> pd.DataFrame:
    rows = []
    drugs = get_drug_nodes(B)

    for drug in drugs:
        dist = compute_drug_disease_distance(B, drug, disease_genes, mode=mode)
        rows.append({"drug": drug, "distance": dist})

    out = pd.DataFrame(rows).sort_values(["distance", "drug"], ascending=[True, True])
    return out


def draw_bipartite_png(B: nx.Graph, out_png: Path) -> None:
    node_sizes = []
    node_colors = []
    for n, d in B.nodes(data=True):
        if d.get("bipartite") == "drug":
            node_colors.append("blue")
            node_sizes.append(300 + 120 * B.degree(n))
        else:
            node_colors.append("red")
            node_sizes.append(220)

    pos = nx.spring_layout(B, seed=42)

    plt.figure(figsize=(10, 8))
    nx.draw(
        B,
        pos,
        node_color=node_colors,
        node_size=node_sizes,
        edge_color="gray",
        width=0.8,
        with_labels=False,
    )
    plt.title("Drug–Gene Bipartite Network")
    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    plt.close()


if __name__ == "__main__":
    ensure_exists(DRUG_GENE_CSV)
    ensure_exists(DISEASE_GENES_TXT)

    df = load_drug_gene_table(DRUG_GENE_CSV)
    drug2genes = build_drug2genes(df)
    if not drug2genes:
        raise ValueError("[ERROR] drug2genes gol — verifică CSV-ul.")
    B = build_bipartite_graph(drug2genes)

    disease_genes = load_disease_genes(DISEASE_GENES_TXT)
    if not disease_genes:
        raise ValueError("[ERROR] disease_genes gol — verifică fișierul .txt.")

    ranking = rank_drugs_by_proximity(B, disease_genes, mode="mean")

    ranking.to_csv(OUT_DRUG_PRIORITY, index=False)
    draw_bipartite_png(B, OUT_NETWORK_PNG)

    print("[INFO] Gata ")
    print(f"  - {OUT_DRUG_PRIORITY}")
    print(f"  - {OUT_NETWORK_PNG}")
