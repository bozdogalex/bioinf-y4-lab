"""Schematic: two triangles joined by one bridge, not inferred expression data.

Community labels are assigned manually. Larger nodes have the largest global
degree; this differs from the within-module weighted ranking in ex01.
"""
from pathlib import Path

import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]


def main():
    graph = nx.Graph()
    graph.add_edges_from([
        ("GeneA", "GeneB"), ("GeneA", "GeneC"), ("GeneB", "GeneC"),
        ("GeneC", "GeneD"),
        ("GeneD", "GeneE"), ("GeneE", "GeneF"), ("GeneD", "GeneF"),
    ])
    manual_modules = {
        "GeneA": 1, "GeneB": 1, "GeneC": 1,
        "GeneD": 2, "GeneE": 2, "GeneF": 2,
    }
    palette = plt.get_cmap("tab10")
    colors = [palette(manual_modules[gene] - 1) for gene in graph]
    degree = dict(graph.degree())
    highlighted = sorted(degree.items(), key=lambda item: (-item[1], item[0]))[:2]
    candidates = {gene for gene, _ in highlighted}
    positions = nx.spring_layout(graph, seed=42)
    fig, ax = plt.subplots(figsize=(6, 5), layout="constrained")
    nx.draw_networkx_edges(graph, positions, alpha=0.4, ax=ax)
    nx.draw_networkx_nodes(graph, positions, ax=ax, node_color=colors,
                           node_size=[500 if gene in candidates else 250 for gene in graph])
    nx.draw_networkx_labels(graph, positions, font_size=10, ax=ax)
    ax.set_title("Două triunghiuri și o punte — module atribuite manual")
    ax.set_axis_off()
    output = ROOT / "data/work/demo/lab07"
    output.mkdir(parents=True, exist_ok=True)
    fig.savefig(output / "demo_network.png", dpi=180)
    plt.close(fig)
    pd.DataFrame(highlighted, columns=["Gene", "Degree"]).to_csv(output / "demo_hubs.csv", index=False)
    print(f"Schemă cu {graph.number_of_nodes()} noduri și {graph.number_of_edges()} muchii.")
    print(f"Noduri evidențiate după gradul global: {highlighted}.")
    print("Modulele sunt atribuite manual; datele nu provin din expresie genică.")
    print(f"Figura și exemplul de grad sunt salvate local în {output}.")


if __name__ == "__main__":
    main()
