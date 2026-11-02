"""Lab 7: inspect and visualize the exact network exported by Lab 6.

Load adjacency.csv + modules.csv + parameters.json. No correlation or
threshold is recomputed here, so the continuation keeps the same graph.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import re

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]


def load_lab06_results(directory: Path) -> tuple[nx.Graph, dict[str, int], dict]:
    adjacency = pd.read_csv(directory / "adjacency.csv", index_col=0)
    adjacency.index = adjacency.index.map(str)
    adjacency.columns = adjacency.columns.map(str)
    if adjacency.empty or adjacency.index.has_duplicates or adjacency.columns.has_duplicates:
        raise ValueError("Matrice de adiacență goală sau identificatori de gene duplicați.")
    if set(adjacency.index) != set(adjacency.columns):
        raise ValueError("Liniile și coloanele trebuie să reprezinte aceleași gene.")
    adjacency = adjacency.loc[adjacency.index, adjacency.index].apply(pd.to_numeric, errors="raise")
    values = adjacency.to_numpy()
    if not np.isfinite(values).all() or (values < 0).any():
        raise ValueError("Ponderile trebuie să fie finite și nenegative.")
    if not np.allclose(values, values.T) or not np.allclose(np.diag(values), 0):
        raise ValueError("Rețeaua trebuie să fie simetrică, fără self-loops.")
    table = pd.read_csv(directory / "modules.csv", dtype={"Gene": str})
    if table.Gene.isna().any() or table.Gene.duplicated().any() or set(table.Gene) != set(adjacency.index):
        raise ValueError("Mapping-ul modulelor trebuie să acopere exact aceleași gene.")
    numeric_modules = pd.to_numeric(table.Module, errors="raise")
    if not np.isfinite(numeric_modules).all() or (numeric_modules % 1 != 0).any() or (numeric_modules < -1).any():
        raise ValueError("Identificatorii modulelor trebuie să fie întregi, cel puțin -1.")
    modules = dict(zip(table.Gene, numeric_modules.astype(int)))
    graph = nx.from_pandas_adjacency(adjacency)
    for gene, degree in graph.degree():
        if (degree == 0) != (modules[gene] == -1):
            raise ValueError("Genele izolate trebuie să aibă modul -1, iar cele conectate un modul detectat.")
    parameters = json.loads((directory / "parameters.json").read_text(encoding="utf-8"))
    actual_counts = {
        "retained_genes": graph.number_of_nodes(),
        "edges": graph.number_of_edges(),
        "isolates": len(list(nx.isolates(graph))),
        "modules": len({module for module in modules.values() if module >= 0}),
    }
    if not isinstance(parameters, dict) or any(parameters.get(key) != count for key, count in actual_counts.items()):
        raise ValueError("parameters.json nu corespunde rețelei și modulelor; încărcați exporturile aceleiași rulări.")
    return graph, modules, parameters


def compute_hubs(graph: nx.Graph, modules: dict[str, int], topk: int = 10) -> pd.DataFrame:
    """Rank connected genes by within-module weighted connectivity.

    Isolates are retained in the network, but never reported as hub candidates.
    """
    rows = []
    for gene in graph:
        if modules[gene] < 0 or graph.degree(gene) == 0:
            continue
        within = sum(attrs["weight"] for other, attrs in graph[gene].items() if modules[other] == modules[gene])
        rows.append({"Gene": gene, "Module": modules[gene], "Degree": graph.degree(gene),
                     "Strength": graph.degree(gene, weight="weight"), "IntramodularStrength": within})
    return pd.DataFrame(rows, columns=["Gene", "Module", "Degree", "Strength", "IntramodularStrength"]).sort_values(
        ["IntramodularStrength", "Degree", "Gene"], ascending=[False, False, True],
    ).head(topk)


def visualize(directory: Path, output_dir: Path, *, seed=42, topk=10) -> dict:
    graph, modules, parameters = load_lab06_results(directory)
    hubs = compute_hubs(graph, modules, topk)
    active = [gene for gene, degree in graph.degree() if degree > 0]
    isolates = [gene for gene, degree in graph.degree() if degree == 0]
    positions = nx.spring_layout(graph.subgraph(active), seed=seed, weight="weight") if active else {}
    # Place isolated genes in an outer ring; they remain visible and exported.
    for i, gene in enumerate(isolates):
        angle = 2 * np.pi * i / max(1, len(isolates))
        positions[gene] = np.array([1.5 * np.cos(angle), 1.5 * np.sin(angle)])
    module_ids = sorted({module for module in modules.values() if module >= 0})
    palette = plt.get_cmap("tab20")
    colors = {module: palette(i % 20) for i, module in enumerate(module_ids)}
    sizes = [40 + 6 * graph.degree(gene) for gene in graph]
    node_colors = [colors.get(modules[gene], "#999999") for gene in graph]
    fig, ax = plt.subplots(figsize=(10, 8), layout="constrained")
    nx.draw_networkx_edges(graph, positions, ax=ax, alpha=0.12, width=0.6, edge_color="#52636b")
    nx.draw_networkx_nodes(graph, positions, ax=ax, node_color=node_colors, node_size=sizes,
                           edgecolors="white", linewidths=0.5)
    labels = {gene: gene for gene in hubs.Gene if graph.degree(gene) > 0}
    nx.draw_networkx_labels(graph, positions, labels=labels, ax=ax, font_size=8)
    for module in module_ids:
        ax.scatter([], [], color=colors[module], label=f"Modul {module}", s=60)
    if isolates:
        ax.scatter([], [], color="#999999", label=f"Izolate (-1): {len(isolates)}", s=60)
    ax.legend(loc="upper left", frameon=False, fontsize=9)
    ax.set_title(f"Rețeaua din Lab 6: {graph.number_of_nodes()} gene, {graph.number_of_edges()} muchii")
    ax.set_axis_off()
    output_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_dir / "network.png", dpi=180)
    plt.close(fig)
    hubs.to_csv(output_dir / "hubs.csv", index=False)
    record = {"lab06_parameters": parameters, "layout": "spring; isolates on outer ring",
              "layout_seed": seed, "topk_labels": topk,
              "node_size": "40 + 6*degree", "hub_ranking": "intramodular weighted connectivity",
              "interpretation": "Layout distances are display coordinates, not expression distances or regulation."}
    (output_dir / "visualization_parameters.json").write_text(json.dumps(record, indent=2), encoding="utf-8")
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--handle", required=True)
    parser.add_argument("--lab06-dir", type=Path)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--topk", type=int, default=10)
    args = parser.parse_args()
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]*", args.handle):
        parser.error("--handle poate conține litere, cifre, _ și -.")
    if args.topk < 1:
        parser.error("--topk trebuie să fie pozitiv.")
    source = args.lab06_dir or ROOT / "data/work" / args.handle / "lab06"
    output = args.output_dir or ROOT / "data/work" / args.handle / "lab07"
    visualize(source, output, seed=args.seed, topk=args.topk)
    print(f"Rețea și hub-uri salvate în {output.resolve()}")


if __name__ == "__main__":
    main()
