"""Student exercise: compare six network definitions on one correlation matrix.

The complete baseline remains ex01_gce_networks.py. Implement compare_networks
below; the driver prepares expression and calculates correlation only once.
No baseline adjacency/module file is written or changed by this exercise.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re

import networkx as nx
import pandas as pd

from ex01_gce_networks import (
    ROOT, DEFAULT_INPUT, read_expression_matrix, log_and_filter,
    correlation_matrix, adjacency_from_correlation,
    graph_from_adjacency, detect_modules_louvain,
)


def compare_networks(corr: pd.DataFrame, *, thresholds=(0.5, 0.7, 0.9),
                     networks=("positive", "unsigned"), seed=42,
                     resolution=1.0) -> pd.DataFrame:
    """TODO: return one row per (network, threshold), using supplied corr.

    Required columns: Network, Threshold, Edges, Density, Isolates, Modules,
    LargestModule. Use weighted retained edges, the supplied seed/resolution
    and the same genes in every configuration. Exclude Module=-1 from module
    counts and sizes; set LargestModule=0 if no module is detected.

    Reuse adjacency_from_correlation, graph_from_adjacency and
    detect_modules_louvain. NetworkX graph methods supply counts/density;
    pandas Series.value_counts or collections.Counter can summarize sizes.
    Do not recalculate correlation, overwrite the baseline, or use the
    simulation's generation_groups.csv to assign modules.
    """
    raise NotImplementedError("Completați funcția compare_networks conform enunțului din README.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--handle", required=True)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--method", choices=["pearson", "spearman"], default="spearman")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--resolution", type=float, default=1.0)
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args()
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]*", args.handle):
        parser.error("--handle poate conține litere, cifre, _ și -.")
    expression = read_expression_matrix(args.input)
    prepared = log_and_filter(expression, variance_threshold=0.05, input_scale="prepared")
    corr = correlation_matrix(prepared, method=args.method)
    try:
        table = compare_networks(corr, seed=args.seed, resolution=args.resolution)
    except NotImplementedError as error:
        parser.exit(1, f"{error}\nBaseline-ul executabil este ex01_gce_networks.py.\n")
    required = {"Network", "Threshold", "Edges", "Density", "Isolates", "Modules", "LargestModule"}
    if not isinstance(table, pd.DataFrame) or not required.issubset(table.columns):
        raise ValueError("Tabelul trebuie să conțină coloanele din contractul compare_networks.")
    expected = {(network, threshold) for network in ("positive", "unsigned") for threshold in (0.5, 0.7, 0.9)}
    if len(table) != 6 or set(zip(table.Network, table.Threshold)) != expected:
        raise ValueError("Tabelul trebuie să conțină exact cele șase configurații cerute.")
    table = table.copy()
    table["Correlation"] = args.method
    table["Weighted"] = True
    table["Seed"] = args.seed
    table["Resolution"] = args.resolution
    table["Genes"] = len(corr)
    output = args.output_dir or ROOT / "data/work" / args.handle / "lab06/comparisons"
    output.mkdir(parents=True, exist_ok=True)
    table.sort_values(["Network", "Threshold"]).to_csv(output / "comparison.csv", index=False)
    params = {"input_sha256": hashlib.sha256(args.input.read_bytes()).hexdigest(),
              "input_scale": "prepared", "variance_threshold": 0.05,
              "samples": prepared.shape[1], "genes": len(corr),
              "correlation": args.method, "seed": args.seed, "resolution": args.resolution,
              "networkx_version": nx.__version__, "weighted": True}
    (output / "comparison_parameters.json").write_text(json.dumps(params, indent=2), encoding="utf-8")
    print(table.to_string(index=False))
    print(f"Comparație salvată în {output}; baseline-ul pentru Lab 7 rămâne separat.")


if __name__ == "__main__":
    main()
