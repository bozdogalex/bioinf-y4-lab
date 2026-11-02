"""Lab 6: runnable baseline for comparing thresholded correlation networks.

This is thresholding + Louvain, not WGCNA/TOM. Default data are explicitly
synthetic prepared expression: no normalization or log transform is reapplied.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import re
import warnings

import networkx as nx
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_INPUT = ROOT / "data/sample/lab06/expression_matrix.csv"
DEFAULT_METADATA = ROOT / "data/sample/lab06/sample_metadata.csv"


def read_expression_matrix(path: Path) -> pd.DataFrame:
    """Read genes on rows, independent samples on columns; reject bad input."""
    with path.open(encoding="utf-8-sig", newline="") as source:
        header = next(csv.reader(source), [])
    sample_names = header[1:]
    if len(sample_names) < 3 or any(not name.strip() for name in sample_names) or len(set(sample_names)) != len(sample_names):
        raise ValueError("Antetul trebuie să conțină cel puțin trei probe cu nume unice și nenule.")
    df = pd.read_csv(path, index_col=0)
    if df.empty or df.shape[1] < 3:
        raise ValueError("Matricea trebuie să conțină gene și cel puțin trei probe.")
    if df.index.has_duplicates or df.columns.has_duplicates or df.index.isna().any():
        raise ValueError("Identificatorii genelor și probelor trebuie să fie unici și nenuli.")
    df.index = df.index.map(str)
    df.columns = df.columns.map(str)
    if any(not gene.strip() for gene in df.index):
        raise ValueError("Identificatorii genelor trebuie să fie nenuli, fără nume formate doar din spații.")
    df = df.apply(pd.to_numeric, errors="raise")
    if not np.isfinite(df.to_numpy()).all():
        raise ValueError("Matricea conține valori lipsă sau nefinite.")
    return df


def log_and_filter(df: pd.DataFrame, variance_threshold: float = 0.05,
                   input_scale: str = "prepared") -> pd.DataFrame:
    """Filter constant/low-variance genes; optionally log normalized counts.

    'normalized_counts' applies log2(x+1), which is NOT sample normalization.
    'prepared' leaves values on their supplied scale.
    """
    if variance_threshold < 0 or not np.isfinite(variance_threshold):
        raise ValueError("Pragul de varianță trebuie să fie finit și nenegativ.")
    if input_scale == "normalized_counts":
        if (df.to_numpy() < 0).any():
            raise ValueError("Valorile pentru log2(x+1) trebuie să fie nenegative.")
        df = np.log2(df + 1)
    elif input_scale != "prepared":
        raise ValueError("input_scale trebuie să fie prepared sau normalized_counts.")
    variance = df.var(axis=1, ddof=1)
    filtered = df.loc[(variance > 0) & (variance >= variance_threshold)].copy()
    if len(filtered) < 2:
        raise ValueError("După filtrare trebuie să rămână cel puțin două gene variabile.")
    return filtered


def correlation_matrix(df: pd.DataFrame, method: str = "spearman") -> pd.DataFrame:
    """Correlate gene rows across samples, retaining the original sign."""
    if method not in {"pearson", "spearman"}:
        raise ValueError("Metoda trebuie să fie pearson sau spearman.")
    corr = df.T.corr(method=method)
    if not np.isfinite(corr.to_numpy()).all():
        raise ValueError("Corelații nefinite: verificați genele constante și valorile lipsă.")
    return corr.clip(-1, 1)


def adjacency_from_correlation(corr: pd.DataFrame, threshold: float,
                               weighted: bool = False,
                               use_abs: bool = False) -> pd.DataFrame:
    """Positive-only r>=threshold or unsigned |r|>=threshold; no self loops.

    Retained weights are r (positive-only) or |r| (unsigned). Positive-only
    thresholding is distinct from WGCNA's signed power transformation.
    """
    if not 0 < threshold <= 1:
        raise ValueError("Pragul trebuie să fie în intervalul (0, 1].")
    strengths = corr.abs() if use_abs else corr
    adjacency = strengths.where(strengths >= threshold, 0.0) if weighted else (strengths >= threshold).astype(float)
    np.fill_diagonal(adjacency.values, 0.0)
    return adjacency


def graph_from_adjacency(adjacency: pd.DataFrame) -> nx.Graph:
    """Preserve isolated genes so the Lab 7 handoff retains all identifiers."""
    return nx.from_pandas_adjacency(adjacency)


def detect_modules_louvain(graph: nx.Graph, seed: int = 42,
                           resolution: float = 1.0) -> dict[str, int]:
    """Louvain on non-isolates; -1 marks genes with no retained edge."""
    if resolution <= 0 or not np.isfinite(resolution):
        raise ValueError("Rezoluția trebuie să fie pozitivă și finită.")
    mapping = {gene: -1 for gene in graph}
    active = graph.subgraph([gene for gene, degree in graph.degree() if degree > 0])
    if active.number_of_edges() == 0:
        return mapping
    communities = nx.community.louvain_communities(active, weight="weight", resolution=resolution, seed=seed)
    for module, members in enumerate(sorted(communities, key=lambda group: min(group))):
        for gene in members:
            mapping[gene] = module
    return mapping


def run_analysis(input_csv: Path, output_dir: Path, *, method="spearman",
                 threshold=0.7, use_abs=False, weighted=True, seed=42,
                 resolution=1.0, variance_threshold=0.05,
                 input_scale="prepared", metadata_csv: Path | None = None,
                 overwrite=False) -> dict:
    expr = read_expression_matrix(input_csv)
    if expr.shape[1] < 15:
        warnings.warn("Mai puțin de 15 probe: folosiți rezultatul doar ca demonstrație.")
    processed = log_and_filter(expr, variance_threshold, input_scale)
    corr = correlation_matrix(processed, method)
    adjacency = adjacency_from_correlation(corr, threshold, weighted, use_abs)
    graph = graph_from_adjacency(adjacency)
    mapping = detect_modules_louvain(graph, seed, resolution)
    modules = pd.DataFrame({"Gene": list(mapping), "Module": list(mapping.values())})
    modules["Isolate"] = modules["Module"].eq(-1)
    edges = pd.DataFrame([
        {"Source": source, "Target": target, "Weight": attrs["weight"], "Correlation": corr.loc[source, target]}
        for source, target, attrs in graph.edges(data=True)
    ], columns=["Source", "Target", "Weight", "Correlation"])
    metadata = None
    if metadata_csv is not None:
        metadata = pd.read_csv(metadata_csv, index_col="Sample")
        metadata.index = metadata.index.map(str)
        if metadata.index.has_duplicates or set(metadata.index) != set(processed.columns):
            raise ValueError("Metadatele trebuie să conțină exact câte o linie pentru fiecare probă.")
        metadata = metadata.loc[processed.columns]
    parameters = {
        "pipeline": "thresholded correlation network + Louvain; not WGCNA",
        "input_file": input_csv.name, "input_sha256": hashlib.sha256(input_csv.read_bytes()).hexdigest(),
        "bundled_synthetic_input": input_csv.resolve() == DEFAULT_INPUT.resolve(),
        "sample_metadata_included": metadata is not None,
        "metadata_sha256": hashlib.sha256(metadata_csv.read_bytes()).hexdigest() if metadata_csv else None,
        "input_scale": input_scale, "correlation": method,
        "network_type": "unsigned_abs" if use_abs else "positive_only",
        "threshold": threshold, "weighted": weighted, "variance_threshold": variance_threshold,
        "seed": seed, "resolution": resolution, "samples": processed.shape[1],
        "input_genes": expr.shape[0], "retained_genes": len(processed),
        "edges": graph.number_of_edges(), "isolates": len(list(nx.isolates(graph))),
        "modules": len({m for m in mapping.values() if m >= 0}),
        "module_minus_one": "isolated gene, not a detected module",
        "versions": {"numpy": np.__version__, "pandas": pd.__version__, "networkx": nx.__version__},
    }
    previous_path = output_dir / "parameters.json"
    if previous_path.exists() and not overwrite:
        try:
            previous = json.loads(previous_path.read_text(encoding="utf-8"))
        except (ValueError, OSError) as error:
            raise ValueError("parameters.json nu poate fi citit; alegeți alt folder sau --overwrite.") from error
        config_keys = (
            "input_sha256", "metadata_sha256", "input_scale", "correlation",
            "network_type", "threshold", "weighted", "variance_threshold",
            "seed", "resolution", "versions",
        )
        if not isinstance(previous, dict) or any(previous.get(key) != parameters[key] for key in config_keys):
            raise ValueError("Folderul conține altă analiză; folosiți --output-dir separat sau --overwrite explicit.")
    output_dir.mkdir(parents=True, exist_ok=True)
    processed.to_csv(output_dir / "processed_expression.csv", index_label="Gene")
    corr.to_csv(output_dir / "correlation.csv", index_label="Gene")
    adjacency.to_csv(output_dir / "adjacency.csv", index_label="Gene")
    modules.to_csv(output_dir / "modules.csv", index=False)
    edges.to_csv(output_dir / "edges.csv", index=False)
    if metadata is not None:
        metadata.to_csv(output_dir / "sample_metadata.csv", index_label="Sample")
    elif (output_dir / "sample_metadata.csv").exists():
        # A successful rerun without metadata must not retain a stale cohort file.
        (output_dir / "sample_metadata.csv").unlink()
    (output_dir / "parameters.json").write_text(json.dumps(parameters, indent=2), encoding="utf-8")
    return parameters


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--handle", required=True, help="Identificator pentru folderul de rezultate")
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--metadata", type=Path)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--overwrite", action="store_true", help="Înlocuiește explicit o analiză existentă cu alți parametri")
    parser.add_argument("--method", choices=["pearson", "spearman"], default="spearman")
    parser.add_argument("--threshold", type=float, default=0.7)
    parser.add_argument("--network", choices=["positive", "unsigned"], default="positive")
    parser.add_argument("--binary", action="store_true", help="Ponderi 0/1 în locul intensităților păstrate")
    parser.add_argument("--variance-threshold", type=float, default=0.05)
    parser.add_argument("--input-scale", choices=["prepared", "normalized_counts"], default="prepared")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--resolution", type=float, default=1.0)
    args = parser.parse_args()
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]*", args.handle):
        parser.error("--handle poate conține litere, cifre, _ și -.")
    out = args.output_dir or ROOT / "data/work" / args.handle / "lab06"
    metadata = args.metadata
    if metadata is None and args.input.resolve() == DEFAULT_INPUT.resolve():
        metadata = DEFAULT_METADATA
    report = run_analysis(
        args.input, out, method=args.method, threshold=args.threshold,
        use_abs=args.network == "unsigned", weighted=not args.binary,
        seed=args.seed, resolution=args.resolution, variance_threshold=args.variance_threshold,
        input_scale=args.input_scale, metadata_csv=metadata, overwrite=args.overwrite,
    )
    print(f"{report['retained_genes']} gene, {report['edges']} muchii, "
          f"{report['modules']} module, {report['isolates']} gene izolate.")
    print(f"Rezultate pentru Lab 7: {out.resolve()}")


if __name__ == "__main__":
    main()
