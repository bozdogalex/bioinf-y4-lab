"""Optional Lab 7: module eigengenes and descriptive treatment associations.

Reuse Lab 6's modules and prepared expression. Standardize genes, calculate
PC1 by SVD, orient toward the mean gene profile, and standardize scores.
The bundled control/tratament simulation has no genuine GO/KEGG annotations.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from ex01_network_viz import ROOT, load_lab06_results


def compute_profiles(expression: pd.DataFrame, modules: dict[str, int]) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return sample-by-module standardized scores and PC1 variance summaries."""
    if (expression.index.has_duplicates or expression.columns.has_duplicates
            or set(expression.index) != set(modules) or expression.shape[1] < 3):
        raise ValueError("Expresia trebuie să conțină aceleași gene și cel puțin trei probe distincte.")
    expression = expression.apply(pd.to_numeric, errors="raise")
    if not np.isfinite(expression.to_numpy()).all():
        raise ValueError("Expresia trebuie să conțină numai valori finite.")
    module_ids = sorted({value for value in modules.values() if value >= 0})
    if not module_ids:
        raise ValueError("Rețeaua nu conține module conectate pentru calculul profilurilor.")
    profiles, rows = {}, []
    for module in module_ids:
        genes = sorted(gene for gene, label in modules.items() if label == module)
        block = expression.loc[genes]
        deviation = block.std(axis=1, ddof=1)
        if (deviation == 0).any():
            raise ValueError("Genele constante trebuie filtrate înainte de analiza modulelor.")
        standardized = block.sub(block.mean(axis=1), axis=0).div(deviation, axis=0).T
        left, singular, _ = np.linalg.svd(standardized.to_numpy(), full_matrices=False)
        score = left[:, 0] * singular[0]
        if np.dot(score, standardized.mean(axis=1).to_numpy()) < 0:
            score = -score
        profiles[f"Module_{module}"] = score / score.std(ddof=1)
        rows.append({"Module": module, "Genes": len(genes),
                     "PC1VarianceExplained": float(singular[0] ** 2 / np.sum(singular ** 2))})
    return pd.DataFrame(profiles, index=expression.columns), pd.DataFrame(rows)


def align_conditions(metadata: pd.DataFrame, samples: pd.Index) -> pd.Series:
    """Align the course's two condition labels to expression sample identifiers."""
    if ("Condition" not in metadata or metadata.index.has_duplicates
            or set(metadata.index) != set(samples)):
        raise ValueError("Metadatele trebuie să identifice exact aceleași probe și coloana Condition.")
    condition = metadata.loc[samples, "Condition"]
    if set(condition) != {"control", "tratament"}:
        raise ValueError("Acest exemplu așteaptă condițiile control și tratament, ambele prezente.")
    return condition


def analyze(directory: Path, output: Path) -> pd.DataFrame:
    _, modules, parameters = load_lab06_results(directory)
    expression_path = directory / "processed_expression.csv"
    metadata_path = directory / "sample_metadata.csv"
    expression = pd.read_csv(expression_path, index_col=0)
    expression.index = expression.index.map(str)
    expression.columns = expression.columns.map(str)
    metadata = pd.read_csv(metadata_path, index_col=0)
    metadata.index = metadata.index.map(str)
    condition = align_conditions(metadata, expression.columns)
    if expression.shape[1] != parameters.get("samples"):
        raise ValueError("Numărul probelor diferă de cel salvat în parameters.json.")
    profiles, summary = compute_profiles(expression, modules)
    treatment = (condition == "tratament").to_numpy(dtype=int)
    summary["TreatmentCorrelation"] = [
        float(np.corrcoef(profiles[f"Module_{module}"], treatment)[0, 1])
        for module in summary.Module
    ]
    summary["ControlSamples"] = int((condition == "control").sum())
    summary["TreatmentSamples"] = int((condition == "tratament").sum())
    fig, ax = plt.subplots(figsize=(9, 4.5), layout="constrained")
    image = ax.imshow(profiles.T, aspect="auto", cmap="RdBu_r", vmin=-2.5, vmax=2.5)
    ax.set_yticks(range(len(profiles.columns)), profiles.columns)
    ax.set_xlabel("Probe în ordinea matricei de expresie")
    ax.set_title("Eigengene-uri standardizate; condițiile sunt în sample_metadata.csv")
    fig.colorbar(image, ax=ax, label="Scor PC1 standardizat")
    output.mkdir(parents=True, exist_ok=True)
    profiles.to_csv(output / "module_profiles.csv", index_label="Sample")
    summary.to_csv(output / "module_summary.csv", index=False)
    fig.savefig(output / "module_profiles.png", dpi=180)
    plt.close(fig)
    record = {
        "lab06_parameters": parameters,
        "expression_sha256": hashlib.sha256(expression_path.read_bytes()).hexdigest(),
        "metadata_sha256": hashlib.sha256(metadata_path.read_bytes()).hexdigest(),
        "method": "per-gene z-score (ddof=1); SVD PC1; sample scores scaled to SD=1",
        "orientation": "PC1 aligned with mean standardized gene expression when nonzero",
        "treatment_encoding": {"control": 0, "tratament": 1},
        "interpretation": "Descriptive association in the discovery samples; no confirmatory p-values",
    }
    (output / "profile_parameters.json").write_text(json.dumps(record, indent=2), encoding="utf-8")
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--handle", required=True)
    parser.add_argument("--lab06-dir", type=Path)
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args()
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]*", args.handle):
        parser.error("--handle poate conține litere, cifre, _ și -.")
    source = args.lab06_dir or ROOT / "data/work" / args.handle / "lab06"
    output = args.output_dir or ROOT / "data/work" / args.handle / "lab07/profiles"
    try:
        summary = analyze(source, output)
    except (ValueError, FileNotFoundError) as error:
        parser.exit(1, f"{error}\nFolosiți exportul complet al variantei de bază din Lab 6.\n")
    print(summary.to_string(index=False))
    print(f"Profiluri și rezumat descriptiv salvate în {output.resolve()}")


if __name__ == "__main__":
    main()
