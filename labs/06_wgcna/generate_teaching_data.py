"""Generate the offline, explicitly synthetic Lab 6 teaching dataset.

Values are continuous, prepared log2-like expression, not RNA-seq counts.
Seed 20261101; 48 independent synthetic samples; 120 invented gene IDs.
The planted groups are generation metadata, not inferred module assignments.
"""
from pathlib import Path
import json

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = ROOT / "data/sample/lab06"


def generate_data(seed=20261101):
    rng = np.random.default_rng(seed)
    samples = [f"S{i:02d}" for i in range(1, 49)]
    genes = [f"G{i:03d}" for i in range(1, 121)]
    treated = np.repeat([0, 1], 24)
    response = 0.9 * (2 * treated - 1) + rng.normal(0, 0.7, 48)
    independent = rng.normal(0, 1, 48)
    response = (response - response.mean()) / response.std(ddof=1)
    independent = (independent - independent.mean()) / independent.std(ddof=1)
    latent = [response, independent, -response]
    values = np.empty((120, 48))
    planted = []
    for i in range(120):
        baseline = rng.uniform(6.5, 9.5)
        residual_sd = rng.uniform(0.45, 0.65)
        if i < 90:
            group = i // 30
            loading = rng.uniform(0.9, 1.1)
            values[i] = baseline + loading * latent[group] + rng.normal(0, residual_sd, 48)
            planted.append(f"P{group + 1}")
        else:
            values[i] = baseline + rng.normal(0, 0.8, 48)
            planted.append("noise")
    expr = pd.DataFrame(values, index=genes, columns=samples)
    expr.index.name = "Gene"
    meta = pd.DataFrame({"Sample": samples, "Condition": np.where(treated, "tratament", "control")})
    truth = pd.DataFrame({"Gene": genes, "PlantedGroup": planted})
    return expr, meta, truth


def main():
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    expr, meta, truth = generate_data()
    expr.to_csv(DATA_DIR / "expression_matrix.csv", float_format="%.8f")
    meta.to_csv(DATA_DIR / "sample_metadata.csv", index=False)
    truth.to_csv(DATA_DIR / "generation_groups.csv", index=False)
    record = {
        "synthetic": True, "seed": 20261101, "genes": 120, "samples": 48,
        "numpy_version": np.__version__, "pandas_version": pd.__version__,
        "unit": "prepared continuous log2-like expression; not raw counts",
        "groups": "P1: G001-G030, P2: G031-G060, P3: G061-G090, noise: G091-G120",
        "latent": "P1 follows standardized response, P2 an independent latent variable, P3 negative response",
        "residual_sd": "uniform 0.45 to 0.65 per planted gene",
        "loading": "uniform 0.9 to 1.1", "baseline": "uniform 6.5 to 9.5",
        "response": "0.9*(2*treatment-1)+N(0,0.7); then standardized",
        "noise_genes_sd": 0.8,
        "warning": "Invented gene IDs: no GO/KEGG enrichment or biological target claims. Groups are not inferred modules.",
    }
    (DATA_DIR / "provenance.json").write_text(json.dumps(record, indent=2), encoding="utf-8")
    print(DATA_DIR)


if __name__ == "__main__":
    main()
