"""
Exercise 10.2 — Identify top SNP–Gene correlations
Rezolvă Task 3 — Cross-Omics Correlation
"""

from pathlib import Path
import pandas as pd
import numpy as np


HANDLE = "Razvann19"

JOINT_CSV = Path(f"labs/10_integrative/submissions/Razvann19/multiomics_concat_Razvann19.csv")
OUT_CSV = Path(f"labs/10_integrative/submissions/{HANDLE}/snp_gene_pairs_{HANDLE}.csv")

def ensure_exists(p: Path) -> None:
    if not p.exists() or not p.is_file():
        raise FileNotFoundError(f"[ERROR] Nu găsesc fișierul: {p}")

if __name__ == "__main__":
    ensure_exists(JOINT_CSV)

    df = pd.read_csv(JOINT_CSV)

    if "SampleID" in df.columns:
        df = df.set_index("SampleID")

    df = df.apply(pd.to_numeric, errors="coerce")

    snp_cols = [c for c in df.columns if c.upper().startswith("SNP")]
    gene_cols = [c for c in df.columns if (c.upper().startswith("GENE") or c.upper().startswith("GENE_") or c.upper().startswith("G")) and not c.upper().startswith("SNP")]

    if not gene_cols:
        gene_cols = [c for c in df.columns if c not in snp_cols]

    if not snp_cols or not gene_cols:
        raise ValueError(
            f"[ERROR] Nu pot separa SNP-uri de gene.\n"
            f"SNP cols: {snp_cols}\nGene cols: {gene_cols}\n"
            f"Verifică numele coloanelor."
        )

    rows = []
    for snp in snp_cols:
        for gene in gene_cols:
            tmp = df[[snp, gene]].dropna()
            if tmp.shape[0] < 2:
                continue
            x = tmp[snp].values
            y = tmp[gene].values

            if np.std(x) == 0 or np.std(y) == 0:
                continue

            r = np.corrcoef(x, y)[0, 1]
            if np.isfinite(r) and abs(r) > 0.5:
                rows.append({"snp": snp, "gene": gene, "correlation": float(r), "n_samples": int(tmp.shape[0])})

    out = pd.DataFrame(rows).sort_values("correlation", key=lambda s: s.abs(), ascending=False)
    out.to_csv(OUT_CSV, index=False)

    print("[INFO] Task 3 gata ")
    print(f"  - SNP features: {len(snp_cols)} | gene features: {len(gene_cols)}")
    print(f"  - pairs kept (|r|>0.5): {len(out)}")
    print(f"  - {OUT_CSV}")