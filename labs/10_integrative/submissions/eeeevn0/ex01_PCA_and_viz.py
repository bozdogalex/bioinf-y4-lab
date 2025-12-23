"""
Exercise 10 — PCA Single-Omics vs Joint

TODO:
- încărcați SNP și Expression
- normalizați fiecare strat (z-score)
- rulați PCA pe:
    1) strat SNP
    2) strat Expression
    3) strat Joint (concat)
- generați 3 figuri PNG
- comparați vizual distribuția probelor
"""

from pathlib import Path
import pandas as pd
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt

HANDLE = "eeeevn0"

SNP_CSV = Path(f"data/work/{HANDLE}/lab10/snp_matrix_{HANDLE}.csv")
EXP_CSV = Path(f"data/work/{HANDLE}/lab10/expression_matrix_{HANDLE}.csv")

OUT_DIR = Path(f"labs/10_integrative/submissions/{HANDLE}")
OUT_DIR.mkdir(parents=True, exist_ok=True)

# incarcare date
df_snp = pd.read_csv(SNP_CSV, index_col=0)
df_exp = pd.read_csv(EXP_CSV, index_col=0)

# aliniere 
common_samples = df_snp.index.intersection(df_exp.index)
df_snp = df_snp.loc[common_samples]
df_exp = df_exp.loc[common_samples]

# normalizare z-score
df_snp_norm = (df_snp - df_snp.mean()) / df_snp.std()
df_exp_norm = (df_exp - df_exp.mean()) / df_exp.std()

# eliminare coloane cu NaN (daca exista)
df_snp_norm = df_snp_norm.dropna(axis=1)
df_exp_norm = df_exp_norm.dropna(axis=1)

# rulare PCA si salvare plot
def run_pca(df, title, out_png):
    pca = PCA(n_components=2)
    proj = pca.fit_transform(df)

    plt.figure(figsize=(5, 4))
    plt.scatter(proj[:, 0], proj[:, 1])
    plt.title(title)
    plt.xlabel("PC1")
    plt.ylabel("PC2")
    plt.tight_layout()
    plt.savefig(out_png)
    plt.close()

run_pca(
    df_snp_norm,
    "PCA — SNP only",
    OUT_DIR / f"pca_snp_{HANDLE}.png"
)

run_pca(
    df_exp_norm,
    "PCA — Expression only",
    OUT_DIR / f"pca_expression_{HANDLE}.png"
)
# concatenare date pentru PCA multi-omics
df_joint = pd.concat([df_snp_norm, df_exp_norm], axis=1)
df_joint.to_csv(OUT_DIR / f"multiomics_concat_{HANDLE}.csv")

run_pca(
    df_joint,
    "PCA — Joint Multi-Omics",
    OUT_DIR / f"pca_joint_{HANDLE}.png"
)
print("Exercise completed. PCA plots saved.")
