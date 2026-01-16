"""
Exercise 10 — PCA Single-Omics vs Joint

Rezolvă:
- Task 2: PCA pe SNPs, Expression, Joint + 3 figuri PNG
Include preprocesarea necesară (din Task 1):
- load SNP/Expression
- alinere probe comune
- z-score pe fiecare strat
- concat pe features
"""

from __future__ import annotations
from pathlib import Path

import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt


HANDLE = "Razvann19"

SNP_CSV = Path(f"labs/10_integrative/submissions/Razvann19/snp_matrix_demo.csv")
EXP_CSV = Path(f"labs/10_integrative/submissions/Razvann19/expression_matrix_demo.csv")

META_CSV = Path(f"data/work/{HANDLE}/lab10/sample_metadata_{HANDLE}.csv")

OUT_DIR = Path(f"labs/10_integrative/submissions/{HANDLE}")
OUT_DIR.mkdir(parents=True, exist_ok=True)

OUT_PCA_SNP = OUT_DIR / f"pca_snp_{HANDLE}.png"
OUT_PCA_EXPR = OUT_DIR / f"pca_expr_{HANDLE}.png"
OUT_PCA_JOINT = OUT_DIR / f"pca_joint_{HANDLE}.png"

OUT_JOINT = OUT_DIR / f"multiomics_concat_{HANDLE}.csv"


def ensure_exists(p: Path) -> None:
    if not p.exists() or not p.is_file():
        raise FileNotFoundError(f"[ERROR] Nu găsesc fișierul: {p}")


def read_matrix(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, index_col=0)
    df.index = df.index.astype(str).str.strip()
    df.columns = df.columns.astype(str).str.strip()
    return df


def to_samples_by_features(df: pd.DataFrame, sample_names: pd.Index) -> pd.DataFrame:
    if len(sample_names.intersection(df.columns)) >= len(sample_names) * 0.8:
        out = df[sample_names].T 
        out.index.name = "SampleID"
        return out

    if len(sample_names.intersection(df.index)) >= len(sample_names) * 0.8:
        out = df.loc[sample_names]  
        out.index.name = "SampleID"
        return out

    raise ValueError(
        "[ERROR] Nu pot detecta orientarea (samples pe rows/cols). "
        "Verifică dacă sample-urile au aceeași convenție între fișiere."
    )


def zscore_features(X: pd.DataFrame) -> pd.DataFrame:
    mu = X.mean(axis=0)
    sd = X.std(axis=0, ddof=0).replace(0, np.nan)
    Z = (X - mu) / sd
    return Z


def load_metadata(meta_path: Path, samples: pd.Index) -> pd.Series | None:

    if not meta_path.exists():
        return None
    meta = pd.read_csv(meta_path)
    if "SampleID" not in meta.columns:
        return None
    label_col = "Condition" if "Condition" in meta.columns else ("Subtype" if "Subtype" in meta.columns else None)
    if label_col is None:
        return None

    meta["SampleID"] = meta["SampleID"].astype(str).str.strip()
    meta[label_col] = meta[label_col].astype(str).str.strip()
    meta = meta.set_index("SampleID")[label_col]

    meta = meta.reindex(samples)
    return meta


def pca_2d(X: pd.DataFrame, n_components: int = 2) -> pd.DataFrame:
    pca = PCA(n_components=n_components, random_state=42)
    coords = pca.fit_transform(X.values)
    out = pd.DataFrame(coords, index=X.index, columns=[f"PC{i+1}" for i in range(n_components)])
    out.attrs["explained"] = pca.explained_variance_ratio_
    return out


def plot_pca(coords: pd.DataFrame, labels: pd.Series | None, out_png: Path, title: str) -> None:
    plt.figure(figsize=(7, 5))

    explained = coords.attrs.get("explained", [np.nan, np.nan])
    xlab = f"PC1 ({explained[0]*100:.1f}%)" if not np.isnan(explained[0]) else "PC1"
    ylab = f"PC2 ({explained[1]*100:.1f}%)" if not np.isnan(explained[1]) else "PC2"

    if labels is None:
        plt.scatter(coords["PC1"], coords["PC2"])
    else:
        for cls in sorted(labels.dropna().unique()):
            idx = labels[labels == cls].index
            plt.scatter(coords.loc[idx, "PC1"], coords.loc[idx, "PC2"], label=str(cls))
        plt.legend(title=labels.name or "Label", loc="best", frameon=True)

    plt.xlabel(xlab)
    plt.ylabel(ylab)
    plt.title(title)
    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    plt.close()


if __name__ == "__main__":
    ensure_exists(SNP_CSV)
    ensure_exists(EXP_CSV)

    df_snp_raw = read_matrix(SNP_CSV)
    df_exp_raw = read_matrix(EXP_CSV)

    samples_snp = pd.Index(df_snp_raw.columns).union(df_snp_raw.index)
    samples_exp = pd.Index(df_exp_raw.columns).union(df_exp_raw.index)

    common = pd.Index(sorted(set(samples_snp).intersection(set(samples_exp))))
    if len(common) < 2:
        raise ValueError(f"[ERROR] Prea puține probe comune: {len(common)}. Verifică fișierele.")

    X_snp = to_samples_by_features(df_snp_raw, common)
    X_exp = to_samples_by_features(df_exp_raw, common)

    X_snp_z = zscore_features(X_snp)
    X_exp_z = zscore_features(X_exp)


    X_joint = pd.concat([X_snp_z, X_exp_z], axis=1)


    X_joint.to_csv(OUT_JOINT)

    labels = load_metadata(META_CSV, X_joint.index)
    if labels is None:
        print("[WARN] Nu am metadata (sample_metadata_<handle>.csv). PCA va fi necolorat (o singură clasă).")
    else:
        print(f"[INFO] Metadata găsit: {META_CSV} (coloană: {labels.name})")

    coords_snp = pca_2d(X_snp_z.fillna(0))
    coords_exp = pca_2d(X_exp_z.fillna(0))
    coords_joint = pca_2d(X_joint.fillna(0))

    plot_pca(coords_snp, labels, OUT_PCA_SNP, "PCA — SNPs (z-score)")
    plot_pca(coords_exp, labels, OUT_PCA_EXPR, "PCA — Expression (z-score)")
    plot_pca(coords_joint, labels, OUT_PCA_JOINT, "PCA — Joint (SNPs + Expression)")

    print("[INFO] Gata ")
    print(f"  - {OUT_JOINT}")
    print(f"  - {OUT_PCA_SNP}")
    print(f"  - {OUT_PCA_EXPR}")
    print(f"  - {OUT_PCA_JOINT}")
