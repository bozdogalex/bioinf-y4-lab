"""
Exercise 10 — PCA Single-Omics vs Joint
Handle: banmepls
"""

from pathlib import Path
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt

HANDLE = "banmepls"

SNP_CSV = Path(f"data/work/{HANDLE}/lab10/snp_matrix_{HANDLE}.csv")
EXP_CSV = Path(f"data/work/{HANDLE}/lab10/expression_matrix_{HANDLE}.csv")

OUT_DIR = Path(f"labs/10_integrative/submissions/{HANDLE}")
OUT_DIR.mkdir(parents=True, exist_ok=True)

OUT_PCA_SNP = OUT_DIR / f"pca_snp_{HANDLE}.png"
OUT_PCA_EXP = OUT_DIR / f"pca_expression_{HANDLE}.png"
OUT_PCA_JOINT = OUT_DIR / f"pca_joint_{HANDLE}.png"


def load_and_align_data(snp_path: Path, exp_path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    df_snp = pd.read_csv(snp_path, index_col=0)
    
    df_exp = pd.read_csv(exp_path, index_col=0)
    
    # Probe comune (coloanele sunt probele)
    common_samples = df_snp.columns.intersection(df_exp.columns)
    if len(common_samples) == 0:
        raise ValueError("Nicio probă comună între SNP și Expression!")
    
    df_snp = df_snp[common_samples]
    df_exp = df_exp[common_samples]
    
    return df_snp, df_exp


def zscore_normalize(df: pd.DataFrame) -> pd.DataFrame:
    """
    Z-score normalization (media 0, deviație std 1) pe coloane
    """
    return (df - df.mean()) / df.std()


def run_pca_and_plot(df: pd.DataFrame, title: str, out_png: Path) -> None:
    """
    Rulează PCA pe 2 componente și salvează scatter plot
    """
    pca = PCA(n_components=2)
    proj = pca.fit_transform(df.T)  # Transpunem ca să fie probe pe rânduri
    
    explained_var = pca.explained_variance_ratio_
    
    plt.figure(figsize=(8, 6))
    plt.scatter(proj[:, 0], proj[:, 1], c="blue", s=50, alpha=0.7)
    plt.title(f"{title}\n(PC1: {explained_var[0]:.3f}, PC2: {explained_var[1]:.3f})")
    plt.xlabel("PC1")
    plt.ylabel("PC2")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    plt.close()


if __name__ == "__main__":
    # Verificare fișiere
    if not SNP_CSV.exists():
        raise FileNotFoundError(f"SNP fișier lipsește: {SNP_CSV}")
    if not EXP_CSV.exists():
        raise FileNotFoundError(f"Expression fișier lipsește: {EXP_CSV}")
    
    df_snp, df_exp = load_and_align_data(SNP_CSV, EXP_CSV)
    
    df_snp_norm = zscore_normalize(df_snp)
    df_exp_norm = zscore_normalize(df_exp)
    
    df_joint = pd.concat([df_snp_norm, df_exp_norm], axis=0)
    
    run_pca_and_plot(df_snp_norm, "PCA - SNP only", OUT_PCA_SNP)
    
    run_pca_and_plot(df_exp_norm, "PCA - Expression only", OUT_PCA_EXP)
    
    run_pca_and_plot(df_joint, "PCA - Joint Multi-Omics", OUT_PCA_JOINT)
    
    print(f"[OK] Figurile salvate în: {OUT_DIR}")