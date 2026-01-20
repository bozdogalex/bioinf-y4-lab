from pathlib import Path
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt

HANDLE = "numipasaaa" # sau "large" in functie de ce ai generat
SNP_CSV = Path(f"data/work/numipasaaa/lab10/snp_matrix_{HANDLE}.csv")
EXP_CSV = Path(f"data/work/numipasaaa/lab10/expression_matrix_{HANDLE}.csv")

OUT_DIR = Path(f"labs/10_integrative/submissions/{HANDLE}")
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Variabila globala pentru a stoca culorile (dictionar SampleID -> Subtype)
sample_colors = {}

def load_and_normalize(path: Path, is_expression=False) -> pd.DataFrame:
    df = pd.read_csv(path, index_col=0)
    
    # Daca e matricea de expresie, salvam Subtype-urile pentru colorare
    if is_expression and 'Subtype' in df.columns:
        global sample_colors
        # Cream un map: SampleID -> Subtype
        sample_colors = df['Subtype'].to_dict()
        
    # Pastram doar datele numerice pentru PCA
    df_numeric = df.select_dtypes(include=['number'])
    
    # Z-score normalization
    df_norm = (df_numeric - df_numeric.mean()) / df_numeric.std()
    
    # Inlocuim NaN cu 0 (in caz ca deviatia standard e 0)
    df_norm = df_norm.fillna(0)
    
    return df_norm

def plot_pca(df: pd.DataFrame, title: str, out_path: Path) -> None:
    # PCA
    pca = PCA(n_components=2)
    proj = pca.fit_transform(df) # Sklearn asteapta (Samples, Features)
    
    # Pregatire culori
    # Default: albastru
    colors = ['blue'] * len(df)
    
    # Daca avem informatii despre Subtype, coloram in functie de ele
    if sample_colors:
        # Extragem subtipurile in ordinea sample-urilor din df
        subtypes = [sample_colors.get(idx, 'Unknown') for idx in df.index]
        unique_subtypes = list(set(subtypes))
        
        # Mapare simpla string -> numar pentru culori
        color_map = {st: i for i, st in enumerate(unique_subtypes)}
        c_values = [color_map[s] for s in subtypes]
        
        plt.figure(figsize=(8, 6))
        scatter = plt.scatter(proj[:, 0], proj[:, 1], c=c_values, cmap='viridis', alpha=0.7)
        plt.legend(handles=scatter.legend_elements()[0], labels=unique_subtypes, title="Subtype")
    else:
        plt.figure(figsize=(8, 6))
        plt.scatter(proj[:, 0], proj[:, 1], c="blue", alpha=0.7)

    plt.title(title)
    plt.xlabel(f"PC1 ({pca.explained_variance_ratio_[0]:.2%} var)")
    plt.ylabel(f"PC2 ({pca.explained_variance_ratio_[1]:.2%} var)")
    plt.tight_layout()
    plt.savefig(out_path)
    plt.close()
    print(f"Salvat: {out_path}")

# --- Main Logic ---

print("loading data...")
# Specificam is_expression=True doar pentru fisierul de expresie care are 'Subtype'
df_snp = load_and_normalize(SNP_CSV, is_expression=False)
df_exp = load_and_normalize(EXP_CSV, is_expression=True)

# Align samples (Intersectia)
common_samples = df_snp.index.intersection(df_exp.index)
print(f"Probe comune: {len(common_samples)}")

df_snp = df_snp.loc[common_samples]
df_exp = df_exp.loc[common_samples]

# Joint dataset (Concatenare pe coloane -> axis=1)
# Atentie: axis=1 pentru a pune Genele langa SNP-uri pentru aceiasi pacienti
df_joint = pd.concat([df_snp, df_exp], axis=1)

# Plotting
plot_pca(df_snp, "PCA on SNP Data", OUT_DIR / f"pca_snp_{HANDLE}.png")
plot_pca(df_exp, "PCA on Expression Data", OUT_DIR / f"pca_expression_{HANDLE}.png")
plot_pca(df_joint, "PCA on Joint Data (SNP + Expr)", OUT_DIR / f"pca_joint_{HANDLE}.png")