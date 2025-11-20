import os
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans, DBSCAN
from sklearn.decomposition import PCA
from scipy.cluster.hierarchy import dendrogram, linkage

handle = "ssmaRRR"
output_dir = Path(f"labs/05_clustering/submissions/{handle}")
output_dir.mkdir(parents=True, exist_ok=True)

# --------------------------------------------------
# 1. Încărcare dataset
# --------------------------------------------------
def load_wdbc_or_fallback():
    """
    Încarcă dataset-ul WDBC de la UCI.
    Dacă nu există conexiune, folosește dataset-ul breast_cancer din sklearn.
    Returnează DataFrame cu Diagnosis (1=M, 0=B) și 30 de features.
    """
    url = "https://archive.ics.uci.edu/ml/machine-learning-databases/breast-cancer-wisconsin/wdbc.data"
    columns = ["ID", "Diagnosis"] + [f"Feature_{i}" for i in range(1, 31)]
    try:
        df = pd.read_csv(url, header=None, names=columns)
        df = df.drop(columns=["ID"])
        df["Diagnosis"] = df["Diagnosis"].map({"M": 1, "B": 0}).astype(int)
        return df
    except Exception:
        from sklearn.datasets import load_breast_cancer
        data = load_breast_cancer()
        X = pd.DataFrame(data.data, columns=data.feature_names)
        rename_map = {col: f"Feature_{i+1}" for i, col in enumerate(X.columns)}
        X = X.rename(columns=rename_map)
        y = pd.Series(data.target, name="Diagnosis")
        y = y.apply(lambda t: 1 if t == 0 else 0)
        df = pd.concat([y, X], axis=1)
        return df

df = load_wdbc_or_fallback()

# --------------------------------------------------
# 2. Preprocesare și standardizare
# --------------------------------------------------
X = df.drop(columns=["Diagnosis"])
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

# --------------------------------------------------
# 3. Hierarchical Clustering
# --------------------------------------------------
Z = linkage(X_scaled, method="average")
plt.figure(figsize=(10, 5))
dendrogram(Z, no_labels=True, count_sort=True)
plt.title("Hierarchical Clustering Dendrogram (average linkage)")
plt.xlabel("Samples")
plt.ylabel("Distance")
plt.tight_layout()
plt.savefig(output_dir / f"hierarchical_{handle}.png", dpi=160)
plt.close()

# --------------------------------------------------
# 4. K-Means Clustering
# --------------------------------------------------
kmeans = KMeans(n_clusters=2, n_init=10, random_state=42)
df["KMeans_Cluster"] = kmeans.fit_predict(X_scaled)

pca = PCA(n_components=2, random_state=42)
X_pca = pca.fit_transform(X_scaled)

plt.figure(figsize=(7, 6))
for label in np.unique(df["KMeans_Cluster"]):
    mask = df["KMeans_Cluster"] == label
    plt.scatter(X_pca[mask, 0], X_pca[mask, 1], label=f"Cluster {label}", s=14)
plt.title("K-Means (K=2) on PCA(2)")
plt.xlabel("PC1")
plt.ylabel("PC2")
plt.legend(fontsize=9)
plt.tight_layout()
plt.savefig(output_dir / f"kmeans_{handle}.png", dpi=160)
plt.close()

# --------------------------------------------------
# 5. DBSCAN Clustering
# --------------------------------------------------
dbscan = DBSCAN(eps=1.5, min_samples=5)
df["DBSCAN_Cluster"] = dbscan.fit_predict(X_scaled)

plt.figure(figsize=(7, 6))
for label in np.unique(df["DBSCAN_Cluster"]):
    mask = df["DBSCAN_Cluster"] == label
    plt.scatter(X_pca[mask, 0], X_pca[mask, 1], label=f"Cluster {label}", s=14)
plt.title("DBSCAN on PCA(2)")
plt.xlabel("PC1")
plt.ylabel("PC2")
plt.legend(fontsize=9)
plt.tight_layout()
plt.savefig(output_dir / f"dbscan_{handle}.png", dpi=160)
plt.close()

# --------------------------------------------------
# 6. Salvare rezultate
# --------------------------------------------------
df[["Diagnosis", "KMeans_Cluster", "DBSCAN_Cluster"]].to_csv(
    output_dir / f"clusters_{handle}.csv", index=False
)

print("Fișiere salvate în:", output_dir)
