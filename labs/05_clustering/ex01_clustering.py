"""
Exercitiu 6 — Clustering pe date de cancer mamar (toy dataset)
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans, DBSCAN
from sklearn.decomposition import PCA
from scipy.cluster.hierarchy import dendrogram, linkage
from pathlib import Path

# ======================
# Handle-ul tau
# ======================
HANDLE = "NistorIuliana03"
output_dir = Path(f"labs/05_clustering/submissions/{HANDLE}")
output_dir.mkdir(parents=True, exist_ok=True)

# ======================
# 1. Load dataset
# ======================
url = "https://archive.ics.uci.edu/ml/machine-learning-databases/breast-cancer-wisconsin/wdbc.data"
columns = ["ID", "Diagnosis"] + [f"Feature_{i}" for i in range(1, 31)]
df = pd.read_csv(url, header=None, names=columns)

# ======================
# 2. Preprocesare
# ======================
df = df.drop(columns=["ID"])
df["Diagnosis"] = df["Diagnosis"].apply(lambda x: 1 if x == "M" else 0)

X = df.drop(columns=["Diagnosis"])

# ======================
# 3. Standardizare
# ======================
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

# ======================
# 4. Hierarchical Clustering
# ======================
Z = linkage(X_scaled, method="average")

plt.figure(figsize=(10, 5))
dendrogram(Z)
plt.title("Hierarchical Clustering Dendrogram")
plt.xlabel("Samples")
plt.ylabel("Distance")
plt.tight_layout()
plt.savefig(output_dir / f"hierarchical_{HANDLE}.png")
plt.close()

# ======================
# 5. K-means clustering
# ======================
kmeans = KMeans(n_clusters=2, random_state=0)
df["KMeans_Cluster"] = kmeans.fit_predict(X_scaled)

pca = PCA(n_components=2)
X_pca = pca.fit_transform(X_scaled)

plt.figure(figsize=(7, 5))
plt.scatter(X_pca[:, 0], X_pca[:, 1], c=df["KMeans_Cluster"], cmap="viridis")
plt.title("K-means (K=2)")
plt.xlabel("PC1")
plt.ylabel("PC2")
plt.tight_layout()
plt.savefig(output_dir / f"kmeans_{HANDLE}.png")
plt.close()

# ======================
# 6. DBSCAN clustering
# ======================
dbscan = DBSCAN(eps=1.5, min_samples=5)
df["DBSCAN_Cluster"] = dbscan.fit_predict(X_scaled)

plt.figure(figsize=(7, 5))
plt.scatter(X_pca[:, 0], X_pca[:, 1], c=df["DBSCAN_Cluster"], cmap="tab10")
plt.title("DBSCAN Clustering")
plt.xlabel("PC1")
plt.ylabel("PC2")
plt.tight_layout()
plt.savefig(output_dir / f"dbscan_{HANDLE}.png")
plt.close()

# ======================
# 7. Salvare rezultate in CSV
# ======================
result_df = df[["Diagnosis", "KMeans_Cluster", "DBSCAN_Cluster"]]
result_df.to_csv(output_dir / f"clusters_{HANDLE}.csv", index=False)

print("\n=== Clustering complete ===")
print("Saved files in:", output_dir)
