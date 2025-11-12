"""
Exercițiul 6 — Clustering pe date de cancer mamar (toy dataset)
Autor: <github_handle>
"""

import pandas as pd
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans, DBSCAN
from scipy.cluster.hierarchy import dendrogram, linkage
from pathlib import Path

if __name__ == "__main__":
    handle = "MarioCld"

    url = "https://archive.ics.uci.edu/ml/machine-learning-databases/breast-cancer-wisconsin/wdbc.data"
    columns = ["ID", "Diagnosis"] + [f"Feature_{i}" for i in range(1, 31)]
    df = pd.read_csv(url, header=None, names=columns)

    df = df.drop(columns=["ID"])
    df["Diagnosis"] = df["Diagnosis"].apply(lambda x: 1 if x == "M" else 0)

    X = df.drop(columns=["Diagnosis"])
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    output_dir = Path(f"labs/05_clustering/submissions/{handle}")
    output_dir.mkdir(parents=True, exist_ok=True)

    linkage_matrix = linkage(X_scaled, method="average")

    plt.figure(figsize=(10, 5))
    dendrogram(linkage_matrix)
    plt.title("Hierarchical Clustering Dendrogram")
    plt.xlabel("Sample index")
    plt.ylabel("Distance")
    plt.tight_layout()
    plt.savefig(output_dir / f"hierarchical_{handle}.png")
    plt.close()

    kmeans = KMeans(n_clusters=2, random_state=42)
    df["KMeans_Cluster"] = kmeans.fit_predict(X_scaled)

    pca = PCA(n_components=2)
    X_pca = pca.fit_transform(X_scaled)

    plt.figure(figsize=(7, 6))
    plt.scatter(X_pca[:, 0], X_pca[:, 1], c=df["KMeans_Cluster"], cmap="viridis", s=40)
    plt.title("K-means Clustering (K=2, PCA vizualizare)")
    plt.xlabel("PC1")
    plt.ylabel("PC2")
    plt.tight_layout()
    plt.savefig(output_dir / f"kmeans_{handle}.png")
    plt.close()

    dbscan = DBSCAN(eps=1.5, min_samples=5)
    df["DBSCAN_Cluster"] = dbscan.fit_predict(X_scaled)

    plt.figure(figsize=(7, 6))
    plt.scatter(X_pca[:, 0], X_pca[:, 1], c=df["DBSCAN_Cluster"], cmap="tab10", s=40)
    plt.title("DBSCAN Clustering (PCA vizualizare)")
    plt.xlabel("PC1")
    plt.ylabel("PC2")
    plt.tight_layout()
    plt.savefig(output_dir / f"dbscan_{handle}.png")
    plt.close()

    # === 7. Salvare rezultate ===
    results = df[["Diagnosis", "KMeans_Cluster", "DBSCAN_Cluster"]]
    results.to_csv(output_dir / f"clusters_{handle}.csv", index=False)
