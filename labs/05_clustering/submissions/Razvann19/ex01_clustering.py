"""
Exercițiu 6 — Clustering pe date de cancer mamar (toy dataset)
"""

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans, DBSCAN

from scipy.cluster.hierarchy import dendrogram, linkage

HANDLE = "Razvann19"

if __name__ == "__main__":
    url = "https://archive.ics.uci.edu/ml/machine-learning-databases/breast-cancer-wisconsin/wdbc.data"
    columns = ["ID", "Diagnosis"] + [f"Feature_{i}" for i in range(1, 31)]
    df = pd.read_csv(url, header=None, names=columns)

    df = df.drop(columns=["ID"])
    df["Diagnosis"] = df["Diagnosis"].apply(lambda x: 1 if x == "M" else 0)
    y = df["Diagnosis"].to_numpy()

    X = df.drop(columns=["Diagnosis"])
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    output_dir = Path(f"labs/05_clustering/submissions/{HANDLE}")
    output_dir.mkdir(parents=True, exist_ok=True)

    Z = linkage(X_scaled, method="average")
    plt.figure(figsize=(10, 5))
    dendrogram(Z, no_labels=True, distance_sort="ascending", color_threshold=None)
    plt.title("Hierarchical Clustering (average linkage) — WDBC")
    plt.xlabel("Samples")
    plt.ylabel("Distance")
    plt.tight_layout()
    plt.savefig(output_dir / f"hierarchical_{HANDLE}.png", dpi=200)
    plt.close()

    # 5) K-means (K=2) + PCA vizualizare
    kmeans = KMeans(n_clusters=2, random_state=42, n_init=10)
    kmeans_labels = kmeans.fit_predict(X_scaled)
    df["KMeans_Cluster"] = kmeans_labels

    pca = PCA(n_components=2, random_state=42)
    X_pca = pca.fit_transform(X_scaled)

    plt.figure(figsize=(7, 6))
    sc = plt.scatter(X_pca[:, 0], X_pca[:, 1], c=kmeans_labels, s=40, alpha=0.8, cmap="viridis")
    plt.title("K-means (K=2) pe PCA — WDBC")
    plt.xlabel("PCA 1")
    plt.ylabel("PCA 2")
    cbar = plt.colorbar(sc)
    cbar.set_label("Cluster")
    plt.tight_layout()
    plt.savefig(output_dir / f"kmeans_{HANDLE}.png", dpi=200)
    plt.close()

    dbscan = DBSCAN(eps=3.0, min_samples=5, n_jobs=-1)
    db_labels = dbscan.fit_predict(X_scaled)
    df["DBSCAN_Cluster"] = db_labels  # -1 = noise

    plt.figure(figsize=(7, 6))
    sc = plt.scatter(X_pca[:, 0], X_pca[:, 1], c=db_labels, s=40, alpha=0.8, cmap="viridis")
    plt.title("DBSCAN pe PCA — WDBC")
    plt.xlabel("PCA 1")
    plt.ylabel("PCA 2")
    cbar = plt.colorbar(sc)
    cbar.set_label("Cluster (-1 = noise)")
    plt.tight_layout()
    plt.savefig(output_dir / f"dbscan_{HANDLE}.png", dpi=200)
    plt.close()

    df[["Diagnosis", "KMeans_Cluster", "DBSCAN_Cluster"]].to_csv(
        output_dir / f"clusters_{HANDLE}.csv", index=False
    )

    print("Saved:",
          f"\n - {output_dir / f'hierarchical_{HANDLE}.png'}",
          f"\n - {output_dir / f'kmeans_{HANDLE}.png'}",
          f"\n - {output_dir / f'dbscan_{HANDLE}.png'}",
          f"\n - {output_dir / f'clusters_{HANDLE}.csv'}")
    unique_db, counts_db = np.unique(db_labels, return_counts=True)
    print("DBSCAN cluster sizes:", dict(zip(unique_db, counts_db)))
