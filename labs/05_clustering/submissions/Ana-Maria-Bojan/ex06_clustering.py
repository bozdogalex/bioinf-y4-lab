# clustering_breast_cancer.py
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from pathlib import Path

# clustering & PCA
from scipy.cluster.hierarchy import dendrogram, linkage
from sklearn.cluster import KMeans, DBSCAN
from sklearn.decomposition import PCA
import numpy as np


RANDOM_STATE = 42
HANDLE = "Ana-Maria-Bojan"  

def main():
    # 1) Încărcați dataset-ul WDBC
    url = "https://archive.ics.uci.edu/ml/machine-learning-databases/breast-cancer-wisconsin/wdbc.data"
    columns = ["ID", "Diagnosis"] + [f"Feature_{i}" for i in range(1, 31)]
    df = pd.read_csv(url, header=None, names=columns)

    # 2) Preprocesare
    df = df.drop(columns=["ID"])
    df["Diagnosis"] = df["Diagnosis"].apply(lambda x: 1 if x == "M" else 0)

    # 3) Standardizare
    X = df.drop(columns=["Diagnosis"])
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    # Director pentru rezultate
    output_dir = Path(f"submissions/{HANDLE}")
    output_dir.mkdir(parents=True, exist_ok=True)

    # 4) Hierarchical Clustering (linkage + dendrogram)
    Z = linkage(X_scaled, method="average")
    plt.figure(figsize=(12, 6))
    dendrogram(Z, truncate_mode="level", p=25, leaf_rotation=90., leaf_font_size=8.)
    plt.title("Hierarchical Clustering Dendrogram (average linkage)")
    plt.xlabel("Sample index (truncated)")
    plt.ylabel("Distance")
    hierarchical_path = output_dir / f"hierarchical_{HANDLE}.png"
    plt.tight_layout()
    plt.savefig(hierarchical_path, dpi=150)
    plt.close()
    print(f"Saved dendrogram: {hierarchical_path}")

    # 5) K-means (K=2) + PCA vizualizare
    kmeans = KMeans(n_clusters=2, random_state=RANDOM_STATE, n_init=10)
    kmeans_labels = kmeans.fit_predict(X_scaled)
    df["KMeans_Cluster"] = kmeans_labels

    pca = PCA(n_components=2, random_state=RANDOM_STATE)
    X_pca = pca.fit_transform(X_scaled)

    plt.figure(figsize=(8, 6))
    # scatter color by cluster
    for cluster in np.unique(kmeans_labels):
        mask = kmeans_labels == cluster
        plt.scatter(X_pca[mask, 0], X_pca[mask, 1], label=f"KMeans {cluster}", alpha=0.7, edgecolor='k', s=40)
    plt.title("KMeans (K=2) clusters projected with PCA")
    plt.xlabel("PC1")
    plt.ylabel("PC2")
    plt.legend()
    kmeans_path = output_dir / f"kmeans_{HANDLE}.png"
    plt.tight_layout()
    plt.savefig(kmeans_path, dpi=150)
    plt.close()
    print(f"Saved kmeans plot: {kmeans_path}")

    # 6) DBSCAN (eps=1.5, min_samples=5) + PCA vizualizare
    dbscan = DBSCAN(eps=1.5, min_samples=5)
    db_labels = dbscan.fit_predict(X_scaled)
    df["DBSCAN_Cluster"] = db_labels

    plt.figure(figsize=(8, 6))
    unique_db = np.unique(db_labels)
    for cluster in unique_db:
        mask = db_labels == cluster
        label = f"Noise" if cluster == -1 else f"Cluster {cluster}"
        plt.scatter(X_pca[mask, 0], X_pca[mask, 1], label=label, alpha=0.7, edgecolor='k', s=40)
    plt.title("DBSCAN clusters projected with PCA")
    plt.xlabel("PC1")
    plt.ylabel("PC2")
    plt.legend()
    dbscan_path = output_dir / f"dbscan_{HANDLE}.png"
    plt.tight_layout()
    plt.savefig(dbscan_path, dpi=150)
    plt.close()
    print(f"Saved DBSCAN plot: {dbscan_path}")

    # 7) Salvare rezultate CSV
    out_csv = output_dir / f"clusters_{HANDLE}.csv"
    df_out = df[["Diagnosis", "KMeans_Cluster", "DBSCAN_Cluster"]].copy()
    df_out.to_csv(out_csv, index=False)
    print(f"Saved clusters CSV: {out_csv}")

if __name__ == "__main__":
    main()
