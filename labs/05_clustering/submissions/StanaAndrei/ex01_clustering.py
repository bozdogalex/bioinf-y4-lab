import pandas as pd
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from pathlib import Path

from scipy.cluster.hierarchy import dendrogram, linkage
from sklearn.cluster import KMeans, DBSCAN
from sklearn.decomposition import PCA

if __name__ == "__main__":
    handle = "StanaAndrei"
    output_dir = Path(f"labs/05_clustering/submissions/{handle}")
    output_dir.mkdir(parents=True, exist_ok=True)

    # TODO 1: Încărcați dataset-ul
    url = "https://archive.ics.uci.edu/ml/machine-learning-databases/breast-cancer-wisconsin/wdbc.data"
    columns = ["ID", "Diagnosis"] + [f"Feature_{i}" for i in range(1, 31)]
    df = pd.read_csv(url, header=None, names=columns)

    # TODO 2: Preprocesare
    # - eliminați coloana ID
    # - transformați Diagnosis: M → 1, B → 0
    df = df.drop(columns=["ID"])
    df["Diagnosis"] = df["Diagnosis"].apply(lambda x: 1 if x == "M" else 0)

    # TODO 3: Standardizare
    X = df.drop(columns=["Diagnosis"])
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    # TODO 4: Hierarchical Clustering
    # - folosiți linkage(X_scaled, method="average")
    # - vizualizați cu dendrogram()
    # - salvați imaginea ca hierarchical_<handle>.png
    Z = linkage(X_scaled, method="average")
    plt.figure(figsize=(14, 6))
    dendrogram(Z, no_labels=True, count_sort=True)
    plt.title("Hierarchical Clustering Dendrogram (average linkage)")
    plt.xlabel("Sample index (no labels shown)")
    plt.ylabel("Distance")
    plt.tight_layout()
    dendro_path = output_dir / f"hierarchical_{handle}.png"
    plt.savefig(dendro_path, dpi=150)
    plt.close()
    print(f"Saved dendrogram to: {dendro_path}")

    # TODO 5: K-means Clustering
    # - aplicați KMeans cu K=2
    # - adăugați etichetele în df["KMeans_Cluster"]
    # - reduceți dimensionalitatea cu PCA(n_components=2)
    # - vizualizați și salvați plotul kmeans_<handle>.png
    kmeans = KMeans(n_clusters=2, random_state=42, n_init=10)
    kmeans_labels = kmeans.fit_predict(X_scaled)
    df["KMeans_Cluster"] = kmeans_labels

    pca = PCA(n_components=2, random_state=42)
    X_pca = pca.fit_transform(X_scaled)

    plt.figure(figsize=(8, 6))
    for label in sorted(df["KMeans_Cluster"].unique()):
        mask = df["KMeans_Cluster"] == label
        plt.scatter(X_pca[mask, 0], X_pca[mask, 1], label=f"Cluster {label}", alpha=0.7, edgecolors="k", linewidths=0.3)
    plt.title("KMeans (K=2) clusters (PCA 2D)")
    plt.xlabel("PC1")
    plt.ylabel("PC2")
    plt.legend()
    plt.tight_layout()
    kmeans_path = output_dir / f"kmeans_{handle}.png"
    plt.savefig(kmeans_path, dpi=150)
    plt.close()
    print(f"Saved KMeans plot to: {kmeans_path}")

    # TODO 6: DBSCAN Clustering
    # - aplicați DBSCAN (ex: eps=1.5, min_samples=5)
    # - adăugați etichetele în df["DBSCAN_Cluster"]
    # - vizualizați și salvați plotul dbscan_<handle>.png
    dbscan = DBSCAN(eps=1.5, min_samples=5)
    dbscan_labels = dbscan.fit_predict(X_scaled)
    df["DBSCAN_Cluster"] = dbscan_labels

    plt.figure(figsize=(8, 6))
    unique_labels = sorted(set(dbscan_labels))
    for label in unique_labels:
        mask = dbscan_labels == label
        if label == -1:
            label_name = "Noise (-1)"
            marker = "x"
        else:
            label_name = f"Cluster {label}"
            marker = "o"
        plt.scatter(X_pca[mask, 0], X_pca[mask, 1], label=label_name, alpha=0.7, marker=marker, edgecolors="k", linewidths=0.3)
    plt.title("DBSCAN clusters (eps=1.5, min_samples=5) (PCA 2D)")
    plt.xlabel("PC1")
    plt.ylabel("PC2")
    plt.legend()
    plt.tight_layout()
    dbscan_path = output_dir / f"dbscan_{handle}.png"
    plt.savefig(dbscan_path, dpi=150)
    plt.close()
    print(f"Saved DBSCAN plot to: {dbscan_path}")

    # TODO 7: Salvare rezultate
    # salvați un CSV cu coloanele ["Diagnosis", "KMeans_Cluster", "DBSCAN_Cluster"]
    result_df = df[["Diagnosis", "KMeans_Cluster", "DBSCAN_Cluster"]].copy()
    csv_path = output_dir / f"clusters_{handle}.csv"
    result_df.to_csv(csv_path, index=False)
    print(f"Saved clusters CSV to: {csv_path}")

    # Summary prints
    print("\nSummary:")
    print("Diagnosis distribution:")
    print(df["Diagnosis"].value_counts())
    print("\nKMeans cluster counts:")
    print(df["KMeans_Cluster"].value_counts())
    print("\nDBSCAN cluster counts (includes -1 for noise):")
    print(df["DBSCAN_Cluster"].value_counts())
