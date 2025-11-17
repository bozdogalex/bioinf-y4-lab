import pandas as pd
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans, DBSCAN
from sklearn.decomposition import PCA
from scipy.cluster.hierarchy import dendrogram, linkage
from pathlib import Path

if __name__ == "__main__":
    # TODO 1: Încărcați dataset-ul
    url = "https://archive.ics.uci.edu/ml/machine-learning-databases/breast-cancer-wisconsin/wdbc.data"
    col_names = ["ID", "Diagnosis"] + [f"feature_{i}" for i in range(1, 31)]
    df = pd.read_csv(url, header=None, names=col_names)

    # TODO 2: Preprocesare
    df = df.drop(columns=["ID"])
    df["Diagnosis"] = df["Diagnosis"].map({"M": 1, "B": 0})

    # TODO 3: Standardizare
    X = df.drop(columns=["Diagnosis"])
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    # Director pentru rezultate
    output_dir = Path("labs/05_clustering/submissions/ameliaPaizs")
    output_dir.mkdir(parents=True, exist_ok=True)

    # TODO 4: Hierarchical Clustering
    Z = linkage(X_scaled, method="average")
    plt.figure(figsize=(10, 6))
    dendrogram(Z)
    plt.savefig(output_dir / "hierarchical_ameliaPaizs.png")
    plt.close()

    # TODO 5: K-means Clustering
    kmeans = KMeans(n_clusters=2, n_init=10)
    df["KMeans_Cluster"] = kmeans.fit_predict(X_scaled)

    pca = PCA(n_components=2)
    X_pca = pca.fit_transform(X_scaled)

    plt.figure(figsize=(8, 6))
    plt.scatter(X_pca[:, 0], X_pca[:, 1], c=df["KMeans_Cluster"], cmap="viridis")
    plt.savefig(output_dir / "kmeans_ameliaPaizs.png")
    plt.close()

    # TODO 6: DBSCAN Clustering
    dbscan = DBSCAN(eps=1.5, min_samples=5)
    df["DBSCAN_Cluster"] = dbscan.fit_predict(X_scaled)

    plt.figure(figsize=(8, 6))
    plt.scatter(X_pca[:, 0], X_pca[:, 1], c=df["DBSCAN_Cluster"], cmap="viridis")
    plt.savefig(output_dir / "dbscan_ameliaPaizs.png")
    plt.close()

    # TODO 7: Salvare rezultate
    df[["Diagnosis", "KMeans_Cluster", "DBSCAN_Cluster"]].to_csv(output_dir / "clusters_ameliaPaizs.csv", index=False)
