import pandas as pd
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans, DBSCAN
from sklearn.decomposition import PCA
from scipy.cluster.hierarchy import linkage, dendrogram
from pathlib import Path

if __name__ == "__main__":
    # --- 1. Încărcare dataset ---
    url = "https://archive.ics.uci.edu/ml/machine-learning-databases/breast-cancer-wisconsin/wdbc.data"
    columns = ["ID", "Diagnosis"] + [f"Feature_{i}" for i in range(1,31)]
    df = pd.read_csv(url, header=None, names=columns)

    # --- 2. Preprocesare ---
    df = df.drop(columns=["ID"])
    df["Diagnosis"] = df["Diagnosis"].apply(lambda x: 1 if x == "M" else 0)

    # --- 3. Standardizare ---
    X = df.drop(columns=["Diagnosis"])
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    # --- Director pentru rezultate ---
    handle = "MadalinaDenisa"
    output_dir = Path(f"labs/05_clustering/submissions/{handle}")
    output_dir.mkdir(parents=True, exist_ok=True)

    # --- 4. Hierarchical Clustering ---
    linked = linkage(X_scaled, method="average")
    plt.figure(figsize=(10,6))
    dendrogram(linked, labels=df.index, leaf_rotation=90)
    plt.title("Hierarchical Clustering Dendrogram")
    plt.savefig(output_dir / f"hierarchical_{handle}.png")
    plt.close()

    # --- 5. K-means Clustering ---
    kmeans = KMeans(n_clusters=2, random_state=42)
    df["KMeans_Cluster"] = kmeans.fit_predict(X_scaled)

    pca = PCA(n_components=2)
    X_pca = pca.fit_transform(X_scaled)
    plt.figure(figsize=(8,6))
    plt.scatter(X_pca[:,0], X_pca[:,1], c=df["KMeans_Cluster"], cmap="viridis", s=50)
    plt.title("K-means Clustering PCA")
    plt.xlabel("PC1")
    plt.ylabel("PC2")
    plt.savefig(output_dir / f"kmeans_{handle}.png")
    plt.close()

    # --- 6. DBSCAN Clustering ---
    dbscan = DBSCAN(eps=1.5, min_samples=5)
    df["DBSCAN_Cluster"] = dbscan.fit_predict(X_scaled)

    plt.figure(figsize=(8,6))
    plt.scatter(X_pca[:,0], X_pca[:,1], c=df["DBSCAN_Cluster"], cmap="tab10", s=50)
    plt.title("DBSCAN Clustering PCA")
    plt.xlabel("PC1")
    plt.ylabel("PC2")
    plt.savefig(output_dir / f"dbscan_{handle}.png")
    plt.close()

    # --- 7. Salvare rezultate ---
    df[["Diagnosis","KMeans_Cluster","DBSCAN_Cluster"]].to_csv(
        output_dir / f"clusters_{handle}.csv", index=False
    )
