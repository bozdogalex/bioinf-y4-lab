"""
Exercițiu 6 — Clustering pe date de cancer mamar (toy dataset)

Instrucțiuni:
1. Încărcați dataset-ul WDBC (breast cancer) de pe UCI Repository.
2. Preprocesați datele: eliminați coloanele irelevante și transformați diagnosticul în valori numerice.
3. Standardizați datele.
4. Implementați și vizualizați clustering-ul folosind:
   - Hierarchical clustering (dendrogramă),
   - K-means (K=2, PCA vizualizare),
   - DBSCAN (PCA vizualizare).
5. Salvați rezultatele în folderul submissions/<handle>/:
   - clusters_<handle>.csv
   - hierarchical_<handle>.png
   - kmeans_<handle>.png
   - dbscan_<handle>.png
"""

import pandas as pd
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from pathlib import Path

# Opțional: puteți importa deja funcțiile necesare
from scipy.cluster.hierarchy import dendrogram, linkage
from sklearn.cluster import KMeans, DBSCAN
from sklearn.decomposition import PCA

if __name__ == "__main__":
    # TODO 1: Încărcați dataset-ul
    url = "https://archive.ics.uci.edu/ml/machine-learning-databases/breast-cancer-wisconsin/wdbc.data"
    columns = ["ID", "Diagnosis"] + [f"Feature_{i}" for i in range(1, 31)]
    df = pd.read_csv(url, header=None, names=columns)

    # TODO 2: Preprocesare
    # - eliminați coloana ID
    df = df.drop(columns=["ID"])
    df["Diagnosis"] = df["Diagnosis"].map({"M": 1, "B": 0})

    # TODO 3: Standardizare
    X = df.drop(columns=["Diagnosis"])
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    # Director pentru rezultate
    output_dir = Path("labs/05_clustering/submissions/numipasaaa")
    output_dir.mkdir(parents=True, exist_ok=True)

    # TODO 4: Hierarchical Clustering
    # - folosiți linkage(X_scaled, method="average")
    # - vizualizați cu dendrogram()
    # - salvați imaginea ca hierarchical_<handle>.png

    linked = linkage(X_scaled, method="average")

    plt.figure(figsize=(12, 8))
    dendrogram(linked, truncate_mode='lastp', p=30)
    plt.title("Hierarchical Clustering Dendrogram")
    plt.xlabel("Samples")
    plt.ylabel("Distance")
    plt.savefig(output_dir / "hierarchical_numipasaaa.png")
    plt.close()

    # TODO 5: K-means Clustering
    # - aplicați KMeans cu K=2
    # - adăugați etichetele în df["KMeans_Cluster"]
    # - reduceți dimensionalitatea cu PCA(n_components=2)
    # - vizualizați și salvați plotul kmeans_<handle>.png

    kmeans = KMeans(n_clusters=2)
    df["KMeans_Cluster"] = kmeans.fit_predict(X_scaled)

    pca = PCA(n_components=2)
    X_pca = pca.fit_transform(X_scaled)

    plt.figure(figsize=(12, 8))
    scatter = plt.scatter(X_pca[:, 0], X_pca[:, 1], c=df["KMeans_Cluster"], cmap="viridis", marker="o")
    plt.title("K-means Clustering (PCA)")
    plt.xlabel("PCA Component 1")
    plt.ylabel("PCA Component 2")
    plt.legend(*scatter.legend_elements(), title="Clusters")
    plt.savefig(output_dir / "kmeans_numipasaaa.png")
    plt.close()
    
    # TODO 6: DBSCAN Clustering
    # - aplicați DBSCAN (ex: eps=1.5, min_samples=5)
    # - adăugați etichetele în df["DBSCAN_Cluster"]
    # - vizualizați și salvați plotul dbscan_<handle>.png


    dbscan = DBSCAN(eps=1.5, min_samples=5)
    df["DBSCAN_Cluster"] = dbscan.fit_predict(X_scaled)

    plt.figure(figsize=(12, 8))
    scatter = plt.scatter(X_pca[:, 0], X_pca[:, 1], c=df["DBSCAN_Cluster"], cmap="viridis", marker="o")
    plt.title("DBSCAN Clustering (PCA)")
    plt.xlabel("PCA Component 1")
    plt.ylabel("PCA Component 2")
    plt.legend(*scatter.legend_elements(), title="Clusters")
    plt.savefig(output_dir / "dbscan_numipasaaa.png")
    plt.close()

    # TODO 7: Salvare rezultate
    # salvați un CSV cu coloanele ["Diagnosis", "KMeans_Cluster", "DBSCAN_Cluster"]
    # în clusters_<handle>.csv
    df[["Diagnosis", "KMeans_Cluster", "DBSCAN_Cluster"]].to_csv(output_dir / "clusters_numipasaaa.csv", index=False)
