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
    # - transformați Diagnosis: M → 1, B → 0
    df = df.drop(columns=["ID"])
    df["Diagnosis"] = df["Diagnosis"].apply(lambda x: 1 if x == "M" else 0)

    # TODO 3: Standardizare
    X = df.drop(columns=["Diagnosis"])
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    output_dir = Path("labs/05_clustering/submissions/marabonea16")
    output_dir.mkdir(parents=True, exist_ok=True)

    # TODO 4: Hierarchical Clustering
    Z = linkage(X_scaled, method="average")
    plt.figure(figsize=(10, 7))
    dendrogram(Z)
    plt.title("Hierarchical Clustering Dendrogram")
    plt.xlabel("Samples")
    plt.ylabel("Distance")
    plt.savefig("labs/05_clustering/submissions/marabonea16/hierarchical_marabonea16.png")
    plt.close()

    # TODO 5: K-means Clustering
    kmeans = KMeans(n_clusters=2, random_state=0)
    kmeans_labels = kmeans.fit_predict(X_scaled)
    df["KMeans_Cluster"] = kmeans_labels
    pca = PCA(n_components=2)
    X_pca = pca.fit_transform(X_scaled)
    plt.figure(figsize=(8, 6))
    plt.scatter(X_pca[:, 0], X_pca[:, 1], c=kmeans_labels, cmap="viridis", s=50, alpha=0.7)
    plt.title("K-means Clustering (K=2) pe date PCA")
    plt.xlabel("PCA 1")
    plt.ylabel("PCA 2")
    plt.colorbar(label="Cluster")
    plt.savefig("labs/05_clustering/submissions/marabonea16/kmeans_marabonea16.png")
    plt.close()

    # TODO 6: DBSCAN Clustering
    dbscan = DBSCAN(eps=1.5, min_samples=5)
    dbscan_labels = dbscan.fit_predict(X_scaled)
    df["DBSCAN_Cluster"] = dbscan_labels
    plt.figure(figsize=(8, 6))
    plt.scatter(X_pca[:, 0], X_pca[:, 1], c=dbscan_labels, cmap="viridis", s=50, alpha=0.7)
    plt.title("DBSCAN Clustering pe date PCA")
    plt.xlabel("PCA 1")
    plt.ylabel("PCA 2")
    plt.colorbar(label="Cluster")
    plt.savefig("labs/05_clustering/submissions/marabonea16/dbscan_marabonea16.png")
    plt.close()

    # TODO 7: Salvare rezultate
    df[["Diagnosis", "KMeans_Cluster", "DBSCAN_Cluster"]].to_csv(output_dir / "clusters_marabonea16.csv", index=False)
