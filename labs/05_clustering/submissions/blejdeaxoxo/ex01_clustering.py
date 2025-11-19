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
    url = "https://archive.ics.uci.edu/ml/machine-learning-databases/breast-cancer-wisconsin/wdbc.data"
    columns = ["ID", "Diagnosis"] + [f"Feature_{i}" for i in range(1, 31)]
    df = pd.read_csv(url, header=None, names=columns)

    df = df.drop(columns=["ID"])
    df["Diagnosis"] = df["Diagnosis"].apply(lambda x: 1 if x == "M" else 0)

    X = df.drop(columns=["Diagnosis"])
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    # Director pentru rezultate
    output_dir = Path("labs/05_clustering/submissions/blejdeaxoxo")
    output_dir.mkdir(parents=True, exist_ok=True)

    Z = linkage(X_scaled, method="average")

    plt.figure(figsize=(12, 6))
    dendrogram(Z, leaf_rotation=90)
    plt.title("Hierarchical Clustering Dendrogram")
    plt.xlabel("Samples")
    plt.ylabel("Distance")
    plt.tight_layout()
    plt.savefig(output_dir / "hierarchical_blejdeaxoxo.png")
    plt.close()

    kmeans = KMeans(n_clusters=2, random_state=42)
    df["KMeans_Cluster"] = kmeans.fit_predict(X_scaled)

    pca = PCA(n_components=2)
    X_pca = pca.fit_transform(X_scaled)

    plt.figure(figsize=(8, 6))
    plt.scatter(X_pca[:, 0], X_pca[:, 1], c=df["KMeans_Cluster"], s=20)
    plt.title("K-Means Clustering (PCA 2D)")
    plt.xlabel("PC1")
    plt.ylabel("PC2")
    plt.tight_layout()
    plt.savefig(output_dir / "kmeans_blejdeaxoxo.png")
    plt.close()

    db = DBSCAN(eps=1.5, min_samples=5)
    df["DBSCAN_Cluster"] = db.fit_predict(X_scaled)

    plt.figure(figsize=(8, 6))
    plt.scatter(X_pca[:, 0], X_pca[:, 1], c=df["DBSCAN_Cluster"], s=20)
    plt.title("DBSCAN Clustering (PCA 2D)")
    plt.xlabel("PC1")
    plt.ylabel("PC2")
    plt.tight_layout()
    plt.savefig(output_dir / "dbscan_blejdeaxoxo.png")
    plt.close()

    df[["Diagnosis", "KMeans_Cluster", "DBSCAN_Cluster"]].to_csv(output_dir / "clusters_blejdeaxoxo.csv",index=False)
