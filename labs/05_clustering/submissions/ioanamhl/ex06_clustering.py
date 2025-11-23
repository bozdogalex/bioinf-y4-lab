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
# from scipy.cluster.hierarchy import dendrogram, linkage
# from sklearn.cluster import KMeans, DBSCAN
# from sklearn.decomposition import PCA

from scipy.cluster.hierarchy import dendrogram, linkage
from sklearn.cluster import KMeans, DBSCAN
from sklearn.decomposition import PCA

import requests
from io import StringIO
import urllib3

urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)


if __name__ == "__main__":

    # TODO 1: Încărcați dataset-ul
    url = "https://archive.ics.uci.edu/ml/machine-learning-databases/breast-cancer-wisconsin/wdbc.data"
    columns = ["ID", "Diagnosis"] + [f"Feature_{i}" for i in range(1, 31)]

    response = requests.get(url, verify=False, timeout=30)
    response.raise_for_status()  

    df = pd.read_csv(StringIO(response.text), header=None, names=columns)

    # TODO 2: Preprocesare
    # - eliminați coloana ID
    df = df.drop(columns=["ID"])
    df["Diagnosis"] = df["Diagnosis"].apply(lambda x: 1 if x == "M" else 0)

    # TODO 3: Standardizare
    X = df.drop(columns=["Diagnosis"])
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

   # Director pentru rezultate
    output_dir = Path(f"labs/05_clustering/submissions/ioanamhl")
    output_dir.mkdir(parents=True, exist_ok=True)

    # TODO 4: Hierarchical Clustering
    # - folosiți linkage(X_scaled, method="average")
    # - vizualizați cu dendrogram()
    # - salvați imaginea ca hierarchical_<handle>.png

    Z = linkage(X_scaled, method="average")
    plt.figure(figsize=(10, 6))
    dendrogram(Z)
    plt.title("Hierarchical Clustering — Dendrogramă")
    plt.xlabel("Sample Index")
    plt.ylabel("Distanță")
    plt.tight_layout()
    plt.savefig(output_dir / f"hierarchical_ioanamhl.png", dpi=150)
    plt.close()

    # TODO 5: K-means Clustering
    # - aplicați KMeans cu K=2
    # - adăugați etichetele în df["KMeans_Cluster"]
    # - reduceți dimensionalitatea cu PCA(n_components=2)
    # - vizualizați și salvați plotul kmeans_<handle>.png
    kmeans = KMeans(n_clusters=2, random_state=0)
    kmeans_labels = kmeans.fit_predict(X_scaled)
    df["KMeans_Cluster"] = kmeans_labels

    pca_kmeans = PCA(n_components=2)
    X_pca_kmeans = pca_kmeans.fit_transform(X_scaled)

    plt.figure(figsize=(8, 6))
    plt.scatter(
        X_pca_kmeans[:, 0],
        X_pca_kmeans[:, 1],
        c=kmeans_labels,
        s=40,
        alpha=0.7,
    )
    plt.title("K-means Clustering (K=2) — PCA Visualization")
    plt.xlabel("PCA 1")
    plt.ylabel("PCA 2")
    plt.tight_layout()
    plt.savefig(output_dir / f"kmeans_ioanamhl.png", dpi=150)
    plt.close()

    # TODO 6: DBSCAN Clustering
    # - aplicați DBSCAN (ex: eps=1.5, min_samples=5)
    # - adăugați etichetele în df["DBSCAN_Cluster"]
    # - vizualizați și salvați plotul dbscan_<handle>.png
    dbscan = DBSCAN(eps=1.5, min_samples=5)
    dbscan_labels = dbscan.fit_predict(X_scaled)
    df["DBSCAN_Cluster"] = dbscan_labels  # -1 = noise

    pca_dbscan = PCA(n_components=2)
    X_pca_dbscan = pca_dbscan.fit_transform(X_scaled)

    plt.figure(figsize=(8, 6))
    plt.scatter(
        X_pca_dbscan[:, 0],
        X_pca_dbscan[:, 1],
        c=dbscan_labels,
        s=40,
        alpha=0.7,
    )
    plt.title("DBSCAN Clustering — PCA Visualization")
    plt.xlabel("PCA 1")
    plt.ylabel("PCA 2")
    plt.tight_layout()
    plt.savefig(output_dir / f"dbscan_ioanamhl.png", dpi=150)
    plt.close()

    # TODO 7: Salvare rezultate
    # salvați un CSV cu coloanele ["Diagnosis", "KMeans_Cluster", "DBSCAN_Cluster"]
    # în clusters_<handle>.csv
    df[["Diagnosis", "KMeans_Cluster", "DBSCAN_Cluster"]].to_csv(
        output_dir / f"clusters_ioanamhl.csv",
        index=False,
    )

    print(f"Toate fișierele au fost generate în: {output_dir}")
