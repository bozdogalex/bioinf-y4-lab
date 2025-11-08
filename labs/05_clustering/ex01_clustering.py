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

if __name__ == "__main__":
    # TODO 1: Încărcați dataset-ul

    # TODO 2: Preprocesare
    # - eliminați coloana ID

    # TODO 3: Standardizare
    X = df.drop(columns=["Diagnosis"])
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    # Director pentru rezultate
    output_dir = Path("labs/05_clustering/submissions/<handle>")
    output_dir.mkdir(parents=True, exist_ok=True)

    # TODO 4: Hierarchical Clustering
    # - folosiți linkage(X_scaled, method="average")
    # - vizualizați cu dendrogram()
    # - salvați imaginea ca hierarchical_<handle>.png

    # TODO 5: K-means Clustering
    # - aplicați KMeans cu K=2
    # - adăugați etichetele în df["KMeans_Cluster"]
    # - reduceți dimensionalitatea cu PCA(n_components=2)
    # - vizualizați și salvați plotul kmeans_<handle>.png

    # TODO 6: DBSCAN Clustering
    # - aplicați DBSCAN (ex: eps=1.5, min_samples=5)
    # - adăugați etichetele în df["DBSCAN_Cluster"]
    # - vizualizați și salvați plotul dbscan_<handle>.png

    # TODO 7: Salvare rezultate
    # salvați un CSV cu coloanele ["Diagnosis", "KMeans_Cluster", "DBSCAN_Cluster"]
    # în clusters_<handle>.csv
