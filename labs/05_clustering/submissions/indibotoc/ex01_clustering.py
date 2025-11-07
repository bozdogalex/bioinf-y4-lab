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
   df.head()
   # TODO 2: Preprocesare
   # - eliminați coloana ID
   # - transformați Diagnosis: M → 1, B → 0
   df = df.drop(columns=["ID"])
   df["Diagnosis"] = df["Diagnosis"].apply(lambda x: 1 if x == "M" else 0)

   # TODO 3: Standardizare
   X = df.drop(columns=["Diagnosis"])
   scaler = StandardScaler()
   X_scaled = scaler.fit_transform(X)

   # Director pentru rezultate
   handle = "indibotoc"
   output_dir = Path("labs/05_clustering/submissions/{handle}")
   output_dir.mkdir(parents=True, exist_ok=True)

   # TODO 4: Hierarchical Clustering
   # - folosiți linkage(X_scaled, method="average")
   # - vizualizați cu dendrogram()
   # - salvați imaginea ca hierarchical_<handle>.png
   Z = linkage(X_scaled, method="average")

   plt.figure(figsize=(12, 5))
   dendrogram(Z, no_labels=True, color_threshold=None)
   plt.title("Hierarchical clustering — Dendrogram")
   plt.xlabel("Samples")
   plt.ylabel("Distance")
   plt.tight_layout()
   plt.savefig(output_dir / f"hierarchical_{handle}.png", dpi=150)
   plt.close()

   # TODO 5: K-means Clustering
   # - aplicați KMeans cu K=2
   # - adăugați etichetele în df["KMeans_Cluster"]
   # - reduceți dimensionalitatea cu PCA(n_components=2)
   # - vizualizați și salvați plotul kmeans_<handle>.png

   kmeans = KMeans(n_clusters=2, random_state=42, n_init="auto")
   kmeans_labels = kmeans.fit_predict(X_scaled)
   df["KMeans_Cluster"] = kmeans_labels

   # PCA 2D
   pca = PCA(n_components=2, random_state=42)
   X_pca = pca.fit_transform(X_scaled)

   plt.figure(figsize=(6, 5))
   plt.scatter(X_pca[:, 0], X_pca[:, 1], c=kmeans_labels, s=20)
   plt.title("K-Means (K=2) — PCA 2D")
   plt.xlabel("PC1")
   plt.ylabel("PC2")
   plt.tight_layout()
   plt.savefig(output_dir / f"kmeans_{handle}.png", dpi=150)
   plt.close()


    # TODO 6: DBSCAN Clustering
    # - aplicați DBSCAN (ex: eps=1.5, min_samples=5)
    # - adăugați etichetele în df["DBSCAN_Cluster"]
    # - vizualizați și salvați plotul dbscan_<handle>.png
   db = DBSCAN(eps=1.5, min_samples=5)
   db_labels = db.fit_predict(X_scaled)
   df["DBSCAN_Cluster"] = db_labels

   plt.figure(figsize=(6, 5))
   plt.scatter(X_pca[:, 0], X_pca[:, 1], c=db_labels, s=20)
   plt.title("DBSCAN — PCA 2D")
   plt.xlabel("PC1")
   plt.ylabel("PC2")
   plt.tight_layout()
   plt.savefig(output_dir / f"dbscan_{handle}.png", dpi=150)
   plt.close()

    # TODO 7: Salvare rezultate
    # salvați un CSV cu coloanele ["Diagnosis", "KMeans_Cluster", "DBSCAN_Cluster"]
    # în clusters_<handle>.csv
   df[["Diagnosis", "KMeans_Cluster", "DBSCAN_Cluster"]].to_csv(
      output_dir / f"clusters_{handle}.csv", index=False
   )
