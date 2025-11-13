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

# TODO: Adăugați importurile necesare
from scipy.cluster.hierarchy import dendrogram, linkage
from sklearn.cluster import KMeans, DBSCAN
from sklearn.decomposition import PCA

if __name__ == "__main__":
    # Definiți handle-ul
    handle = "RAZVAN24RR"

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

    # Director pentru rezultate
    output_dir = Path(f"labs/05_clustering/submissions/{handle}")
    output_dir.mkdir(parents=True, exist_ok=True)

    # TODO 4: Hierarchical Clustering
    # - folosiți linkage(X_scaled, method="average")
    # - vizualizați cu dendrogram()
    # - salvați imaginea ca hierarchical_<handle>.png
    print("Se rulează Hierarchical Clustering...")
    plt.figure(figsize=(12, 8))
    Z = linkage(X_scaled, method="average")
    dendrogram(Z, leaf_rotation=90., leaf_font_size=8.,
               labels=df.index,
               truncate_mode='level', p=5) # Trunchiem pentru vizibilitate
    plt.title(f'Dendrogramă Hierarchical Clustering - {handle}')
    plt.xlabel('Index Eșantion')
    plt.ylabel('Distanță')
    plt.tight_layout()
    plt.savefig(output_dir / f"hierarchical_{handle}.png")
    plt.close() # Închidem figura pentru a elibera memoria

    # TODO 5: K-means Clustering
    # - aplicați KMeans cu K=2
    # - adăugați etichetele în df["KMeans_Cluster"]
    # - reduceți dimensionalitatea cu PCA(n_components=2)
    # - vizualizați și salvați plotul kmeans_<handle>.png
    print("Se rulează K-means Clustering...")
    kmeans = KMeans(n_clusters=2, random_state=42, n_init=10)
    kmeans_labels = kmeans.fit_predict(X_scaled)
    df["KMeans_Cluster"] = kmeans_labels

    # Aplicăm PCA pentru vizualizare (se putea face o singură dată)
    pca = PCA(n_components=2)
    X_pca = pca.fit_transform(X_scaled)

    plt.figure(figsize=(8, 6))
    scatter = plt.scatter(X_pca[:, 0], X_pca[:, 1], c=kmeans_labels, cmap='viridis', alpha=0.7)
    plt.title(f'K-means Clustering (K=2) - Vizualizare PCA - {handle}')
    plt.xlabel('Componenta Principală 1')
    plt.ylabel('Componenta Principală 2')
    plt.legend(handles=scatter.legend_elements()[0], labels=['Cluster 0', 'Cluster 1'])
    plt.savefig(output_dir / f"kmeans_{handle}.png")
    plt.close()

    # TODO 6: DBSCAN Clustering
    # - aplicați DBSCAN (ex: eps=1.5, min_samples=5)
    # - adăugați etichetele în df["DBSCAN_Cluster"]
    # - vizualizați și salvați plotul dbscan_<handle>.png
    print("Se rulează DBSCAN Clustering...")
    # Parametrii pot necesita ajustare, folosim cei sugerați
    dbscan = DBSCAN(eps=1.5, min_samples=5) 
    dbscan_labels = dbscan.fit_predict(X_scaled)
    df["DBSCAN_Cluster"] = dbscan_labels

    plt.figure(figsize=(8, 6))
    plt.scatter(X_pca[:, 0], X_pca[:, 1], c=dbscan_labels, cmap='plasma', alpha=0.7)
    plt.title(f'DBSCAN Clustering - Vizualizare PCA - {handle}')
    plt.xlabel('Componenta Principală 1')
    plt.ylabel('Componenta Principală 2')
    # Adăugăm o bară de culoare pentru a vedea etichetele (inclusiv zgomotul -1)
    plt.colorbar(label='Etichetă Cluster DBSCAN')
    plt.savefig(output_dir / f"dbscan_{handle}.png")
    plt.close()

    # TODO 7: Salvare rezultate
    # salvați un CSV cu coloanele ["Diagnosis", "KMeans_Cluster", "DBSCAN_Cluster"]
    # în clusters_<handle>.csv
    print("Se salvează rezultatele...")
    output_csv_path = output_dir / f"clusters_{handle}.csv"
    df[["Diagnosis", "KMeans_Cluster", "DBSCAN_Cluster"]].to_csv(output_csv_path, index=False)

    print(f"Execuție finalizată. Rezultatele au fost salvate în {output_dir}")