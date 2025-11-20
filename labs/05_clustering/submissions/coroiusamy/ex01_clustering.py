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

# Importurile necesare pentru clustering și vizualizare
from scipy.cluster.hierarchy import dendrogram, linkage
from sklearn.cluster import KMeans, DBSCAN
from sklearn.decomposition import PCA

if __name__ == "__main__":
    
    github_handle = "coroiusamy" 

    # TODO 1: Încărcați dataset-ul
    url = "https://archive.ics.uci.edu/ml/machine-learning-databases/breast-cancer-wisconsin/wdbc.data"

    columns = ["ID", "Diagnosis"] + [f"Feature_{i}" for i in range(1, 31)]
    df = pd.read_csv(url, header=None, names=columns)

    # TODO 2: Preprocesare
    # - eliminați coloana ID
    df = df.drop(columns=["ID"])
    # diagnosticul (M=1, B=0)
    df["Diagnosis"] = df["Diagnosis"].apply(lambda x: 1 if x == "M" else 0)

    # TODO 3: Standardizare
    # Separăm features (X) de eticheta reală (Diagnosis)
    X = df.drop(columns=["Diagnosis"])
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    # Director pentru rezultate
    output_dir = Path(f"labs/05_clustering/submissions/{github_handle}")
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"Datele au fost încărcate și standardizate. Se salvează rezultatele în: {output_dir}")

    # --- TODO 4: Hierarchical Clustering ---
    print("Se rulează Hierarchical Clustering...")
    Z = linkage(X_scaled, method="average")

    plt.figure(figsize=(15, 7))
    dendrogram(Z, leaf_rotation=90., leaf_font_size=8.,
               truncate_mode='lastp', p=30)
    plt.title(f"Dendrogramă Hierarchical Clustering (handle: {github_handle})")
    plt.ylabel("Distanță (Average Linkage)")
    plt.xlabel("Index eșantion sau mărime cluster (dacă e trunchiat)")

    plt.savefig(output_dir / f"hierarchical_{github_handle}.png")
    plt.close() # Închidem figura pentru a elibera memoria

    # --- TODO 5: K-means Clustering ---
    print("Se rulează K-means (K=2)...")
    # KMeans cu K=2
    kmeans = KMeans(n_clusters=2, random_state=42, n_init=10) 
    kmeans_labels = kmeans.fit_predict(X_scaled)
    
    df["KMeans_Cluster"] = kmeans_labels

    # PCA(n_components=2)
    pca = PCA(n_components=2)
    X_pca = pca.fit_transform(X_scaled)

    # plotul kmeans_<handle>.png
    plt.figure(figsize=(8, 6))
    plt.scatter(X_pca[:, 0], X_pca[:, 1], c=kmeans_labels, cmap="viridis", s=50, alpha=0.7)
    plt.title(f"K-means Clustering (K=2) pe date PCA (handle: {github_handle})")
    plt.xlabel("Componenta Principală 1")
    plt.ylabel("Componenta Principală 2")
    plt.colorbar(label="Cluster")
    plt.savefig(output_dir / f"kmeans_{github_handle}.png")
    plt.close()

    # --- TODO 6: DBSCAN Clustering ---
    print("Se rulează DBSCAN...")
    #  DBSCAN (eps=1.5, min_samples=5)
    dbscan = DBSCAN(eps=1.5, min_samples=5)
    dbscan_labels = dbscan.fit_predict(X_scaled)

    # etichetele în df["DBSCAN_Cluster"]
    df["DBSCAN_Cluster"] = dbscan_labels

    # Reutilizăm X_pca calculat anterior
    plt.figure(figsize=(8, 6))
    plt.scatter(X_pca[:, 0], X_pca[:, 1], c=dbscan_labels, cmap="plasma", s=50, alpha=0.7)
    plt.title(f"DBSCAN Clustering (eps=1.5, min=5) pe date PCA (handle: {github_handle})")
    plt.xlabel("Componenta Principală 1")
    plt.ylabel("Componenta Principală 2")
    plt.colorbar(label="Cluster ( -1 = Zgomot )")
    plt.savefig(output_dir / f"dbscan_{github_handle}.png")
    plt.close()

    # --- TODO 7: Salvare rezultate ---
    print("Se salvează etichetele clusterelor în CSV...")
    # - CSV cu coloanele ["Diagnosis", "KMeans_Cluster", "DBSCAN_Cluster"]
    output_csv_path = output_dir / f"clusters_{github_handle}.csv"
    df[["Diagnosis", "KMeans_Cluster", "DBSCAN_Cluster"]].to_csv(output_csv_path, index=False)

    print(f"\n--- Script finalizat ---")
    print(f"Toate fișierele au fost salvate în: {output_dir.resolve()}")