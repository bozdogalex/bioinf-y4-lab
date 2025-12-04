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
import numpy as np

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

    # Director pentru rezultate
    output_dir = Path(".")
    output_dir.mkdir(parents=True, exist_ok=True)

    # TODO 4: Hierarchical Clustering
    # - folosiți linkage(X_scaled, method="average")
    # - vizualizați cu dendrogram()
    # - salvați imaginea ca hierarchical_<handle>.png
    
    print("[...] Aplicare Hierarchical Clustering...")
    linkage_matrix = linkage(X_scaled, method="average")
    
    plt.figure(figsize=(12, 6))
    dendrogram(linkage_matrix, truncate_mode='lastp', p=30)
    plt.title("Dendrogramă - Hierarchical Clustering (Breast Cancer Dataset)")
    plt.xlabel("Index eșantion (sau cluster)")
    plt.ylabel("Distanță")
    plt.tight_layout()
    plt.savefig(output_dir / "hierarchical_AlexTGoCreative.png", dpi=300)
    plt.close()
    print(f"[OK] Salvat: hierarchical_AlexTGoCreative.png")

    # TODO 5: K-means Clustering
    # - aplicați KMeans cu K=2
    # - adăugați etichetele în df["KMeans_Cluster"]
    # - reduceți dimensionalitatea cu PCA(n_components=2)
    # - vizualizați și salvați plotul kmeans_<handle>.png
    
    print("[...] Aplicare K-means Clustering...")
    kmeans = KMeans(n_clusters=2, random_state=42, n_init=10)
    df["KMeans_Cluster"] = kmeans.fit_predict(X_scaled)
    
    # Reducere dimensionalitate pentru vizualizare
    pca = PCA(n_components=2)
    X_pca = pca.fit_transform(X_scaled)
    
    plt.figure(figsize=(10, 6))
    scatter = plt.scatter(X_pca[:, 0], X_pca[:, 1], 
                         c=df["KMeans_Cluster"], 
                         cmap='viridis', 
                         alpha=0.6, 
                         edgecolors='k', 
                         s=50)
    plt.colorbar(scatter, label='Cluster')
    plt.title("K-means Clustering (K=2) - Vizualizare PCA")
    plt.xlabel(f"PC1 ({pca.explained_variance_ratio_[0]:.2%} varianță)")
    plt.ylabel(f"PC2 ({pca.explained_variance_ratio_[1]:.2%} varianță)")
    plt.tight_layout()
    plt.savefig(output_dir / "kmeans_AlexTGoCreative.png", dpi=300)
    plt.close()
    print(f"[OK] Salvat: kmeans_AlexTGoCreative.png")

    # TODO 6: DBSCAN Clustering
    # - aplicați DBSCAN (ex: eps=1.5, min_samples=5)
    # - adăugați etichetele în df["DBSCAN_Cluster"]
    # - vizualizați și salvați plotul dbscan_<handle>.png
    
    print("[...] Aplicare DBSCAN Clustering...")
    dbscan = DBSCAN(eps=1.5, min_samples=5)
    df["DBSCAN_Cluster"] = dbscan.fit_predict(X_scaled)
    
    n_clusters = len(set(df["DBSCAN_Cluster"])) - (1 if -1 in df["DBSCAN_Cluster"] else 0)
    n_noise = list(df["DBSCAN_Cluster"]).count(-1)
    
    plt.figure(figsize=(10, 6))
    scatter = plt.scatter(X_pca[:, 0], X_pca[:, 1], 
                         c=df["DBSCAN_Cluster"], 
                         cmap='plasma', 
                         alpha=0.6, 
                         edgecolors='k', 
                         s=50)
    plt.colorbar(scatter, label='Cluster (-1 = noise)')
    plt.title(f"DBSCAN Clustering - Vizualizare PCA\n{n_clusters} clustere, {n_noise} puncte zgomot")
    plt.xlabel(f"PC1 ({pca.explained_variance_ratio_[0]:.2%} varianță)")
    plt.ylabel(f"PC2 ({pca.explained_variance_ratio_[1]:.2%} varianță)")
    plt.tight_layout()
    plt.savefig(output_dir / "dbscan_AlexTGoCreative.png", dpi=300)
    plt.close()
    print(f"[OK] Salvat: dbscan_AlexTGoCreative.png")

    # TODO 7: Salvare rezultate
    # salvați un CSV cu coloanele ["Diagnosis", "KMeans_Cluster", "DBSCAN_Cluster"]
    # în clusters_<handle>.csv
    
    results_df = df[["Diagnosis", "KMeans_Cluster", "DBSCAN_Cluster"]]
    results_df.to_csv(output_dir / "clusters_AlexTGoCreative.csv", index=False)
    print(f"[OK] Salvat: clusters_AlexTGoCreative.csv")
    
    # Afișare statistici
    print("\n=== Statistici Clustering ===")
    print(f"Total eșantioane: {len(df)}")
    print(f"\nDiagnostic real:")
    print(df["Diagnosis"].value_counts())
    print(f"\nK-means clustere:")
    print(df["KMeans_Cluster"].value_counts())
    print(f"\nDBSCAN clustere:")
    print(df["DBSCAN_Cluster"].value_counts())
    
    # Comparare cu diagnosticul real
    from sklearn.metrics import adjusted_rand_score
    kmeans_ari = adjusted_rand_score(df["Diagnosis"], df["KMeans_Cluster"])
    dbscan_ari = adjusted_rand_score(df["Diagnosis"], df["DBSCAN_Cluster"])
    
    print(f"\n=== Adjusted Rand Index (vs Diagnosis) ===")
    print(f"K-means ARI: {kmeans_ari:.3f}")
    print(f"DBSCAN ARI: {dbscan_ari:.3f}")
   
