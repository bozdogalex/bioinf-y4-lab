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
import numpy as np
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans, DBSCAN
from sklearn.decomposition import PCA
from scipy.cluster.hierarchy import dendrogram, linkage
from pathlib import Path

# Configurare pentru afișarea mai frumoasă a ploturilor
plt.style.use('default')
plt.rcParams['figure.figsize'] = (10, 6)

if __name__ == "__main__":
    # TODO 1: Încărcați dataset-ul
    print("Încărcare dataset...")
    url = "https://archive.ics.uci.edu/ml/machine-learning-databases/breast-cancer-wisconsin/wdbc.data"
    column_names = ['ID', 'Diagnosis'] + [f'Feature_{i}' for i in range(1, 31)]
    df = pd.read_csv(url, header=None, names=column_names)
    
    print(f"Dimensiune dataset: {df.shape}")
    print(f"Primele 5 rânduri:\n{df.head()}")
    print(f"\nDistribuția diagnosticelor:\n{df['Diagnosis'].value_counts()}")

    # TODO 2: Preprocesare
    # - eliminați coloana ID
    # - transformați diagnosticul în valori numerice (M=1, B=0)
    print("\nPreprocesare date...")
    df_processed = df.drop(columns=['ID'])
    df_processed['Diagnosis_Numeric'] = df_processed['Diagnosis'].map({'M': 1, 'B': 0})
    
    # Separă features de target
    X = df_processed.drop(columns=['Diagnosis', 'Diagnosis_Numeric'])
    y = df_processed['Diagnosis_Numeric']
    
    print(f"Features shape: {X.shape}")
    print(f"Target shape: {y.shape}")

    # TODO 3: Standardizare
    print("\nStandardizare date...")
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    print(f"Date standardizate shape: {X_scaled.shape}")

    # Director pentru rezultate
    handle = "Botoaca-Florentina-Veronica"  
    output_dir = Path(f"labs/05_clustering/submissions/{handle}")
    output_dir.mkdir(parents=True, exist_ok=True)
    print(f"\nDirector output: {output_dir}")

    # TODO 4: Hierarchical Clustering
    print("\nAplicare Hierarchical Clustering...")
    # Folosim un subset pentru dendrogramă mai clară (primile 30 de probe)
    n_samples_dendro = min(30, X_scaled.shape[0])
    Z = linkage(X_scaled[:n_samples_dendro], method="average", metric='euclidean')
    
    plt.figure(figsize=(12, 8))
    dendrogram(Z, 
               labels=df['Diagnosis'].values[:n_samples_dendro],
               leaf_rotation=90,
               leaf_font_size=10)
    plt.title(f'Hierarchical Clustering Dendrogram\n(Primile {n_samples_dendro} probe) - {handle}')
    plt.xlabel('Probe (etichete: Diagnosis)')
    plt.ylabel('Distanță')
    plt.tight_layout()
    plt.savefig(output_dir / f'hierarchical_{handle}.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("Dendrogramă salvată!")

    # TODO 5: K-means Clustering
    print("\nAplicare K-means Clustering...")
    kmeans = KMeans(n_clusters=2, random_state=42, n_init=10)
    kmeans_labels = kmeans.fit_predict(X_scaled)
    df_processed['KMeans_Cluster'] = kmeans_labels
    
    # Vizualizare cu PCA
    pca = PCA(n_components=2)
    X_pca = pca.fit_transform(X_scaled)
    
    plt.figure(figsize=(10, 8))
    scatter = plt.scatter(X_pca[:, 0], X_pca[:, 1], 
                         c=kmeans_labels, cmap='viridis', alpha=0.7)
    plt.colorbar(scatter, label='Cluster')
    plt.title(f'K-means Clustering (K=2) - PCA Vizualizare - {handle}')
    plt.xlabel(f'Componenta Principală 1 ({pca.explained_variance_ratio_[0]:.2%} varianță)')
    plt.ylabel(f'Componenta Principală 2 ({pca.explained_variance_ratio_[1]:.2%} varianță)')
    
    # Adăugăm centroizii proiectați în spațiul PCA
    centers_pca = pca.transform(kmeans.cluster_centers_)
    plt.scatter(centers_pca[:, 0], centers_pca[:, 1], 
               c='red', marker='X', s=200, label='Centroizi')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(output_dir / f'kmeans_{handle}.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("Vizualizare K-means salvată!")

    # TODO 6: DBSCAN Clustering
    print("\nAplicare DBSCAN Clustering...")
    dbscan = DBSCAN(eps=1.5, min_samples=5)
    dbscan_labels = dbscan.fit_predict(X_scaled)
    df_processed['DBSCAN_Cluster'] = dbscan_labels
    
    # Vizualizare DBSCAN cu PCA
    plt.figure(figsize=(10, 8))
    scatter = plt.scatter(X_pca[:, 0], X_pca[:, 1], 
                         c=dbscan_labels, cmap='Set3', alpha=0.7)
    plt.colorbar(scatter, label='Cluster')
    plt.title(f'DBSCAN Clustering - PCA Vizualizare - {handle}')
    plt.xlabel(f'Componenta Principală 1 ({pca.explained_variance_ratio_[0]:.2%} varianță)')
    plt.ylabel(f'Componenta Principală 2 ({pca.explained_variance_ratio_[1]:.2%} varianță)')
    
    # Highlight outliers (cluster = -1)
    outliers_mask = dbscan_labels == -1
    if np.any(outliers_mask):
        plt.scatter(X_pca[outliers_mask, 0], X_pca[outliers_mask, 1], 
                   c='red', marker='x', s=100, label='Outliers')
        plt.legend()
    
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(output_dir / f'dbscan_{handle}.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("Vizualizare DBSCAN salvată!")

    # Analiză comparativă
    print("\n=== ANALIZĂ COMPARATIVĂ ===")
    print(f"K-means clusters: {np.unique(kmeans_labels)}")
    print(f"DBSCAN clusters: {np.unique(dbscan_labels)}")
    print(f"Număr de outliers DBSCAN: {np.sum(dbscan_labels == -1)}")
    
    # Comparație cu diagnosticul real
    from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score
    
    ari_kmeans = adjusted_rand_score(y, kmeans_labels)
    ari_dbscan = adjusted_rand_score(y, dbscan_labels) if len(np.unique(dbscan_labels)) > 1 else 0
    
    print(f"\nAdjusted Rand Score vs diagnostic real:")
    print(f"K-means: {ari_kmeans:.3f}")
    print(f"DBSCAN: {ari_dbscan:.3f}")

    # TODO 7: Salvare rezultate
    print("\nSalvare rezultate...")
    results_df = pd.DataFrame({
        'Diagnosis': df_processed['Diagnosis'],
        'Diagnosis_Numeric': df_processed['Diagnosis_Numeric'],
        'KMeans_Cluster': df_processed['KMeans_Cluster'],
        'DBSCAN_Cluster': df_processed['DBSCAN_Cluster']
    })
    
    results_df.to_csv(output_dir / f'clusters_{handle}.csv', index=False)
    
    # Salvare statistici
    stats_df = pd.DataFrame({
        'Metrică': ['K-means ARI', 'DBSCAN ARI', 'K-means Clusters', 'DBSCAN Clusters', 'Outliers DBSCAN'],
        'Valoare': [ari_kmeans, ari_dbscan, len(np.unique(kmeans_labels)), 
                   len(np.unique(dbscan_labels)), np.sum(dbscan_labels == -1)]
    })
    stats_df.to_csv(output_dir / f'clustering_stats_{handle}.csv', index=False)
    
    print("Rezultate salvate!")
    print(f"\nFișiere generate în: {output_dir}")
    print(f"✓ hierarchical_{handle}.png")
    print(f"✓ kmeans_{handle}.png") 
    print(f"✓ dbscan_{handle}.png")
    print(f"✓ clusters_{handle}.csv")
    print(f"✓ clustering_stats_{handle}.csv")
    
    # Afișare summary
    print(f"\n=== SUMMARY ===")
    print(f"Total probe: {len(results_df)}")
    print(f"Diagnostic M (Malign): {sum(results_df['Diagnosis'] == 'M')}")
    print(f"Diagnostic B (Benign): {sum(results_df['Diagnosis'] == 'B')}")
    print(f"K-means vs real ARI: {ari_kmeans:.3f}")