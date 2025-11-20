
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans, DBSCAN
from sklearn.decomposition import PCA
from scipy.cluster.hierarchy import dendrogram, linkage
from pathlib import Path


def load_breast_cancer_data():
    """
    Încarcă dataset-ul WDBC de pe UCI Repository.
    """
    print("Loading Breast Cancer Wisconsin Diagnostic dataset...")
    url = "https://archive.ics.uci.edu/ml/machine-learning-databases/breast-cancer-wisconsin/wdbc.data"
    columns = ["ID", "Diagnosis"] + [f"Feature_{i}" for i in range(1, 31)]
    
    df = pd.read_csv(url, header=None, names=columns)
    print(f"Dataset loaded: {df.shape[0]} samples, {df.shape[1]} columns")
    
    return df


def preprocess_data(df):
    """
    Preprocesează datele: elimină ID, transformă Diagnosis în numeric.
    """
    print("\nPreprocessing data...")
    
    # Elimină coloana ID (irelevantă pentru clustering)
    df = df.drop(columns=["ID"])
    
    # Transformă Diagnosis: M (Malignant) = 1, B (Benign) = 0
    df["Diagnosis"] = df["Diagnosis"].apply(lambda x: 1 if x == "M" else 0)
    
    print(f"Diagnosis distribution:")
    print(f"  Malignant (1): {(df['Diagnosis'] == 1).sum()}")
    print(f"  Benign (0): {(df['Diagnosis'] == 0).sum()}")
    
    return df


def standardize_features(df):
    """
    Standardizează features (mean=0, std=1).
    """
    print("\nStandardizing features...")
    X = df.drop(columns=["Diagnosis"])
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    
    print(f"Features standardized: {X_scaled.shape}")
    
    return X_scaled, df


def hierarchical_clustering(X_scaled, output_dir, handle):
    """
    Aplică Hierarchical Clustering și salvează dendrograma.
    """
    print("\n" + "="*70)
    print("HIERARCHICAL CLUSTERING")
    print("="*70)
    
    # Calculare linkage matrix (average linkage)
    print("Computing linkage matrix (average method)...")
    linkage_matrix = linkage(X_scaled, method='average')
    
    # Vizualizare dendrogramă
    plt.figure(figsize=(12, 8))
    dendrogram(linkage_matrix, 
               truncate_mode='lastp',  # Afișează doar ultimele p clustere
               p=30,  # Număr de clustere afișate
               show_leaf_counts=True,
               leaf_font_size=10)
    
    plt.title('Hierarchical Clustering Dendrogram (Average Linkage)', fontsize=14)
    plt.xlabel('Sample Index (or Cluster Size)', fontsize=12)
    plt.ylabel('Distance', fontsize=12)
    plt.tight_layout()
    
    # Salvare
    output_file = output_dir / f"hierarchical_{handle}.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"✓ Dendrogram saved to: {output_file}")


def kmeans_clustering(X_scaled, df, output_dir, handle):
    """
    Aplică K-means Clustering (K=2) și vizualizează cu PCA.
    """
    print("\n" + "="*70)
    print("K-MEANS CLUSTERING (K=2)")
    print("="*70)
    
    # K-means
    print("Running K-means with K=2...")
    kmeans = KMeans(n_clusters=2, random_state=42, n_init=10)
    kmeans_labels = kmeans.fit_predict(X_scaled)
    
    df["KMeans_Cluster"] = kmeans_labels
    
    print(f"Cluster distribution:")
    print(f"  Cluster 0: {(kmeans_labels == 0).sum()}")
    print(f"  Cluster 1: {(kmeans_labels == 1).sum()}")
    
    # PCA pentru vizualizare
    print("Applying PCA for 2D visualization...")
    pca = PCA(n_components=2)
    X_pca = pca.fit_transform(X_scaled)
    
    print(f"Explained variance: {pca.explained_variance_ratio_}")
    
    # Vizualizare
    plt.figure(figsize=(10, 8))
    scatter = plt.scatter(X_pca[:, 0], X_pca[:, 1], 
                         c=kmeans_labels, 
                         cmap='viridis', 
                         s=50, 
                         alpha=0.7,
                         edgecolors='k',
                         linewidth=0.5)
    
    plt.title('K-means Clustering (K=2) - PCA Visualization', fontsize=14)
    plt.xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.2%} variance)', fontsize=12)
    plt.ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.2%} variance)', fontsize=12)
    plt.colorbar(scatter, label='Cluster')
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    
    # Salvare
    output_file = output_dir / f"kmeans_{handle}.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"✓ K-means plot saved to: {output_file}")
    
    return X_pca


def dbscan_clustering(X_scaled, df, X_pca, output_dir, handle):
    """
    Aplică DBSCAN Clustering și vizualizează cu PCA.
    """
    print("\n" + "="*70)
    print("DBSCAN CLUSTERING")
    print("="*70)
    
    # DBSCAN
    print("Running DBSCAN (eps=1.5, min_samples=5)...")
    dbscan = DBSCAN(eps=1.5, min_samples=5)
    dbscan_labels = dbscan.fit_predict(X_scaled)
    
    df["DBSCAN_Cluster"] = dbscan_labels
    
    # Statistici
    n_clusters = len(set(dbscan_labels)) - (1 if -1 in dbscan_labels else 0)
    n_noise = list(dbscan_labels).count(-1)
    
    print(f"Number of clusters: {n_clusters}")
    print(f"Number of noise points: {n_noise}")
    print(f"Cluster distribution:")
    for label in set(dbscan_labels):
        count = (dbscan_labels == label).sum()
        if label == -1:
            print(f"  Noise (-1): {count}")
        else:
            print(f"  Cluster {label}: {count}")
    
    # Vizualizare
    plt.figure(figsize=(10, 8))
    
    # Colorare specială pentru noise (-1)
    unique_labels = set(dbscan_labels)
    colors = plt.cm.Spectral(np.linspace(0, 1, len(unique_labels)))
    
    for k, col in zip(unique_labels, colors):
        if k == -1:
            # Noise - negru
            col = 'k'
            marker = 'x'
        else:
            marker = 'o'
        
        class_member_mask = (dbscan_labels == k)
        xy = X_pca[class_member_mask]
        plt.scatter(xy[:, 0], xy[:, 1], 
                   c=[col], 
                   marker=marker,
                   s=50 if k != -1 else 30,
                   alpha=0.7,
                   edgecolors='k' if k != -1 else None,
                   linewidth=0.5,
                   label=f'Cluster {k}' if k != -1 else 'Noise')
    
    plt.title(f'DBSCAN Clustering - PCA Visualization\n({n_clusters} clusters, {n_noise} noise points)', 
              fontsize=14)
    plt.xlabel('PC1', fontsize=12)
    plt.ylabel('PC2', fontsize=12)
    plt.legend(loc='best')
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    
    # Salvare
    output_file = output_dir / f"dbscan_{handle}.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"✓ DBSCAN plot saved to: {output_file}")


def save_results(df, output_dir, handle):
    """
    Salvează rezultatele clustering în CSV.
    """
    print("\n" + "="*70)
    print("SAVING RESULTS")
    print("="*70)
    
    # Selectare coloane relevante
    results_df = df[["Diagnosis", "KMeans_Cluster", "DBSCAN_Cluster"]]
    
    # Salvare CSV
    output_file = output_dir / f"clusters_{handle}.csv"
    results_df.to_csv(output_file, index=False)
    
    print(f"✓ Cluster assignments saved to: {output_file}")
    print(f"\nFirst 10 rows:")
    print(results_df.head(10))


def main():
    handle = "IrisDanila"
    
    print("="*70)
    print("LAB 5 - CLUSTERING ANALYSIS ON BREAST CANCER DATA")
    print("="*70)
    
    # Director pentru rezultate
    output_dir = Path(f"labs/05_clustering/submissions/{handle}")
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # 1. Încărcare date
    df = load_breast_cancer_data()
    
    # 2. Preprocesare
    df = preprocess_data(df)
    
    # 3. Standardizare
    X_scaled, df = standardize_features(df)
    
    # 4. Hierarchical Clustering
    hierarchical_clustering(X_scaled, output_dir, handle)
    
    # 5. K-means Clustering
    X_pca = kmeans_clustering(X_scaled, df, output_dir, handle)
    
    # 6. DBSCAN Clustering
    dbscan_clustering(X_scaled, df, X_pca, output_dir, handle)
    
    # 7. Salvare rezultate
    save_results(df, output_dir, handle)
    
    print("\n" + "="*70)
    print("CLUSTERING ANALYSIS COMPLETE!")
    print("="*70)
    print(f"\nGenerated files in {output_dir}:")
    print(f"  - clusters_{handle}.csv")
    print(f"  - hierarchical_{handle}.png")
    print(f"  - kmeans_{handle}.png")
    print(f"  - dbscan_{handle}.png")
    print("="*70)


if __name__ == "__main__":
    main()
