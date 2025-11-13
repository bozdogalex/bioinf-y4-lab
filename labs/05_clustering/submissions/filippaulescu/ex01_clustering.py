"""
Exercițiu 6 — Clustering pe date de cancer mamar (toy dataset)

Implementeaza clustering-ul pe dataset-ul WDBC folosind metodele:
- Hierarchical Clustering (Dendrograma)
- K-Means (K=2)
- DBSCAN
Vizualizarea se face folosind PCA 2D pentru K-Means si DBSCAN.
"""

import pandas as pd
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans, DBSCAN
from scipy.cluster.hierarchy import dendrogram, linkage
from pathlib import Path
import numpy as np

# --- Setări ---
HANDLE_NAME = "filippaulescu"
OUTPUT_DIR = Path(f"labs/05_clustering/submissions/{HANDLE_NAME}")
DATA_URL = "https://archive.ics.uci.edu/ml/machine-learning-databases/breast-cancer-wisconsin/wdbc.data"
COLUMN_NAMES = ["ID", "Diagnosis"] + [f"Feature_{i}" for i in range(1, 31)]


def load_and_preprocess_data(url, columns):
    """Încarcă dataset-ul și efectuează preprocesarea de bază."""
    print("1. Incarcare si Preprocesare date...")
    df = pd.read_csv(url, header=None, names=columns)
    
    # Elimină coloana ID
    df = df.drop(columns=["ID"])
    
    # Transformă Diagnosis: M (Malign) -> 1, B (Benign) -> 0
    # Folosim .map pentru o citire mai clara
    diagnosis_map = {"M": 1, "B": 0}
    df["Diagnosis"] = df["Diagnosis"].map(diagnosis_map)
    
    return df

def standardize_features(df):
    """Standardizeaza coloanele de feature-uri (coloana 'Diagnosis' este exclusa)."""
    print("2. Standardizare date...")
    # Selectam doar coloanele numerice de feature-uri
    feature_cols = [col for col in df.columns if col.startswith("Feature")]
    X = df[feature_cols].values
    
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    
    return X_scaled

def plot_and_save_dendrogram(X_scaled, output_path):
    """Implementeaza si salveaza vizualizarea Hierarchical Clustering (Dendrograma)."""
    print("3. Ruleaza Hierarchical Clustering (Dendrograma)...")
    # Calculul matricii de legatura (linkage matrix)
    # Metoda "average" (UPGMA) este o alegere frecventa
    Z = linkage(X_scaled, method="average", metric='euclidean') 
    
    plt.figure(figsize=(15, 6))
    dendrogram(
        Z, 
        no_labels=True, 
        color_threshold=None, # Afiseaza toate linkurile
        leaf_rotation=90.
    )
    plt.title(f"Hierarchical Clustering (Average Linkage) - {HANDLE_NAME}")
    plt.xlabel("Esantioane (Samples)")
    plt.ylabel("Distanta Euclidiana")
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.tight_layout()
    plt.savefig(output_path, dpi=180)
    plt.close()
    print(f"   -> Dendrograma salvata: {output_path}")
    
    return Z

def perform_and_plot_clustering(df, X_scaled, method_name, model, output_path, pca_components):
    """Aplica un model de clustering si vizualizeaza rezultatele folosind PCA."""
    print(f"4. Ruleaza {method_name}...")
    
    # 4a. Aplicare Clustering
    labels = model.fit_predict(X_scaled)
    df[f"{method_name}_Cluster"] = labels
    
    # 4b. Vizualizare (folosind PCA)
    plt.figure(figsize=(7, 6))
    
    # Daca PCA nu a fost inca facut, il facem acum.
    # Folosim aceleasi componente principale (PC1, PC2) pentru toate ploturile.
    if pca_components is None:
        pca = PCA(n_components=2, random_state=42)
        pca_components = pca.fit_transform(X_scaled)
        
    # Numarul de clustere determinate (excluzand zgomotul -1 de la DBSCAN)
    n_clusters = len(np.unique(labels[labels != -1]))
    
    # Setarea culorilor in functie de etichete (labels)
    plt.scatter(
        pca_components[:, 0], 
        pca_components[:, 1], 
        c=labels, 
        s=30, 
        cmap='viridis' if n_clusters > 1 else 'gray',
        alpha=0.8
    )
    
    # Titlu specific
    title = f"{method_name} Clustering (PCA 2D)"
    if method_name == "KMeans":
        title += f" (K={model.n_clusters})"
    elif method_name == "DBSCAN":
        title += f" (Clusters: {n_clusters})"
        
    plt.title(title)
    plt.xlabel("Principal Component 1 (PC1)")
    plt.ylabel("Principal Component 2 (PC2)")
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.tight_layout()
    plt.savefig(output_path, dpi=180)
    plt.close()
    print(f"   -> Plotul {method_name} salvat: {output_path}")
    
    return df, pca_components

if __name__ == "__main__":
    # Creare director de iesire
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    
    # Pasul 1 & 2: Incarcare si Preprocesare
    df_original = load_and_preprocess_data(DATA_URL, COLUMN_NAMES)
    df_working = df_original.copy()
    
    # Pasul 3: Standardizare
    X_scaled = standardize_features(df_working)
    
    # PCA se face o singura data pentru vizualizare K-Means si DBSCAN
    pca = PCA(n_components=2, random_state=42)
    X_pca = pca.fit_transform(X_scaled)
    
    # ----------------------------------------------------
    # 3. Hierarchical Clustering
    # ----------------------------------------------------
    dendrogram_path = OUTPUT_DIR / f"hierarchical_{HANDLE_NAME}.png"
    plot_and_save_dendrogram(X_scaled, dendrogram_path)
    
    # ----------------------------------------------------
    # 4. K-Means Clustering (K=2)
    # ----------------------------------------------------
    kmeans_model = KMeans(n_clusters=2, random_state=42, n_init="auto")
    kmeans_path = OUTPUT_DIR / f"kmeans_{HANDLE_NAME}.png"
    df_working, X_pca = perform_and_plot_clustering(
        df_working, 
        X_scaled, 
        "KMeans", 
        kmeans_model, 
        kmeans_path, 
        X_pca
    )
    
    # ----------------------------------------------------
    # 5. DBSCAN Clustering (Hyperparametri: eps=1.5, min_samples=5)
    # ----------------------------------------------------
    # Nota: Acesti hiperparametri sunt specifici pentru datele standardizate
    dbscan_model = DBSCAN(eps=1.5, min_samples=5) 
    dbscan_path = OUTPUT_DIR / f"dbscan_{HANDLE_NAME}.png"
    df_working, X_pca = perform_and_plot_clustering(
        df_working, 
        X_scaled, 
        "DBSCAN", 
        dbscan_model, 
        dbscan_path, 
        X_pca
    )
    
    # ----------------------------------------------------
    # 6. Salvare Rezultate Finale
    # ----------------------------------------------------
    output_cols = ["Diagnosis", "KMeans_Cluster", "DBSCAN_Cluster"]
    results_path = OUTPUT_DIR / f"clusters_{HANDLE_NAME}.csv"
    
    # Ne asiguram ca nu salvam coloanele 'Feature_i' inutile
    if all(col in df_working.columns for col in output_cols):
        df_working[output_cols].to_csv(results_path, index=False)
        print(f"\n5. Rezultatele clustering salvate: {results_path}")
    else:
        print("\nEroare: Coloanele de cluster nu au fost gasite in DataFrame pentru salvare.")

    print("\n Executie finalizata pentru filippaulescu.")