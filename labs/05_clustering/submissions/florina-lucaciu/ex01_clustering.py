import pandas as pd
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from pathlib import Path

from scipy.cluster.hierarchy import dendrogram, linkage
from sklearn.cluster import KMeans, DBSCAN
from sklearn.decomposition import PCA

if __name__ == "__main__":
    # schimba handle-ul daca e nevoie
    handle = "florina-lucaciu" #

    url = "https://archive.ics.uci.edu/ml/machine-learning-databases/breast-cancer-wisconsin/wdbc.data"
    columns = ["ID", "Diagnosis"] + [f"Feature_{i}" for i in range(1, 31)]
    df = pd.read_csv(url, header=None, names=columns)

    # Stergem coloana ID
    df = df.drop(columns=["ID"])
    # Convertim Diagnosis: Malign (M) -> 1, Benign (B) -> 0
    df["Diagnosis"] = df["Diagnosis"].apply(lambda x: 1 if x == "M" else 0)

    # Separam datele (X) de etichete (Diagnosis)
    X = df.drop(columns=["Diagnosis"])
    scaler = StandardScaler()
    # Scalam datele
    X_scaled = scaler.fit_transform(X)

    # Directorul unde salvam fisierele
    output_dir = Path(f"labs/05_clustering/submissions/{handle}")
    output_dir.mkdir(parents=True, exist_ok=True)
    print(f"Salvam rezultatele in: {output_dir.resolve()}")

    print("Se ruleaza Hierarchical Clustering...")
    # Calculam legaturile (linkage) intre puncte
    Z = linkage(X_scaled, method="average")

    # Cream dendrograma
    plt.figure(figsize=(12, 8))
    dendrogram(Z)
    plt.title(f'Dendrograma Hierarchical Clustering ({handle})')
    plt.xlabel('Index esantion')
    plt.ylabel('Distanta (Average Linkage)')
    
    # Salvam imaginea
    plt.savefig(output_dir / f"hierarchical_{handle}.png")
    print(f"Salvat: hierarchical_{handle}.png")
    plt.close() # Inchidem plotul

    print("Se ruleaza K-Means...")
    # Aplicam K-Means cu 2 clustere
    kmeans = KMeans(n_clusters=2, random_state=42, n_init=10)
    kmeans_labels = kmeans.fit_predict(X_scaled)
    
    # Adaugam clusterele K-Means in tabel
    df["KMeans_Cluster"] = kmeans_labels

    # Reducem datele la 2 dimensiuni (PCA) pentru vizualizare
    pca = PCA(n_components=2)
    X_pca = pca.fit_transform(X_scaled)

    # Vizualizam clusterele K-Means
    plt.figure(figsize=(8, 6))
    scatter = plt.scatter(X_pca[:, 0], X_pca[:, 1], c=kmeans_labels, cmap='viridis', alpha=0.7)
    plt.title(f'K-Means Clustering (K=2) - PCA View ({handle})')
    plt.xlabel('Componenta Principala 1')
    plt.ylabel('Componenta Principala 2')
    plt.legend(handles=scatter.legend_elements()[0], labels=['Cluster 0', 'Cluster 1'])
    
    # Salvam imaginea
    plt.savefig(output_dir / f"kmeans_{handle}.png")
    print(f"Salvat: kmeans_{handle}.png")
    plt.close()

    print("Se ruleaza DBSCAN...")
    # Aplicam DBSCAN
    # (eps=5 si min_samples=10 sunt alese pentru aceste date scalate)
    dbscan = DBSCAN(eps=5, min_samples=10) 
    dbscan_labels = dbscan.fit_predict(X_scaled)

    # Adaugam clusterele DBSCAN in tabel
    df["DBSCAN_Cluster"] = dbscan_labels

    # Vizualizam clusterele DBSCAN (folosim acelasi X_pca)
    plt.figure(figsize=(8, 6))
    scatter_db = plt.scatter(X_pca[:, 0], X_pca[:, 1], c=dbscan_labels, cmap='plasma', alpha=0.7)
    
    # Cream legenda (eticheta -1 inseamna "zgomot" sau "outlier")
    unique_labels = sorted(list(set(dbscan_labels)))
    legend_labels = [f'Cluster {l}' if l != -1 else 'Zgomot (Noise)' for l in unique_labels]
    
    plt.title(f'DBSCAN Clustering - PCA View ({handle})')
    plt.xlabel('Componenta Principala 1')
    plt.ylabel('Componenta Principala 2')
    plt.legend(handles=scatter_db.legend_elements()[0], labels=legend_labels)

    # Salvam imaginea
    plt.savefig(output_dir / f"dbscan_{handle}.png")
    print(f"Salvat: dbscan_{handle}.png")
    plt.close()

    # Salvam un CSV doar cu diagnoza reala si clusterele gasite
    output_csv_path = output_dir / f"clusters_{handle}.csv"
    df_results = df[["Diagnosis", "KMeans_Cluster", "DBSCAN_Cluster"]]
    df_results.to_csv(output_csv_path, index=False)
    print(f"Salvat: clusters_{handle}.csv")
    
    # comparatie rapida in terminal
    print("\n--- Analiza rapida a clusterelor ---")
    print("\nK-Means vs Diagnoza Reala:")
    print(pd.crosstab(df['Diagnosis'], df['KMeans_Cluster']))
    
    print("\nDBSCAN vs Diagnoza Reala:")
    print(pd.crosstab(df['Diagnosis'], df['DBSCAN_Cluster']))
    print("(Nota: -1 inseamna zgomot/outlier in DBSCAN)")
    
    print("\nFinalizat!")