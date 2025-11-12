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
from scipy.cluster.hierarchy import dendrogram, linkage
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from pathlib import Path
from sklearn.cluster import KMeans, DBSCAN
from sklearn.decomposition import PCA

# Opțional: puteți importa deja funcțiile necesare
# from scipy.cluster.hierarchy import dendrogram, linkage
# from sklearn.cluster import KMeans, DBSCAN
# 

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
    output_dir = Path("labs/05_clustering/submissions/MariusJalba")
    output_dir.mkdir(parents=True, exist_ok=True)

    # TODO 4: Hierarchical Clustering
    Z = linkage(X_scaled, method="average")
    plt.figure(figsize=(50,30))
    dendrogram(Z)
    out_path = "labs/05_clustering/submissions/MariusJalba/hierarchical_MariusJalba.png"
    plt.tight_layout()
    plt.savefig(out_path, dpi=200)
    plt.close()

    # TODO 5: K-means Clustering
    # - aplicați KMeans cu K=2
    # - adăugați etichetele în df["KMeans_Cluster"]
    # - reduceți dimensionalitatea cu PCA(n_components=2)
    # - vizualizați și salvați plotul kmeans_<handle>.png
    kmeans = KMeans(n_clusters=2, init='k-means++', n_init='auto', max_iter=300, random_state=42)
    kmeans_labels = kmeans.fit_predict(X_scaled)
    df["KMeans_Cluster"] = kmeans_labels
    pca = PCA(n_components=2, random_state=42)
    X_pca = pca.fit_transform(X_scaled)
    df["PC1"] = X_pca[:,0]
    df["PC2"] = X_pca[:,1]
    plt.figure(figsize=(10, 8))
    plt.scatter(df["PC1"], df["PC2"], c=df["KMeans_Cluster"])
    plt.title("KMeans (K=2) on WDBC – PCA(2D)")
    plt.xlabel("PC1")
    plt.ylabel("PC2")
    plt.tight_layout()
    out_scatter =  "labs/05_clustering/submissions/MariusJalba/kmeans_MariusJalba.png"
    plt.savefig(out_scatter, dpi=200)
    plt.close()


    # TODO 6: DBSCAN Clustering
    # - aplicați DBSCAN (ex: eps=1.5, min_samples=5)
    # - adăugați etichetele în df["DBSCAN_Cluster"]
    # - vizualizați și salvați plotul dbscan_<handle>.png
    db = DBSCAN(eps=1.5, min_samples=5).fit(X_scaled)
    labels = db.labels_
    df["DBSCAN_Cluster"] = labels
    coords = PCA(n_components=2, random_state=42).fit_transform(X_scaled) #aplicam PCA sa putem reprezenta toate dimensiunile in doar 2 axe
    plt.scatter(coords[:,0], coords[:,1], c=df["DBSCAN_Cluster"], cmap="rainbow", s=10)
    plt.title("DBSCAN Clustering (PCA 2D)")
    plt.xlabel("PC1")
    plt.ylabel("PC2")
    plt.tight_layout()
    plt.savefig("labs/05_clustering/submissions/MariusJalba/dbscan_MariusJalba.png")
    plt.close()

    # TODO 7: Salvare rezultate
    # salvați un CSV cu coloanele ["Diagnosis", "KMeans_Cluster", "DBSCAN_Cluster"]
    # în clusters_<handle>.csv
    path = "labs/05_clustering/submissions/MariusJalba/clusters_MariusJalba.csv"
    df[["Diagnosis", "KMeans_Cluster","DBSCAN_Cluster"]].to_csv(path, index=False)

    print("Toate fisierele de output au fost salvate in submissions/MariusJalba") 
    
