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

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from sklearn.cluster import DBSCAN, KMeans
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from scipy.cluster.hierarchy import dendrogram, linkage

HANDLE = "olezhka0809"

DATA_PATH = Path("/workspaces/bioinf-y4-lab/data/work/{HANDLE}/lab05/wdbc.data") 


def load_wdbc(path: Path) -> pd.DataFrame:
    """
    Încarcă dataset-ul WDBC de pe UCI.
    Structura clasică (wdbc.data) este:
    ID, Diagnosis(M/B), 30 de features numerice.
    """

    if not path.exists():
        raise FileNotFoundError(f"Nu găsesc fișierul cu date WDBC: {path.resolve()}")

    # UCI WDBC are de obicei 32 de coloane: 1 ID, 1 Diagnosis, 30 features
    colnames = [
        "ID",
        "Diagnosis",
    ]
    # 30 de features: le numim generic f1..f30 dacă nu avem header
    colnames += [f"F{i}" for i in range(1, 31)]

    # fișier fără header -> header=None
    df = pd.read_csv(path, header=None, names=colnames)
    return df


if __name__ == "__main__":
    # ======================
    # 1) Încărcare dataset
    # ======================
    df = load_wdbc(DATA_PATH)
    print(f"[INFO] Date încărcate: {df.shape[0]} rânduri, {df.shape[1]} coloane")

    # ======================
    # 2) Preprocesare
    # ======================
    # - eliminăm coloana ID (nu e feature biologic)
    df = df.drop(columns=["ID"])

    # - transformăm diagnosticul în numeric: M=1 (malign), B=0 (benign)
    #   (dacă în fișier literele sunt mici, poți pune .str.upper() înainte)
    df["Diagnosis"] = df["Diagnosis"].map({"M": 1, "B": 0})

    # Verificăm că nu avem valori lipsă la diagnostic
    if df["Diagnosis"].isna().any():
        raise ValueError("Există valori nenumerice în coloana Diagnosis după mapare!")

    # ======================
    # 3) Standardizare features
    # ======================
    # X = toate coloanele numerice fără Diagnosis (eticheta "reală")
    X = df.drop(columns=["Diagnosis"])
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    # Director pentru rezultate
    output_dir = Path(f"/workspaces/bioinf-y4-lab/labs/05_clustering/submissions/{HANDLE}")
    output_dir.mkdir(parents=True, exist_ok=True)

    # ======================
    # 4) Hierarchical clustering + dendrogramă
    # ======================
    print("[INFO] Hierarchical clustering (average linkage)...")
    Z = linkage(X_scaled, method="average")

    plt.figure(figsize=(12, 6))
    dendrogram(
        Z,
        truncate_mode="lastp",  # ca să nu fie prea haotic (poți scoate dacă vrei tot)
        p=30,                   # ultimele 30 de clustere
        leaf_rotation=90.,
        leaf_font_size=8.,
        show_contracted=True,
    )
    plt.title("Hierarchical Clustering (average linkage) — WDBC")
    plt.xlabel("Cluster / Samples")
    plt.ylabel("Distanță")

    hierarchical_png = output_dir / f"hierarchical_{HANDLE}.png"
    plt.tight_layout()
    plt.savefig(hierarchical_png, dpi=300)
    plt.close()
    print(f"[OK] Dendrogramă salvată în: {hierarchical_png}")

    # ======================
    # 5) K-means (K=2) + PCA vizualizare 2D
    # ======================
    print("[INFO] KMeans clustering (K=2)...")
    kmeans = KMeans(n_clusters=2, random_state=42, n_init=10)
    kmeans_labels = kmeans.fit_predict(X_scaled)
    df["KMeans_Cluster"] = kmeans_labels

    # PCA 2D pentru vizualizare
    pca = PCA(n_components=2, random_state=42)
    X_pca = pca.fit_transform(X_scaled)

    plt.figure(figsize=(8, 6))
    scatter = plt.scatter(
        X_pca[:, 0],
        X_pca[:, 1],
        c=kmeans_labels,
        cmap="viridis",
        s=10,
        alpha=0.8,
    )
    plt.title("KMeans (K=2) pe WDBC — Proiecție PCA 2D")
    plt.xlabel("PC1")
    plt.ylabel("PC2")
    plt.colorbar(scatter, label="KMeans cluster")

    kmeans_png = output_dir / f"kmeans_{HANDLE}.png"
    plt.tight_layout()
    plt.savefig(kmeans_png, dpi=300)
    plt.close()
    print(f"[OK] Plot KMeans salvat în: {kmeans_png}")

    # ======================
    # 6) DBSCAN + PCA vizualizare 2D
    # ======================
    print("[INFO] DBSCAN clustering...")
    dbscan = DBSCAN(eps=1.5, min_samples=5)  # valori "toy", le poți ajusta
    db_labels = dbscan.fit_predict(X_scaled)
    df["DBSCAN_Cluster"] = db_labels

    plt.figure(figsize=(8, 6))
    scatter_db = plt.scatter(
        X_pca[:, 0],
        X_pca[:, 1],
        c=db_labels,
        cmap="tab10",
        s=10,
        alpha=0.8,
    )
    plt.title("DBSCAN pe WDBC — Proiecție PCA 2D")
    plt.xlabel("PC1")
    plt.ylabel("PC2")
    plt.colorbar(scatter_db, label="DBSCAN cluster (-1 = noise)")

    dbscan_png = output_dir / f"dbscan_{HANDLE}.png"
    plt.tight_layout()
    plt.savefig(dbscan_png, dpi=300)
    plt.close()
    print(f"[OK] Plot DBSCAN salvat în: {dbscan_png}")

    # ======================
    # 7) Salvare rezultate
    # ======================
    clusters_csv = output_dir / f"clusters_{HANDLE}.csv"
    df[["Diagnosis", "KMeans_Cluster", "DBSCAN_Cluster"]].to_csv(
        clusters_csv,
        index=False,
    )
    print(f"[OK] Rezultatele de clustering salvate în: {clusters_csv}")

    print("\n[SUMMARY]")
    print(df[["Diagnosis", "KMeans_Cluster", "DBSCAN_Cluster"]].head())