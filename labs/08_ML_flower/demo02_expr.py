"""
Demo 2 — Supervised + Unsupervised pe expresie genică (subset GTEx / TP53-like)

Scop:
- Supervised: Random Forest pentru clasificarea probelor (Label = tip țesut / status)
- Unsupervised: PCA + KMeans pe aceleași probe
- Comparație: clustere vs etichete (crosstab)

Fișier de intrare:
- labs/08_ml/data/expression_demo.csv
  - rânduri: probe
  - coloane: gene + ultima coloană Label
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report, confusion_matrix
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelEncoder

DATA_CSV = Path("../data/sample/tissue_gene_expression_demo.csv")
OUT_DIR = Path("labs/08_ml/demo_outputs")
OUT_DIR.mkdir(parents=True, exist_ok=True)

OUT_CONFUSION = OUT_DIR / "demo_expr_confusion.png"
OUT_FEATIMP = OUT_DIR / "demo_expr_featimportance.png"
OUT_PCA_SUP = OUT_DIR / "demo_expr_pca_supervised.png"
OUT_PCA_UNSUP = OUT_DIR / "demo_expr_pca_unsupervised.png"


def main() -> None:
    print(f"[INFO] Loading demo data from {DATA_CSV}")
    df = pd.read_csv(DATA_CSV)
    print(df.head())

    X = df.iloc[:, :-1]
    y = df.iloc[:, -1]

    le = LabelEncoder()
    y_enc = le.fit_transform(y)
    classes = le.classes_
    print("[INFO] Classes:", list(classes))

    # --- Supervised: RF -----------------------------------------------------
    X_train, X_test, y_train, y_test = train_test_split(
        X, y_enc, test_size=0.2, random_state=42, stratify=y_enc
    )

    rf = RandomForestClassifier(
        n_estimators=200,
        random_state=42,
        n_jobs=-1
    )
    rf.fit(X_train, y_train)
    y_pred = rf.predict(X_test)

    print("\n=== Classification Report (RF) ===")
    print(classification_report(y_test, y_pred, target_names=classes))

    cm = confusion_matrix(y_test, y_pred)
    plt.figure(figsize=(5, 4))
    sns.heatmap(
        cm,
        annot=True,
        fmt="d",
        cmap="Blues",
        xticklabels=classes,
        yticklabels=classes,
    )
    plt.xlabel("Predicted")
    plt.ylabel("Actual")
    plt.title("Demo — RF confusion matrix")
    plt.tight_layout()
    plt.savefig(OUT_CONFUSION, dpi=200)
    plt.show()

    # Feature importance
    importances = rf.feature_importances_
    df_imp = pd.DataFrame(
        {"Gene": X.columns, "Importance": importances}
    ).sort_values("Importance", ascending=False)
    plt.figure(figsize=(8, 5))
    df_imp.head(20).plot(kind="barh", x="Gene", y="Importance")
    plt.gca().invert_yaxis()
    plt.title("Demo — Top 20 genes by RF importance")
    plt.tight_layout()
    plt.savefig(OUT_FEATIMP, dpi=200)
    plt.show()

    # --- PCA + KMeans (Unsupervised) --------------------------------------
    pca = PCA(n_components=2, random_state=42)
    X_pca = pca.fit_transform(X.values)

    kmeans = KMeans(n_clusters=len(classes), random_state=42, n_init="auto")
    clusters = kmeans.fit_predict(X.values)

    # PCA plot colorat după eticheta reală (supervised view)
    plt.figure(figsize=(6, 5))
    for idx, cls in enumerate(classes):
        mask = (y_enc == idx)
        plt.scatter(
            X_pca[mask, 0],
            X_pca[mask, 1],
            label=f"Label: {cls}",
            s=40,
            alpha=0.8,
        )
    plt.title("PCA — puncte colorate după eticheta reală")
    plt.xlabel("PC1")
    plt.ylabel("PC2")
    plt.legend()
    plt.tight_layout()
    plt.savefig(OUT_PCA_SUP, dpi=200)
    plt.show()

    # PCA plot colorat după cluster (unsupervised view)
    plt.figure(figsize=(6, 5))
    for cl in np.unique(clusters):
        mask = (clusters == cl)
        plt.scatter(
            X_pca[mask, 0],
            X_pca[mask, 1],
            label=f"Cluster {cl}",
            s=40,
            alpha=0.8,
        )
    plt.title("PCA — puncte colorate după cluster (KMeans)")
    plt.xlabel("PC1")
    plt.ylabel("PC2")
    plt.legend()
    plt.tight_layout()
    plt.savefig(OUT_PCA_UNSUP, dpi=200)
    plt.show()

    # Crosstab clustere vs etichete reale
    df_res = pd.DataFrame({"Label": y, "Cluster": clusters})
    ctab = pd.crosstab(df_res["Label"], df_res["Cluster"])
    print("\nCrosstab Label vs Cluster:")
    print(ctab)


if __name__ == "__main__":
    main()
    print("[INFO] Demo finished.")
