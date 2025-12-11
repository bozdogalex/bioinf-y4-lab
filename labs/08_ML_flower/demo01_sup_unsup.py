"""
Demo 1 — Supervised vs Unsupervised pe date sintetice 2D

Scop:
- să vedem clar diferența între:
  - clasificare supravegheată (Random Forest cu etichete)
  - clustering nesupravegheat (KMeans fără etichete)
- să comparăm clusterele cu etichetele reale (crosstab).
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.datasets import make_classification
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score
from sklearn.model_selection import train_test_split

OUT_DIR = Path("labs/08_ml/demo_outputs")
OUT_DIR.mkdir(parents=True, exist_ok=True)
OUT_FIG_SUP = OUT_DIR / "demo_synthetic_supervised.png"
OUT_FIG_UNSUP = OUT_DIR / "demo_synthetic_unsupervised.png"


def main() -> None:
    # 1) Date sintetice 2D cu 3 clase reale
    X, y = make_classification(
        n_samples=300,
        n_features=2,
        n_redundant=0,
        n_informative=2,
        n_clusters_per_class=1,
        n_classes=3,
        random_state=42,
    )

    # 2) Supervised: RF pe train/test
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.3, random_state=42, stratify=y
    )

    rf = RandomForestClassifier(n_estimators=100, random_state=42)
    rf.fit(X_train, y_train)
    y_pred = rf.predict(X_test)
    acc = accuracy_score(y_test, y_pred)
    print(f"[Supervised] Test accuracy: {acc:.3f}")

    # vizualizare supervised
    cmap_bold = ["#FF0000", "#00AA00", "#0000FF"]
    plt.figure(figsize=(6, 5))
    for cls in np.unique(y_train):
        mask = y_train == cls
        plt.scatter(
            X_train[mask, 0],
            X_train[mask, 1],
            c=cmap_bold[cls],
            label=f"Class {cls}",
            edgecolor="k",
            s=40,
        )
    plt.title("Supervised — puncte colorate după eticheta reală (train)")
    plt.xlabel("Feature 1")
    plt.ylabel("Feature 2")
    plt.legend()
    plt.tight_layout()
    plt.savefig(OUT_FIG_SUP, dpi=200)
    plt.show()

    # 3) Unsupervised: KMeans pe toate punctele (fără etichete)
    kmeans = KMeans(n_clusters=3, random_state=42, n_init="auto")
    clusters = kmeans.fit_predict(X)

    # vizualizare unsupervised
    cmap_unsup = ["#FFAAAA", "#AAFFAA", "#AAAFFF"]
    plt.figure(figsize=(6, 5))
    for cl in np.unique(clusters):
        mask = clusters == cl
        plt.scatter(
            X[mask, 0],
            X[mask, 1],
            c=cmap_unsup[cl],
            label=f"Cluster {cl}",
            edgecolor="k",
            s=40,
        )
    plt.title("Unsupervised — KMeans clusters (fără etichete)")
    plt.xlabel("Feature 1")
    plt.ylabel("Feature 2")
    plt.legend()
    plt.tight_layout()
    plt.savefig(OUT_FIG_UNSUP, dpi=200)
    plt.show()

    # 4) Comparație clustere vs etichete reale
    df = pd.DataFrame({"TrueLabel": y, "Cluster": clusters})
    ctab = pd.crosstab(df["TrueLabel"], df["Cluster"])
    print("\nCrosstab TrueLabel vs Cluster:")
    print(ctab)


if __name__ == "__main__":
    main()
