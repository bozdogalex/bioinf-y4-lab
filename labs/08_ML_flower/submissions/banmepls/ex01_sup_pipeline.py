"""
Exercise 8 — Supervised ML pipeline pentru expresie genică (Random Forest)
"""

from __future__ import annotations
from pathlib import Path
from typing import Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.cluster import KMeans
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report, confusion_matrix
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelEncoder

# --------------------------
# Config — completați cu valorile voastre
# --------------------------
HANDLE = "banmepls"

DATA_CSV = Path(f"data/work/{HANDLE}/lab08/expression_matrix_{HANDLE}.csv")

TEST_SIZE = 0.2
RANDOM_STATE = 42
N_ESTIMATORS = 200
TOPK_FEATURES = 20  # opțional, pentru extensie

OUT_DIR = Path(f"labs/08_ML_flower/submissions/{HANDLE}")
OUT_DIR.mkdir(parents=True, exist_ok=True)

OUT_CONFUSION = OUT_DIR / f"confusion_rf_{HANDLE}.png"
OUT_REPORT = OUT_DIR / f"classification_report_{HANDLE}.txt"
OUT_FEATIMP = OUT_DIR / f"feature_importance_{HANDLE}.csv"
OUT_CLUSTER_CROSSTAB = OUT_DIR / f"cluster_crosstab_{HANDLE}.csv"


# --------------------------
# Utils
# --------------------------
def ensure_exists(path: Path) -> None:
    if not path.is_file():
        raise FileNotFoundError(f"Fișierul lipsește: {path}")


def load_dataset(path: Path) -> Tuple[pd.DataFrame, pd.Series]:
    df = pd.read_csv(path)
    
    label_col = df.columns[-1]
    X = df.drop(columns=[label_col])
    y = df[label_col]
    
    X = X.apply(pd.to_numeric, errors='coerce').dropna(axis=1, how='any')
    
    return X, y


def encode_labels(y: pd.Series) -> Tuple[np.ndarray, LabelEncoder]:
    le = LabelEncoder()
    y_encoded = le.fit_transform(y)
    return y_encoded, le


def train_random_forest(
    X_train: pd.DataFrame,
    y_train: np.ndarray,
    n_estimators: int,
    random_state: int,
) -> RandomForestClassifier:
    rf = RandomForestClassifier(n_estimators=N_ESTIMATORS, random_state=RANDOM_STATE, n_jobs=-1)
    rf.fit(X_train, y_train)
    return rf


def evaluate_model(
    model: RandomForestClassifier,
    X_test: pd.DataFrame,
    y_test: np.ndarray,
    label_encoder: LabelEncoder,
    out_png: Path,
    out_txt: Path,
) -> None:
    y_pred = model.predict(X_test)
    
    unique_labels = np.unique(np.concatenate([y_test, y_pred])).astype(int)
    present_classes = []
    for i in unique_labels:
        if 0 <= i < len(le.classes_):
            present_classes.append(str(le.classes_[i]))
        else:
            present_classes.append(f"Class_{i}")
    
    report = classification_report(
        y_test, y_pred,
        labels=unique_labels,
        target_names=present_classes,
        zero_division=0
    )
    
    print("\n=== Classification Report ===")
    print(report)
    
    with open(out_txt, "w", encoding="utf-8") as f:
        f.write(report)
    
    cm = confusion_matrix(y_test, y_pred, labels=unique_labels)
    plt.figure(figsize=(12, 10))
    sns.heatmap(cm, annot=True, fmt="d", cmap="Blues",
                xticklabels=present_classes, yticklabels=present_classes)
    plt.xlabel("Predicted")
    plt.ylabel("Actual")
    plt.title("Confusion Matrix - Random Forest")
    plt.xticks(rotation=90, fontsize=8)
    plt.yticks(rotation=0, fontsize=8)
    plt.tight_layout()
    plt.savefig(out_png, dpi=200, bbox_inches='tight')
    plt.close()


def compute_feature_importance(
    model: RandomForestClassifier,
    feature_names: pd.Index,
    out_csv: Path,
) -> pd.DataFrame:
    importances = model.feature_importances_
    df_imp = pd.DataFrame({
        "Gene": feature_names,
        "Importance": importances
    }).sort_values("Importance", ascending=False)
    
    df_imp.to_csv(out_csv, index=False)
    return df_imp


def run_kmeans_and_crosstab(
    X: pd.DataFrame,
    y: np.ndarray,
    label_encoder: LabelEncoder,
    n_clusters: int,
    out_csv: Path,
) -> None:
    n_classes = len(le.classes_)
    n_clusters = max(2, min(n_classes, len(np.unique(y_enc))))
    
    kmeans = KMeans(n_clusters=n_clusters, random_state=RANDOM_STATE, n_init="auto")
    clusters = kmeans.fit_predict(X)
    
    df = pd.DataFrame({
        "Label": le.inverse_transform(y_enc),
        "Cluster": clusters
    })
    ctab = pd.crosstab(df["Label"], df["Cluster"])
    ctab.to_csv(out_csv)
    
    print("\nCrosstab Label vs Cluster (KMeans):")
    print(ctab)


# --------------------------
# Main
# --------------------------
if __name__ == "__main__":
    ensure_exists(DATA_CSV)

    X, y = load_dataset(DATA_CSV)
    y_enc, le = encode_labels(y)
    
    X_train, X_test, y_train, y_test = train_test_split(
        X, y_enc,
        test_size=TEST_SIZE,
        random_state=RANDOM_STATE,
        shuffle=True
    )

    rf = train_random_forest(X_train, y_train, N_ESTIMATORS, RANDOM_STATE)
    evaluate_model(rf, X_test, y_test, le, OUT_CONFUSION, OUT_REPORT)

    feat_imp = compute_feature_importance(rf, X.columns, OUT_FEATIMP)

    n_classes = len(le.classes_)
    run_kmeans_and_crosstab(
        X=X,
        y=y_enc,
        label_encoder=le,
        n_clusters=n_classes,
        out_csv=OUT_CLUSTER_CROSSTAB
    )

    print(f"\n[OK] Totul gata! Rezultate în: {OUT_DIR}")
