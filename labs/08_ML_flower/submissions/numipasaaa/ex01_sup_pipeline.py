"""
Exercise 8 — Supervised ML pipeline pentru expresie genică (Random Forest)

TODO-uri principale:
- Încărcați matricea de expresie (ex. subset TP53 / GTEx) pentru HANDLE-ul vostru
- Separați features (gene) și label (ultima coloană)
- Encodați etichetele
- Împărțiți în train/test
- Antrenați un RandomForestClassifier (model de bază)
- Evaluați: classification_report + matrice de confuzie (salvate)
- Calculați importanța trăsăturilor și salvați în CSV
- (Opțional) Aplicați KMeans pe X și comparați clustere vs etichete reale
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
HANDLE = "numipasaaa"

DATA_CSV = Path(f"data/work/{HANDLE}/lab06/fisier_mic2.csv")

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
    """
    TODO:
    - verificați că fișierul de input există
    - dacă nu, ridicați FileNotFoundError cu un mesaj clar
    """
    if not path.is_file():
        raise FileNotFoundError(f"Nu am găsit fișierul: {path}")
    pass


def load_dataset(path: Path) -> Tuple[pd.DataFrame, pd.Series]:
    """
    TODO:
    - citiți CSV-ul cu pandas
    - presupunem că ultima coloană este label-ul (ex. 'Label')
    - X = toate coloanele mai puțin ultima
    - y = ultima coloană
    """
    df = pd.read_csv(path, index_col=0)
    X = df.iloc[:, :-1]
    y = df.iloc[:, -1]
    
    return X, y


def drop_rare_classes(
    X: pd.DataFrame,
    y: pd.Series,
    min_count: int = 2,
) -> Tuple[pd.DataFrame, pd.Series, pd.Index]:
    """
    Remove classes with fewer than `min_count` samples to enable stratified splits.
    Returns filtered X, y and the index of dropped label names.
    """
    counts = y.value_counts()
    keep_labels = counts[counts >= min_count].index
    drop_labels = counts[counts < min_count].index
    if len(drop_labels) > 0:
        mask = y.isin(keep_labels)
        return X.loc[mask], y.loc[mask], drop_labels
    return X, y, pd.Index([])


def encode_labels(y: pd.Series) -> Tuple[np.ndarray, LabelEncoder]:
    """
    TODO:
    - folosiți LabelEncoder pentru a converti etichetele string în valori numerice
    - returnați y_encoded și encoder-ul
    """
    le = LabelEncoder()
    y_enc = le.fit_transform(y)
    return y_enc, le


def train_random_forest(
    X_train: pd.DataFrame,
    y_train: np.ndarray,
    n_estimators: int,
    random_state: int,
) -> RandomForestClassifier:
    """
    TODO:
    - inițializați RandomForestClassifier
    - antrenați modelul pe X_train, y_train
    - returnați modelul
    """
    rf = RandomForestClassifier(
        n_estimators=n_estimators,
        random_state=random_state,
        n_jobs=-1
    )

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
    """
    TODO:
    - calculați predicțiile pe X_test
    - generați classification_report (print + save)
    - calculați matricea de confuzie și salvați imaginea .png
    """
    y_pred = model.predict(X_test)

    # Ensure classification_report/CM use consistent label indices and names
    # Cast class names to strings to satisfy sklearn's formatting
    target_names = [str(cls) for cls in label_encoder.classes_]
    labels = np.arange(len(target_names))

    report = classification_report(
        y_test,
        y_pred,
        labels=labels,
        target_names=target_names,
        zero_division=0,
    )
    print(report)
    out_txt.write_text(report)
    
    cm = confusion_matrix(y_test, y_pred, labels=labels)
    plt.figure(figsize=(5, 4))
    sns.heatmap(
        cm,
        annot=True,
        fmt="d",
        cmap="Blues",
        xticklabels=target_names,
        yticklabels=target_names,
    )
    plt.xlabel("Predicted")
    plt.ylabel("Actual")
    plt.title("Random Forest — confusion matrix")
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close()
    pass


def compute_feature_importance(
    model: RandomForestClassifier,
    feature_names: pd.Index,
    out_csv: Path,
) -> pd.DataFrame:
    """
    TODO:
    - extrageți feature_importances_ din model
    - construiți un DataFrame cu coloane: 'Feature', 'Importance'
    - sortați descrescător
    - salvați în CSV și returnați DataFrame-ul
    """
    importances = model.feature_importances_
    df_imp = pd.DataFrame(
        {"Feature": feature_names, "Importance": importances}
    ).sort_values("Importance", ascending=False)
    df_imp.to_csv(out_csv, index=False)
    return df_imp


def run_kmeans_and_crosstab(
    X: pd.DataFrame,
    y: np.ndarray,
    label_encoder: LabelEncoder,
    n_clusters: int,
    out_csv: Path,
) -> None:
    """
    (Opțional) TODO:
    - rulați KMeans cu n_clusters egal cu numărul de clase
    - construiți un crosstab între eticheta reală și cluster
    - salvați crosstab-ul în CSV
    """
    kmeans = KMeans(n_clusters=n_clusters, random_state=RANDOM_STATE, n_init="auto")
    clusters = kmeans.fit_predict(X.values)
    
    df = pd.DataFrame(
        {"Label": label_encoder.inverse_transform(y), "Cluster": clusters}
    )
    ctab = pd.crosstab(df["Label"], df["Cluster"])
    ctab.to_csv(out_csv)
    print("Crosstab Label vs Cluster:")
    print(ctab)
    pass


# --------------------------
# Main
# --------------------------
if __name__ == "__main__":
    # TODO 1: verificați fișierul de input
    ensure_exists(DATA_CSV)

    # TODO 2: încărcați datele (X, y)
    X, y = load_dataset(DATA_CSV)

    # Drop classes with fewer than 2 samples to allow stratified split
    X, y, dropped = drop_rare_classes(X, y, min_count=2)
    if len(dropped) > 0:
        print(f"[WARN] Dropped classes with <2 samples: {list(dropped)}")

    # TODO 3: encodați etichetele și împărțiți în train/test
    y_enc, le = encode_labels(y)
    try:
        X_train, X_test, y_train, y_test = train_test_split(
            X, y_enc,
            test_size=TEST_SIZE,
            random_state=RANDOM_STATE,
            stratify=y_enc,
        )
    except ValueError as e:
        print(f"[WARN] Stratified split failed: {e}. Falling back to non-stratified split.")
        X_train, X_test, y_train, y_test = train_test_split(
            X, y_enc,
            test_size=TEST_SIZE,
            random_state=RANDOM_STATE,
            stratify=None,
        )

    # TODO 4: antrenați modelul RF și evaluați
    rf = train_random_forest(X_train, y_train, N_ESTIMATORS, RANDOM_STATE)
    evaluate_model(rf, X_test, y_test, le, OUT_CONFUSION, OUT_REPORT)

    # TODO 5: calculați importanța trăsăturilor
    feat_imp_df = compute_feature_importance(rf, X.columns, OUT_FEATIMP)

    # TODO 6 (opțional): rulați KMeans și salvați crosstab-ul
    n_classes = len(le.classes_)
    run_kmeans_and_crosstab(X, y_enc, le, n_clusters=n_classes, out_csv=OUT_CLUSTER_CROSSTAB)

