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
import gc 

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
HANDLE = "coroiusamy"

DATA_CSV = Path(f"data/work/{HANDLE}/lab08/expression_matrix_{HANDLE}.csv")

TEST_SIZE = 0.2
RANDOM_STATE = 42
N_ESTIMATORS = 50 
TOPK_FEATURES = 20  

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


def load_dataset(path: Path) -> Tuple[pd.DataFrame, pd.Series]:
    """
    TODO:
    - citiți CSV-ul cu pandas
    - presupunem că ultima coloană este label-ul (ex. 'Label')
    - X = toate coloanele mai puțin ultima
    - y = ultima coloană
    """
    print("[INFO] Loading dataset (Robusted for strings)...")
    
    # 1. Header only
    df_head = pd.read_csv(path, nrows=0)
    all_cols = df_head.columns.tolist()
    
    # 2. Selectie 
    MAX_COLS = 5000
    target_col = all_cols[-1]
    feature_cols = all_cols[:MAX_COLS]
    use_cols = feature_cols + [target_col]
    
    # 3. Citire
    df = pd.read_csv(path, usecols=use_cols)
    
    # 4. Selectăm DOAR coloanele numerice pentru X
    y = df.iloc[:, -1]
    X_temp = df.iloc[:, :-1]
    
    # select_dtypes pastreaza doar numerele, aruncand ID-urile text
    X = X_temp.select_dtypes(include=[np.number]).astype("float32")
    
    del df
    del X_temp
    gc.collect()
    
    print(f"[INFO] Dataset loaded. X shape: {X.shape}")
    return X, y


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
    print(f"[INFO] Training RF (n_jobs=1)...")
    rf = RandomForestClassifier(
        n_estimators=n_estimators,
        random_state=random_state,
        n_jobs=1
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

    target_names = label_encoder.classes_
    
    present_labels = np.unique(np.concatenate((y_test, y_pred)))
    
    # Ensure target names are strings to avoid TypeError: object has no len()
    target_names_present = [str(target_names[i]) for i in present_labels]

    report = classification_report(
        y_test, 
        y_pred, 
        labels=present_labels, 
        target_names=target_names_present,
        zero_division=0 
    )
    print(report)
    out_txt.write_text(report)

    cm = confusion_matrix(y_test, y_pred, labels=present_labels)
    plt.figure(figsize=(6, 5))
    sns.heatmap(
        cm,
        annot=True,
        fmt="d",
        cmap="Blues",
        xticklabels=target_names_present,
        yticklabels=target_names_present,
    )
    plt.xlabel("Predicted")
    plt.ylabel("Actual")
    plt.title("Random Forest — confusion matrix")
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close()


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


# --------------------------
# Main
# --------------------------
if __name__ == "__main__":
    # TODO 1: verificați fișierul de input
    ensure_exists(DATA_CSV)

    # TODO 2: încărcați datele (X, y)
    X, y = load_dataset(DATA_CSV)

    # FIX: Filtrare clase rare (mai puțin de 2 probe)
    print("[INFO] Filtering rare classes...")
    class_counts = y.value_counts()
    valid_classes = class_counts[class_counts > 1].index
    mask = y.isin(valid_classes)

    if not mask.all():
        print(f"[WARN] Eliminam {(~mask).sum()} probe cu clase unice (prea puține pentru stratify).")
        X = X[mask]
        y = y[mask]

    # TODO 3: encodați etichetele și împărțiți în train/test
    y_enc, le = encode_labels(y)
    X_train, X_test, y_train, y_test = train_test_split(
        X, y_enc,
        test_size=TEST_SIZE,
        random_state=RANDOM_STATE,
        stratify=y_enc,
    )

    # TODO 4: antrenați modelul RF și evaluați
    rf = train_random_forest(X_train, y_train, N_ESTIMATORS, RANDOM_STATE)
    evaluate_model(rf, X_test, y_test, le, OUT_CONFUSION, OUT_REPORT)

    # TODO 5: calculați importanța trăsăturilor
    feat_imp_df = compute_feature_importance(rf, X.columns, OUT_FEATIMP)

    # TODO 6 (opțional): rulați KMeans și salvați crosstab-ul
    n_classes = len(le.classes_)
    run_kmeans_and_crosstab(X, y_enc, le, n_clusters=n_classes, out_csv=OUT_CLUSTER_CROSSTAB)

    print("[INFO] Skeleton Exercise 8 — completați TODO-urile și rulați scriptul.")