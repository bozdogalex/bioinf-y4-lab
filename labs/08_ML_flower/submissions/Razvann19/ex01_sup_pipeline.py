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
HANDLE = "Razvann19"

DATA_CSV = Path("labs/08_ML_flower/submissions/Razvann19/expression_matrix.csv")

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
    Verificați că fișierul de input există.
    Dacă nu, ridicați FileNotFoundError cu un mesaj clar.
    """
    if not path.is_file():
        raise FileNotFoundError(f"[ERROR] Nu am găsit fișierul de input: {path}")
    print(f"[INFO] Găsit fișierul de input: {path}")


def load_dataset(path: Path) -> Tuple[pd.DataFrame, pd.Series]:
    """
    - Citim CSV-ul cu pandas
    - Presupunem că există o coloană 'Label' cu clase (low / medium / high)
    - X = toate coloanele numerice, fără 'Label'
    - y = coloana 'Label'
    """
    print(f"[INFO] Încărcăm dataset-ul din {path}")
    df = pd.read_csv(path)
    print("[INFO] Shape dataframe:", df.shape)
    print("[INFO] Coloane:", list(df.columns))

    # y = eticheta
    if "Label" not in df.columns:
        raise ValueError("Nu există coloana 'Label' în dataset!")

    y = df["Label"]

    # X = toate coloanele în afară de Label
    X_all = df.drop(columns=["Label"])

    # păstrăm DOAR coloanele numerice ca features
    X = X_all.select_dtypes(include=[np.number])

    print("[INFO] Număr features numerice:", X.shape[1])
    print("[INFO] Primele etichete unice:", y.unique())
    return X, y



def encode_labels(y: pd.Series) -> Tuple[np.ndarray, LabelEncoder]:
    """
    Folosim LabelEncoder pentru a converti etichetele string în valori numerice.
    Returnăm y_encoded și encoder-ul.
    """
    le = LabelEncoder()
    y_enc = le.fit_transform(y)
    print("[INFO] Clase encodate:", list(le.classes_))
    return y_enc, le


def train_random_forest(
    X_train: pd.DataFrame,
    y_train: np.ndarray,
    n_estimators: int,
    random_state: int,
) -> RandomForestClassifier:
    """
    - Inițializăm RandomForestClassifier
    - Antrenăm modelul pe X_train, y_train
    - Returnăm modelul
    """
    print("[INFO] Antrenăm RandomForestClassifier...")
    rf = RandomForestClassifier(
        n_estimators=n_estimators,
        random_state=random_state,
        n_jobs=-1,
    )
    rf.fit(X_train, y_train)
    print("[INFO] Antrenare RF finalizată.")
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
    - Calculăm predicțiile pe X_test
    - Generăm classification_report (print + save)
    - Calculăm matricea de confuzie și salvăm imaginea .png
    """
    print("[INFO] Evaluăm modelul pe setul de test...")
    y_pred = model.predict(X_test)

    target_names = label_encoder.classes_
    report = classification_report(y_test, y_pred, target_names=target_names)
    print("\n=== Classification Report (Random Forest) ===")
    print(report)
    out_txt.write_text(report)
    print(f"[INFO] Salvat classification_report în {out_txt}")

    cm = confusion_matrix(y_test, y_pred)
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
    print(f"[INFO] Salvat matricea de confuzie în {out_png}")


def compute_feature_importance(
    model: RandomForestClassifier,
    feature_names: pd.Index,
    out_csv: Path,
) -> pd.DataFrame:
    """
    - Extragem feature_importances_ din model
    - Construim un DataFrame cu coloane: 'Feature', 'Importance'
    - Sortăm descrescător
    - Salvăm în CSV și returnăm DataFrame-ul
    """
    print("[INFO] Calculăm importanța trăsăturilor (feature_importances_)...")
    importances = model.feature_importances_
    df_imp = pd.DataFrame(
        {"Feature": feature_names, "Importance": importances}
    ).sort_values("Importance", ascending=False)

    df_imp.to_csv(out_csv, index=False)
    print(f"[INFO] Salvat feature importance în {out_csv}")
    return df_imp


def run_kmeans_and_crosstab(
    X: pd.DataFrame,
    y: np.ndarray,
    label_encoder: LabelEncoder,
    n_clusters: int,
    out_csv: Path,
) -> None:
    """
    (Opțional):
    - Rulăm KMeans cu n_clusters egal cu numărul de clase
    - Construim un crosstab între eticheta reală și cluster
    - Salvăm crosstab-ul în CSV
    """
    print(f"[INFO] Rulăm KMeans cu n_clusters={n_clusters}...")
    kmeans = KMeans(n_clusters=n_clusters, random_state=RANDOM_STATE, n_init="auto")
    clusters = kmeans.fit_predict(X.values)

    df = pd.DataFrame(
        {"Label": label_encoder.inverse_transform(y), "Cluster": clusters}
    )
    ctab = pd.crosstab(df["Label"], df["Cluster"])
    ctab.to_csv(out_csv)
    print("[INFO] Crosstab Label vs Cluster:")
    print(ctab)
    print(f"[INFO] Salvat crosstab-ul în {out_csv}")


# --------------------------
# Main
# --------------------------
if __name__ == "__main__":
    # TODO 1: verificați fișierul de input
    ensure_exists(DATA_CSV)

    # TODO 2: încărcați datele (X, y)
    X, y = load_dataset(DATA_CSV)

    # TODO 3: encodați etichetele și împărțiți în train/test
    y_enc, le = encode_labels(y)
    X_train, X_test, y_train, y_test = train_test_split(
        X,
        y_enc,
        test_size=TEST_SIZE,
        random_state=RANDOM_STATE,
        stratify=y_enc,
    )

    # TODO 4: antrenați modelul RF și evaluați
    rf = train_random_forest(X_train, y_train, N_ESTIMATORS, RANDOM_STATE)
    evaluate_model(rf, X_test, y_test, le, OUT_CONFUSION, OUT_REPORT)

    # TODO 5: calculați importanța trăsăturilor
    feat_imp_df = compute_feature_importance(rf, X.columns, OUT_FEATIMP)

    # (Opțional) top-k features afișate în consolă
    if TOPK_FEATURES is not None and TOPK_FEATURES > 0:
        print(f"\n[INFO] Top {TOPK_FEATURES} features după importanță:")
        print(feat_imp_df.head(TOPK_FEATURES))

    # TODO 6 (opțional): rulați KMeans și salvați crosstab-ul
    n_classes = len(le.classes_)
    run_kmeans_and_crosstab(X, y_enc, le, n_clusters=n_classes, out_csv=OUT_CLUSTER_CROSSTAB)

    print("[INFO] Exercise 8 terminat cu succes.")
