"""
Exercise 8b — Logistic Regression vs Random Forest pe expresie genică

Scop:
- să antrenăm și să comparăm două modele:
  - Logistic Regression (multiclass, liniar)
  - Random Forest (non-liniar, bazat pe arbori)
- să vedem dacă performanța și erorile sunt similare sau diferite

TODO:
- Încărcați expresia pentru HANDLE
- Împărțiți în X (gene) și y (Label)
- Encodați etichetele
- Împărțiți în train/test
- Scalați features pentru logistic regression
- Antrenați RF și Logistic Regression
- Comparați classification_report pentru ambele modele
"""

from __future__ import annotations
from pathlib import Path
from typing import Tuple

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import classification_report
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelEncoder, StandardScaler

# --------------------------
# Config
# --------------------------
HANDLE = "numipasaaa"

DATA_CSV = Path(f"data/work/{HANDLE}/lab06/fisier_mic2.csv")

TEST_SIZE = 0.2
RANDOM_STATE = 42
N_ESTIMATORS = 200
MAX_ITER_LOGREG = 1000

OUT_DIR = Path(f"labs/08_ML_flower/submissions/{HANDLE}")
OUT_DIR.mkdir(parents=True, exist_ok=True)

OUT_REPORT_TXT = OUT_DIR / f"rf_vs_logreg_report_{HANDLE}.txt"


# --------------------------
# Utils
# --------------------------
def ensure_exists(path: Path) -> None:
    """
    TODO:
    - verificați că fișierul există
    - dacă nu, ridicați excepție
    """
    if not path.is_file():
        raise FileNotFoundError(f"Nu am găsit fișierul: {path}")
    pass


def load_dataset(path: Path) -> Tuple[pd.DataFrame, pd.Series]:
    """
    TODO:
    - citiți CSV cu pandas
    - X = toate coloanele mai puțin ultima
    - y = ultima coloană (Label)
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
    - folosiți LabelEncoder pentru a obține y_enc
    """
    le = LabelEncoder()
    y_enc = le.fit_transform(y)
    return y_enc, le


def train_models(
    X_train: pd.DataFrame,
    y_train: np.ndarray,
) -> Tuple[RandomForestClassifier, LogisticRegression, StandardScaler]:
    """
    TODO:
    - antrenați două modele:
      - RandomForestClassifier
      - LogisticRegression (cu scaling înainte)
    - întoarceți (rf, logreg, scaler)
    """
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    
    rf = RandomForestClassifier(
        n_estimators=N_ESTIMATORS,
        random_state=RANDOM_STATE,
        n_jobs=-1
    )
    rf.fit(X_train, y_train)
    
    logreg = LogisticRegression(
        multi_class="multinomial",
        max_iter=MAX_ITER_LOGREG,
        n_jobs=-1
    )
    logreg.fit(X_train_scaled, y_train)
    
    return rf, logreg, scaler


def compare_models(
    rf: RandomForestClassifier,
    logreg: LogisticRegression,
    scaler: StandardScaler,
    X_test: pd.DataFrame,
    y_test: np.ndarray,
    label_encoder: LabelEncoder,
    out_txt: Path,
) -> None:
    """
    TODO:
    - calculați predicții pentru ambele modele
    - generați classification_report pentru RF și pentru Logistic Regression
    - scrieți într-un singur fișier .txt, cu secțiuni separate
    """
    X_test_scaled = scaler.transform(X_test)
    
    y_pred_rf = rf.predict(X_test)
    y_pred_logreg = logreg.predict(X_test_scaled)

    # Ensure classification_report/CM use consistent label indices and names
    # Cast class names to strings to satisfy sklearn's formatting
    target_names = [str(cls) for cls in label_encoder.classes_]
    labels = np.arange(len(target_names))
    
    report_rf = classification_report(y_test, y_pred_rf, labels=labels, target_names=target_names)
    report_logreg = classification_report(y_test, y_pred_logreg, labels=labels, target_names=target_names)
    
    print("=== Random Forest ===")
    print(report_rf)
    print("\n=== Logistic Regression ===")
    print(report_logreg)
    
    combined = (
        "=== Random Forest ===\n"
        + report_rf
        + "\n\n=== Logistic Regression ===\n"
        + report_logreg
    )
    out_txt.write_text(combined)
    pass


# --------------------------
# Main
# --------------------------
if __name__ == "__main__":
    # TODO 1: verificați fișierul
    ensure_exists(DATA_CSV)

    # TODO 2: încărcați X, y
    X, y = load_dataset(DATA_CSV)

    # Drop classes with fewer than 2 samples to allow stratified split
    X, y, dropped = drop_rare_classes(X, y, min_count=2)
    if len(dropped) > 0:
        print(f"[WARN] Dropped classes with <2 samples: {list(dropped)}")

    # TODO 3: encodați etichetele și împărțiți în train/test
    y_enc, le = encode_labels(y)
    X_train, X_test, y_train, y_test = train_test_split(
        X, y_enc,
        test_size=TEST_SIZE,
        random_state=RANDOM_STATE,
        stratify=y_enc,
    )

    # TODO 4: antrenați ambele modele
    rf, logreg, scaler = train_models(X_train, y_train)

    # TODO 5: comparați modelele și salvați raportul
    compare_models(rf, logreg, scaler, X_test, y_test, le, OUT_REPORT_TXT)

