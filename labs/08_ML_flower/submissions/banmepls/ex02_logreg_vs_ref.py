"""
Exercise 8b — Logistic Regression vs Random Forest pe expresie genică
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
HANDLE = "banmepls"

DATA_CSV = Path(f"data/work/{HANDLE}/lab08/expression_matrix_{HANDLE}.csv")

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
    if not path.is_file():
        raise FileNotFoundError(f"Fișierul lipsește: {path}")


def load_dataset(path: Path) -> Tuple[pd.DataFrame, pd.Series]:
    df = pd.read_csv(path)
    label_col = df.columns[-1]
    X = df.drop(columns=[label_col])
    y = df[label_col]
    
    # Convertim X în numeric (important!)
    X = X.apply(pd.to_numeric, errors='coerce').dropna(axis=1, how='any')
    
    if X.empty:
        raise ValueError("Nu au rămas coloane numerice în X!")
    
    print(f"Shape X: {X.shape} | Shape y: {y.shape}")
    return X, y


def encode_labels(y: pd.Series) -> Tuple[np.ndarray, LabelEncoder]:
    le = LabelEncoder()
    y_encoded = le.fit_transform(y)
    print("Clase detectate:", list(le.classes_))
    return y_encoded, le


def train_models(
    X_train: pd.DataFrame,
    y_train: np.ndarray,
) -> Tuple[RandomForestClassifier, LogisticRegression, StandardScaler]:
    # Random Forest (nu are nevoie de scaling)
    rf = RandomForestClassifier(
        n_estimators=N_ESTIMATORS,
        random_state=RANDOM_STATE,
        n_jobs=-1
    )
    rf.fit(X_train, y_train)
    
    # Logistic Regression (are nevoie de scaling)
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    
    logreg = LogisticRegression(
        multi_class="multinomial",  # implicit în viitor, dar explicit acum
        max_iter=MAX_ITER_LOGREG,
        random_state=RANDOM_STATE,
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
    # Predicții RF (fără scaling)
    y_pred_rf = rf.predict(X_test)
    
    # Predicții LogReg (cu scaling)
    X_test_scaled = scaler.transform(X_test)
    y_pred_logreg = logreg.predict(X_test_scaled)
    
    # Ia clasele unice prezente în y_test și predicții
    unique_labels = np.unique(np.concatenate([y_test, y_pred_rf, y_pred_logreg])).astype(int)
    
    # Construim target_names doar pentru clasele reale prezente
    present_classes = []
    for i in unique_labels:
        if 0 <= i < len(label_encoder.classes_):
            present_classes.append(str(label_encoder.classes_[i]))  # forțăm string
        else:
            present_classes.append(f"Class_{i}")
    
    # Raport RF
    report_rf = classification_report(
        y_test, y_pred_rf,
        labels=unique_labels,
        target_names=present_classes,
        zero_division=0
    )
    
    # Raport LogReg
    report_logreg = classification_report(
        y_test, y_pred_logreg,
        labels=unique_labels,
        target_names=present_classes,
        zero_division=0
    )
    
    # Afișare în consolă
    print("=== Random Forest ===")
    print(report_rf)
    print("\n=== Logistic Regression ===")
    print(report_logreg)
    
    # Salvăm totul într-un fișier
    combined = (
        "=== Random Forest ===\n"
        + report_rf
        + "\n\n=== Logistic Regression ===\n"
        + report_logreg
    )
    with open(out_txt, "w", encoding="utf-8") as f:
        f.write(combined)
    
    print(f"\nRaport comparativ salvat în: {out_txt}")


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
    
    rf, logreg, scaler = train_models(X_train, y_train)
    
    compare_models(rf, logreg, scaler, X_test, y_test, le, OUT_REPORT_TXT)
    
    print(f"\n[OK] Totul gata! Rezultate în: {OUT_DIR}")
