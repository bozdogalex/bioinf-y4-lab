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
HANDLE = "Razvann19"

DATA_CSV = Path(f"labs/08_ML_flower/submissions/Razvann19/expression_matrix.csv")

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
    - verificați că fișierul există
    - dacă nu, ridicați excepție
    """
    if not path.is_file():
        raise FileNotFoundError(f"[ERROR] Nu am găsit fișierul: {path}")
    print(f"[INFO] Găsit fișierul de input: {path}")


def load_dataset(path: Path) -> Tuple[pd.DataFrame, pd.Series]:
    """
    - citiți CSV cu pandas
    - presupunem că există o coloană 'Label' pe ultima poziție
    - X = toate coloanele numerice, fără 'Label'
    - y = coloana 'Label'
    """
    print(f"[INFO] Încărcăm dataset-ul din {path}")
    df = pd.read_csv(path)
    print("[INFO] Shape dataframe:", df.shape)
    print("[INFO] Coloane:", list(df.columns))

    # y = ultima coloană (Label)
    y = df.iloc[:, -1]

    # Scoatem Label din features
    X_all = df.iloc[:, :-1]

    # Păstrăm DOAR coloanele numerice (genele de expresie),
    # aruncăm orice coloană cu string-uri (ex: IDs, nume gene)
    X = X_all.select_dtypes(include=[np.number])

    print("[INFO] Număr features numerice:", X.shape[1])
    print("[INFO] Primele etichete unice:", y.unique())
    return X, y



def encode_labels(y: pd.Series) -> Tuple[np.ndarray, LabelEncoder]:
    """
    - folosiți LabelEncoder pentru a obține y_enc
    """
    le = LabelEncoder()
    y_enc = le.fit_transform(y)
    print("[INFO] Clase encodate:", list(le.classes_))
    return y_enc, le


def train_models(
    X_train: pd.DataFrame,
    y_train: np.ndarray,
) -> Tuple[RandomForestClassifier, LogisticRegression, StandardScaler]:
    """
    - antrenați două modele:
      - RandomForestClassifier
      - LogisticRegression (cu scaling înainte)
    - întoarceți (rf, logreg, scaler)
    """
    print("[INFO] Antrenăm modelele RF și Logistic Regression...")

    # scaler pentru Logistic Regression
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)

    # Random Forest pe features ne-scalate
    rf = RandomForestClassifier(
        n_estimators=N_ESTIMATORS,
        random_state=RANDOM_STATE,
        n_jobs=-1,
    )
    rf.fit(X_train, y_train)

    # Logistic Regression pe features scalate
    logreg = LogisticRegression(
        multi_class="multinomial",
        max_iter=MAX_ITER_LOGREG,
        n_jobs=-1,
    )
    logreg.fit(X_train_scaled, y_train)

    print("[INFO] Antrenare modele finalizată.")
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
    - calculați predicții pentru ambele modele
    - generați classification_report pentru RF și pentru Logistic Regression
    - scrieți într-un singur fișier .txt, cu secțiuni separate
    """
    print("[INFO] Comparam modelele pe setul de test...")

    # scaling pentru Logistic Regression
    X_test_scaled = scaler.transform(X_test)

    # predicții
    y_pred_rf = rf.predict(X_test)
    y_pred_logreg = logreg.predict(X_test_scaled)

    target_names = label_encoder.classes_

    report_rf = classification_report(y_test, y_pred_rf, target_names=target_names)
    report_logreg = classification_report(y_test, y_pred_logreg, target_names=target_names)

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
    print(f"[INFO] Salvat raport comparativ în {out_txt}")


# --------------------------
# Main
# --------------------------
if __name__ == "__main__":
    # TODO 1: verificați fișierul
    ensure_exists(DATA_CSV)

    # TODO 2: încărcați X, y
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

    # TODO 4: antrenați ambele modele
    rf, logreg, scaler = train_models(X_train, y_train)

    # TODO 5: comparați modelele și salvați raportul
    compare_models(rf, logreg, scaler, X_test, y_test, le, OUT_REPORT_TXT)

    print("[INFO] Exercise 8b — RF vs Logistic Regression rulat cu succes.")
