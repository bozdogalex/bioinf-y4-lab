"""
Demo — Logistic Regression pe expresie genică (multiclass)

Scop:
- să arătăm un model liniar clasic (Logistic Regression) pe același tip de date
- să comparăm mental cu Random Forest (non-liniar, bazat pe arbori)
- să vedem classification_report + matrice de confuzie

Fișier de intrare:
- labs/08_ml/data/expression_demo.csv
  - rânduri: probe
  - coloane: gene (G1, ..., Gn) + ultima coloană 'Label'
"""

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import classification_report, confusion_matrix
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelEncoder, StandardScaler

DATA_CSV = Path("../data/sample/tissue_gene_expression_demo.csv")
OUT_DIR = Path("labs/08_ml/demo_outputs")
OUT_DIR.mkdir(parents=True, exist_ok=True)

OUT_CONFUSION = OUT_DIR / "demo_logreg_confusion.png"


def main() -> None:
    print(f"[INFO] Loading demo data from {DATA_CSV}")
    df = pd.read_csv(DATA_CSV)
    print(df.head())

    # X = gene, y = Label
    X = df.iloc[:, :-1]
    y = df.iloc[:, -1]

    # Encodare etichete
    le = LabelEncoder()
    y_enc = le.fit_transform(y)
    classes = le.classes_
    print("[INFO] Classes:", list(classes))

    # Împărțire train/test
    X_train, X_test, y_train, y_test = train_test_split(
        X, y_enc,
        test_size=0.2,
        random_state=42,
        stratify=y_enc,
    )

    # Scale features – logistic regression se comportă mai bine cu scaling
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)

    # Logistic Regression multiclass
    logreg = LogisticRegression(
        multi_class="multinomial",
        max_iter=1000,
        n_jobs=-1
    )
    logreg.fit(X_train_scaled, y_train)

    y_pred = logreg.predict(X_test_scaled)

    # Raport de clasificare
    print("\n=== Classification Report (Logistic Regression) ===")
    report = classification_report(y_test, y_pred, target_names=classes)
    print(report)

    # Matr. confuzie
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
    plt.title("Logistic Regression — confusion matrix (expression_demo)")
    plt.tight_layout()
    plt.savefig(OUT_CONFUSION, dpi=200)
    plt.show()

    print(f"[INFO] Saved confusion matrix to {OUT_CONFUSION}")


if __name__ == "__main__":
    main()
    print("[INFO] Logistic Regression demo finished.")
