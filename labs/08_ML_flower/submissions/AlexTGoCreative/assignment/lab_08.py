"""
Lab 08 — Machine Learning pe date omice (Supervised, Unsupervised, Semi-Supervised)

Pipeline complet pentru:
- Task 1: Pregatirea datelor
- Task 2: Supervised ML (Random Forest)
- Task 3: Logistic Regression (bonus +1p)
- Task 4: Unsupervised ML (PCA + KMeans)
- Task 5: Semi-Supervised Learning
- Bonus: Comparatie PCA inainte si dupa eliminarea genelor cu varianta mica
"""

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import LabelEncoder, StandardScaler
from sklearn.model_selection import train_test_split
from sklearn.metrics import (
    classification_report,
    confusion_matrix,
    accuracy_score,
    f1_score,
)
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans

HANDLE = "AlexTGoCreative"
DATA_PATH = Path("/workspaces/bioinf-y4-lab/data/work/AlexTGoCreative/lab08/expression_matrix_1.csv")
OUT_DIR = Path("assignment")
OUT_DIR.mkdir(parents=True, exist_ok=True)
OUT_CLASSIFICATION_REPORT = OUT_DIR / f"classification_report_{HANDLE}.txt"
OUT_CONFUSION_RF = OUT_DIR / f"confusion_rf_{HANDLE}.png"
OUT_FEATURE_IMPORTANCE = OUT_DIR / f"feature_importance_{HANDLE}.csv"
OUT_CLUSTER_CROSSTAB = OUT_DIR / f"cluster_crosstab_{HANDLE}.csv"
OUT_SUP_VS_UNSUP_SCATTER = OUT_DIR / f"sup_vs_unsup_scatter_{HANDLE}.png"
OUT_LOGREG_REPORT = OUT_DIR / f"rf_vs_logreg_report_{HANDLE}.txt"
OUT_PCA_BONUS = OUT_DIR / f"pca_bonus_variance_{HANDLE}.png"
OUT_SEMI_SUPERVISED_REPORT = OUT_DIR / f"semi_supervised_report_{HANDLE}.txt"

RANDOM_STATE = 42


def load_and_prepare_data(data_path: Path) -> tuple:
    """
    Task 1: Pregatirea datelor
    - Incarca datele
    - Separa X (gene) si y (Label)
    - Encodeaza etichetele cu LabelEncoder
    - Imparte in train/test stratificat
    """
    
    df = pd.read_csv(data_path)
    print(f" Date incarcate: {df.shape[0]} probe x {df.shape[1]} coloane")
    print(f" Coloane: {list(df.columns[:5])} ... {list(df.columns[-3:])}")
    
    # Separare X si y
    # Prima coloana este sample_id, ultima este Label
    X = df.iloc[:, 1:-1]  # Genele (excludem sample_id si Label)
    y = df["Label"]
    
    print(f"\n Dimensiuni X (features): {X.shape}")
    print(f" Distributia claselor:")
    print(y.value_counts())
    
    # Encodare etichete
    le = LabelEncoder()
    y_encoded = le.fit_transform(y)
    classes = le.classes_
    print(f"\n Clase encodate: {dict(zip(classes, range(len(classes))))}")
    
    # Train/test split
    X_train, X_test, y_train, y_test = train_test_split(
        X, y_encoded,
        test_size=0.2,
        random_state=RANDOM_STATE,
        stratify=y_encoded
    )
    
    print(f"\n Train set: {X_train.shape[0]} probe")
    print(f" Test set: {X_test.shape[0]} probe")
    
    return X, y_encoded, X_train, X_test, y_train, y_test, le, classes


def supervised_random_forest(
    X_train: pd.DataFrame,
    X_test: pd.DataFrame,
    y_train: np.ndarray,
    y_test: np.ndarray,
    classes: np.ndarray,
    feature_names: list
) -> RandomForestClassifier:
    """
    Task 2: Supervised ML - Random Forest
    - Antrenare model
    - Generare classification_report
    - Matrice de confuzie + salvare PNG
    - Feature importances + salvare CSV
    """
    print("\n" + "=" * 70)
    print("TASK 2: SUPERVISED ML - RANDOM FOREST")
    print("=" * 70)
    
    # Antrenare Random Forest
    rf = RandomForestClassifier(
        n_estimators=200,
        max_depth=10,
        min_samples_split=5,
        random_state=RANDOM_STATE,
        n_jobs=-1
    )
    rf.fit(X_train, y_train)
    
    # Predictii
    y_pred = rf.predict(X_test)
    
    # Classification report
    report = classification_report(y_test, y_pred, target_names=classes)
    print(report)
    
    with open(OUT_CLASSIFICATION_REPORT, "w") as f:
        f.write("=" * 60 + "\n")
        f.write("CLASSIFICATION REPORT - RANDOM FOREST\n")
        f.write("=" * 60 + "\n\n")
        f.write(report)
        f.write(f"\n\nAccuracy: {accuracy_score(y_test, y_pred):.4f}")
        f.write(f"\nF1-Score (weighted): {f1_score(y_test, y_pred, average='weighted'):.4f}")
    print(f"[SAVED] {OUT_CLASSIFICATION_REPORT}")

    cm = confusion_matrix(y_test, y_pred)
    plt.figure(figsize=(8, 6))
    sns.heatmap(
        cm,
        annot=True,
        fmt="d",
        cmap="Blues",
        xticklabels=classes,
        yticklabels=classes,
        annot_kws={"size": 14}
    )
    plt.xlabel("Predicted", fontsize=12)
    plt.ylabel("Actual", fontsize=12)
    plt.title("Confusion Matrix - Random Forest", fontsize=14)
    plt.tight_layout()
    plt.savefig(OUT_CONFUSION_RF, dpi=300)
    plt.close()
    print(f"[SAVED] {OUT_CONFUSION_RF}")
    
    # Feature importances
    importances = rf.feature_importances_
    df_importance = pd.DataFrame({
        "Gene": feature_names,
        "Importance": importances
    }).sort_values("Importance", ascending=False)
    
    df_importance.to_csv(OUT_FEATURE_IMPORTANCE, index=False)
    print(f"[SAVED] {OUT_FEATURE_IMPORTANCE}")
    
    print("\n=== Top 10 Gene dupa Importanta ===")
    print(df_importance.head(10).to_string(index=False))
    
    return rf


def logistic_regression_comparison(
    X_train: pd.DataFrame,
    X_test: pd.DataFrame,
    y_train: np.ndarray,
    y_test: np.ndarray,
    classes: np.ndarray,
    rf_accuracy: float
) -> None:
    """
    Task 3 (Bonus +1p): Logistic Regression cu scaling
    - Comparatie cu Random Forest
    - Discutie modele liniare vs non-liniare
    """
    
    # Scaling pentru Logistic Regression
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)
    
    # Antrenare Logistic Regression
    logreg = LogisticRegression(
        max_iter=1000,
        random_state=RANDOM_STATE,
        solver="lbfgs"
    )
    logreg.fit(X_train_scaled, y_train)
    
    # Predictii
    y_pred_lr = logreg.predict(X_test_scaled)
    
    # Classification report
    report_lr = classification_report(y_test, y_pred_lr, target_names=classes)
    lr_accuracy = accuracy_score(y_test, y_pred_lr)
    lr_f1 = f1_score(y_test, y_pred_lr, average='weighted')

    
    print("\n=== COMPARATIE RF vs Logistic Regression ===")
    print(f"Random Forest Accuracy: {rf_accuracy:.4f}")
    print(f"Logistic Regression Accuracy: {lr_accuracy:.4f}")
    print(f"Diferenta: {rf_accuracy - lr_accuracy:.4f}")
    
    with open(OUT_LOGREG_REPORT, "w") as f:
        f.write("=" * 60 + "\n")
        f.write("COMPARATIE: RANDOM FOREST vs LOGISTIC REGRESSION\n")
        f.write("=" * 60 + "\n\n")
        f.write("LOGISTIC REGRESSION:\n")
        f.write(report_lr)
        f.write(f"\nAccuracy LR: {lr_accuracy:.4f}")
        f.write(f"\nF1-Score LR (weighted): {lr_f1:.4f}")
        f.write(f"\n\n{'='*60}\n")
        f.write("COMPARATIE SUMARA:\n")
        f.write(f"Random Forest Accuracy: {rf_accuracy:.4f}\n")
        f.write(f"Logistic Regression Accuracy: {lr_accuracy:.4f}\n")
        f.write(f"Diferenta (RF - LR): {rf_accuracy - lr_accuracy:.4f}\n")
        f.write("\n" + "=" * 60 + "\n")
    print(f"[SAVED] {OUT_LOGREG_REPORT}")


def unsupervised_pca_kmeans(
    X: pd.DataFrame,
    y_encoded: np.ndarray,
    classes: np.ndarray
) -> tuple:
    """
    Task 4: Unsupervised ML - PCA + KMeans
    - PCA pentru vizualizare
    - KMeans cu 2-4 clustere
    - Scatter PCA colorat dupa cluster
    - Crosstab: Label X Cluster
    """
    
    # Standardizare pentru PCA
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    
    # PCA
    pca = PCA(n_components=2, random_state=RANDOM_STATE)
    X_pca = pca.fit_transform(X_scaled)
    
    print(f" PCA - Varianta explicata: {pca.explained_variance_ratio_.sum():.2%}")
    print(f"  - PC1: {pca.explained_variance_ratio_[0]:.2%}")
    print(f"  - PC2: {pca.explained_variance_ratio_[1]:.2%}")
    
    # KMeans cu diferite numere de clustere
    n_clusters = 3  # Folosim 3 (numarul de clase)
    kmeans = KMeans(n_clusters=n_clusters, random_state=RANDOM_STATE, n_init="auto")
    clusters = kmeans.fit_predict(X_scaled)
    
    # Crosstab: Label × Cluster
    df_crosstab = pd.DataFrame({
        "TrueLabel": [classes[i] for i in y_encoded],
        "Cluster": clusters
    })
    crosstab = pd.crosstab(df_crosstab["TrueLabel"], df_crosstab["Cluster"])
    print("\n=== Crosstab: True Label × Cluster ===")
    print(crosstab)
    
    crosstab.to_csv(OUT_CLUSTER_CROSSTAB)
    
    # Vizualizare: Supervised vs Unsupervised scatter
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    # Plot 1: PCA colorat dupa eticheta reala (Supervised view)
    colors_sup = ["#E74C3C", "#27AE60", "#3498DB"]
    for i, cls in enumerate(classes):
        mask = y_encoded == i
        axes[0].scatter(
            X_pca[mask, 0],
            X_pca[mask, 1],
            c=colors_sup[i],
            label=cls,
            alpha=0.7,
            edgecolor="k",
            s=60
        )
    axes[0].set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]:.1%})", fontsize=11)
    axes[0].set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]:.1%})", fontsize=11)
    axes[0].set_title("PCA - Colorate dupa Eticheta Reala (Supervised)", fontsize=12)
    axes[0].legend(title="Label")
    axes[0].grid(alpha=0.3)
    
    # Plot 2: PCA colorat dupa cluster (Unsupervised view)
    colors_unsup = ["#9B59B6", "#F39C12", "#1ABC9C"]
    for cl in range(n_clusters):
        mask = clusters == cl
        axes[1].scatter(
            X_pca[mask, 0],
            X_pca[mask, 1],
            c=colors_unsup[cl],
            label=f"Cluster {cl}",
            alpha=0.7,
            edgecolor="k",
            s=60
        )
    axes[1].set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]:.1%})", fontsize=11)
    axes[1].set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]:.1%})", fontsize=11)
    axes[1].set_title("PCA - Colorate dupa Cluster KMeans (Unsupervised)", fontsize=12)
    axes[1].legend(title="Cluster")
    axes[1].grid(alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(OUT_SUP_VS_UNSUP_SCATTER, dpi=300)
    plt.close()
    print(f"[SAVED] {OUT_SUP_VS_UNSUP_SCATTER}")
    
    return X_pca, pca, clusters


def semi_supervised_experiment(
    X: pd.DataFrame,
    y_encoded: np.ndarray,
    classes: np.ndarray,
    unlabeled_fraction: float = 0.4
) -> None:
    """
    Task 5: Semi-Supervised Learning (mini-experiment)
    - Marcheaza 30-50% din etichete ca 'unknown'
    - Antreneaza RF doar pe date etichetate
    - Genereaza pseudo-etichete pentru datele neetichetate
    - Reantreneaza pe setul complet
    """
    print("\n" + "=" * 70)
    print("TASK 5: SEMI-SUPERVISED LEARNING")
    print("=" * 70)
    
    np.random.seed(RANDOM_STATE)
    n_samples = len(y_encoded)
    n_unlabeled = int(n_samples * unlabeled_fraction)
    
    # Cream o copie a etichetelor
    y_semi = y_encoded.copy()
    
    # Alegem indici random pentru a fi "unlabeled"
    unlabeled_indices = np.random.choice(n_samples, size=n_unlabeled, replace=False)
    labeled_indices = np.array([i for i in range(n_samples) if i not in unlabeled_indices])
    
    print(f" Total probe: {n_samples}")
    print(f" Probe etichetate: {len(labeled_indices)} ({100-unlabeled_fraction*100:.0f}%)")
    print(f" Probe neetichetate: {n_unlabeled} ({unlabeled_fraction*100:.0f}%)")
    
    # Split pentru evaluare finala
    X_train_eval, X_test_eval, y_train_eval, y_test_eval = train_test_split(
        X, y_encoded,
        test_size=0.2,
        random_state=RANDOM_STATE,
        stratify=y_encoded
    )
    
    # Model 1: Doar pe date etichetate (baseline)
    X_labeled = X.iloc[labeled_indices]
    y_labeled = y_encoded[labeled_indices]
    
    rf_baseline = RandomForestClassifier(
        n_estimators=100,
        random_state=RANDOM_STATE,
        n_jobs=-1
    )
    rf_baseline.fit(X_labeled, y_labeled)
    y_pred_baseline = rf_baseline.predict(X_test_eval)
    acc_baseline = accuracy_score(y_test_eval, y_pred_baseline)
    
    print(f"\n[BASELINE] Accuracy (antrenat doar pe {100-unlabeled_fraction*100:.0f}% date): {acc_baseline:.4f}")
    
    # Generare pseudo-etichete
    X_unlabeled = X.iloc[unlabeled_indices]
    pseudo_labels = rf_baseline.predict(X_unlabeled)
    
    # Combinam datele etichetate cu cele pseudo-etichetate
    X_combined = pd.concat([X_labeled, X_unlabeled], axis=0)
    y_combined = np.concatenate([y_labeled, pseudo_labels])
    
    # Model 2: Antrenat pe setul combinat (pseudo-labeling)
    rf_pseudo = RandomForestClassifier(
        n_estimators=100,
        random_state=RANDOM_STATE,
        n_jobs=-1
    )
    rf_pseudo.fit(X_combined, y_combined)
    y_pred_pseudo = rf_pseudo.predict(X_test_eval)
    acc_pseudo = accuracy_score(y_test_eval, y_pred_pseudo)
    
    print(f"[PSEUDO-LABELING] Accuracy (dupa pseudo-labeling): {acc_pseudo:.4f}")
    print(f"[DIFERENTA] {acc_pseudo - acc_baseline:+.4f}")
    
    # Model 3: Referinta - antrenat pe toate datele etichetate (upper bound)
    rf_full = RandomForestClassifier(
        n_estimators=100,
        random_state=RANDOM_STATE,
        n_jobs=-1
    )
    rf_full.fit(X_train_eval, y_train_eval)
    y_pred_full = rf_full.predict(X_test_eval)
    acc_full = accuracy_score(y_test_eval, y_pred_full)
    
    print(f"[FULL DATA] Accuracy (100% date etichetate): {acc_full:.4f}")
    
    # Salvare raport
    with open(OUT_SEMI_SUPERVISED_REPORT, "w") as f:
        f.write("=" * 60 + "\n")
        f.write("SEMI-SUPERVISED LEARNING EXPERIMENT\n")
        f.write("=" * 60 + "\n\n")
        f.write(f"Procent date neetichetate: {unlabeled_fraction*100:.0f}%\n")
        f.write(f"Probe etichetate: {len(labeled_indices)}\n")
        f.write(f"Probe neetichetate: {n_unlabeled}\n\n")
        f.write("REZULTATE:\n")
        f.write("-" * 40 + "\n")
        f.write(f"Baseline (doar date etichetate): {acc_baseline:.4f}\n")
        f.write(f"Dupa Pseudo-labeling: {acc_pseudo:.4f}\n")
        f.write(f"Full data (upper bound): {acc_full:.4f}\n")
        f.write(f"\nImbunatatire pseudo-labeling: {acc_pseudo - acc_baseline:+.4f}\n")

    print(f"[SAVED] {OUT_SEMI_SUPERVISED_REPORT}")


def bonus_pca_variance_comparison(X: pd.DataFrame) -> None:
    """
    Bonus: Comparatie PCA inainte si dupa eliminarea a 10 gene cu varianta mica
    """
    
    # Calculam varianta pentru fiecare gena
    variances = X.var()
    variances_sorted = variances.sort_values()
    
    # Gene de eliminat (cele 10 cu varianta mica)
    genes_to_remove = variances_sorted.head(10).index.tolist()
    
    # Datele fara genele cu varianta mica
    X_reduced = X.drop(columns=genes_to_remove)
    
    print(f"\n Dimensiuni originale: {X.shape}")
    print(f" Dimensiuni dupa eliminare: {X_reduced.shape}")
    
    # Standardizare
    scaler = StandardScaler()
    X_scaled_full = scaler.fit_transform(X)
    X_scaled_reduced = scaler.fit_transform(X_reduced)
    
    # PCA pe ambele seturi
    pca_full = PCA(n_components=2, random_state=RANDOM_STATE)
    X_pca_full = pca_full.fit_transform(X_scaled_full)
    
    pca_reduced = PCA(n_components=2, random_state=RANDOM_STATE)
    X_pca_reduced = pca_reduced.fit_transform(X_scaled_reduced)
    
    # Vizualizare comparativa
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    # Plot 1: PCA cu toate genele
    axes[0].scatter(
        X_pca_full[:, 0],
        X_pca_full[:, 1],
        c='steelblue',
        alpha=0.6,
        edgecolor="k",
        s=50
    )
    axes[0].set_xlabel(f"PC1 ({pca_full.explained_variance_ratio_[0]:.1%})", fontsize=11)
    axes[0].set_ylabel(f"PC2 ({pca_full.explained_variance_ratio_[1]:.1%})", fontsize=11)
    axes[0].set_title(f"PCA - Toate genele ({X.shape[1]} gene)\n"
                      f"Varianta explicata totala: {pca_full.explained_variance_ratio_.sum():.1%}", 
                      fontsize=12)
    axes[0].grid(alpha=0.3)
    
    # Plot 2: PCA fara genele cu varianta mica
    axes[1].scatter(
        X_pca_reduced[:, 0],
        X_pca_reduced[:, 1],
        c='coral',
        alpha=0.6,
        edgecolor="k",
        s=50
    )
    axes[1].set_xlabel(f"PC1 ({pca_reduced.explained_variance_ratio_[0]:.1%})", fontsize=11)
    axes[1].set_ylabel(f"PC2 ({pca_reduced.explained_variance_ratio_[1]:.1%})", fontsize=11)
    axes[1].set_title(f"PCA - Fara 10 gene cu varianta mica ({X_reduced.shape[1]} gene)\n"
                      f"Varianta explicata totala: {pca_reduced.explained_variance_ratio_.sum():.1%}", 
                      fontsize=12)
    axes[1].grid(alpha=0.3)
    
    plt.suptitle("Comparatie PCA inainte si dupa eliminarea genelor cu varianta mica", 
                 fontsize=14, y=1.02)
    plt.tight_layout()
    plt.savefig(OUT_PCA_BONUS, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"[SAVED] {OUT_PCA_BONUS}")
    
    print(f"\n Varianta explicata (toate genele): {pca_full.explained_variance_ratio_.sum():.2%}")
    print(f" Varianta explicata (fara 10 gene): {pca_reduced.explained_variance_ratio_.sum():.2%}")


def main() -> None:
    """
    Functia principala care ruleaza toate task-urile.
    """
    
    # Task 1: Pregatirea datelor
    X, y_encoded, X_train, X_test, y_train, y_test, le, classes = load_and_prepare_data(DATA_PATH)
    
    # Task 2: Supervised ML - Random Forest
    rf_model = supervised_random_forest(
        X_train, X_test, y_train, y_test, classes, list(X.columns)
    )
    rf_accuracy = accuracy_score(y_test, rf_model.predict(X_test))
    
    # Task 3 (Bonus): Logistic Regression
    logistic_regression_comparison(
        X_train, X_test, y_train, y_test, classes, rf_accuracy
    )
    
    # Task 4: Unsupervised ML - PCA + KMeans
    X_pca, pca, clusters = unsupervised_pca_kmeans(X, y_encoded, classes)
    
    # Task 5: Semi-Supervised Learning
    semi_supervised_experiment(X, y_encoded, classes, unlabeled_fraction=0.4)
    
    # Bonus: Comparatie PCA
    bonus_pca_variance_comparison(X)


if __name__ == "__main__":
    main()
