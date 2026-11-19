"""Lab 8 reference workflow. Generated artifacts stay in ignored data/work."""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sklearn
from sklearn.base import clone
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.dummy import DummyClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, balanced_accuracy_score, classification_report, f1_score, ConfusionMatrixDisplay
from sklearn.model_selection import StratifiedKFold, cross_val_score, train_test_split
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.semi_supervised import SelfTrainingClassifier

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "data/sample/lab06"
SEED = 42
LABELS = ["Control", "Treated"]


def load_dataset(source: Path = SOURCE):
    """Align by sample ID and transpose genes x samples. Exclude identifiers.

    Input is continuous simulated expression on an analysis-ready scale.
    Raw RNA-seq counts would need a different preparation.
    """
    expr = pd.read_csv(source / "expression_matrix.csv", index_col=0)
    meta = pd.read_csv(source / "sample_metadata.csv")
    if expr.index.has_duplicates or expr.columns.has_duplicates:
        raise ValueError("Expression gene and sample identifiers must be unique.")
    if not {"Sample", "Condition"}.issubset(meta.columns):
        raise ValueError("Metadata needs Sample and Condition columns.")
    if meta.Sample.isna().any() or meta.Sample.duplicated().any():
        raise ValueError("Metadata sample identifiers must be present and unique.")
    meta = meta.set_index("Sample")
    if set(expr.columns) != set(meta.index):
        raise ValueError("Expression and metadata must contain exactly the same sample IDs.")
    X = expr.T.apply(pd.to_numeric, errors="raise").astype(float)
    if X.shape[1] == 0 or not np.isfinite(X.to_numpy()).all():
        raise ValueError("Expression must have finite numeric features.")
    y = meta.loc[X.index, "Condition"].map({"control": "Control", "tratament": "Treated"})
    if y.isna().any() or set(y) != set(LABELS) or y.value_counts().min() < 8:
        raise ValueError("Expected control/tratament labels with at least 8 samples per class.")
    y.name = "Condition"
    return X, y


def split_ids(X, y):
    return train_test_split(X.index.to_numpy(), test_size=0.25, random_state=SEED, stratify=y)


def development_folds(X, y):
    return list(StratifiedKFold(3, shuffle=True, random_state=SEED).split(X, y))


def models():
    return {
        "Dummy": DummyClassifier(strategy="most_frequent"),
        "Logistic": make_pipeline(StandardScaler(), LogisticRegression(C=1.0, max_iter=2000, random_state=SEED)),
        "Forest": RandomForestClassifier(n_estimators=200, min_samples_leaf=2, max_features="sqrt", random_state=SEED, n_jobs=1),
    }


def scores(y, pred):
    return {"accuracy": float(accuracy_score(y, pred)), "balanced_accuracy": float(balanced_accuracy_score(y, pred)),
            "macro_f1": float(f1_score(y, pred, labels=LABELS, average="macro", zero_division=0))}


def output_dir(handle):
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]{0,38}", handle):
        raise ValueError("Use a handle containing letters, numbers, underscores or hyphens.")
    out = ROOT / "data/work" / handle / "lab08"
    out.mkdir(parents=True, exist_ok=True)
    return out


def write_json(path, obj):
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def prepare(handle):
    X, y = load_dataset()
    train, test = split_ids(X, y)
    out = output_dir(handle)
    X.to_csv(out / "features.csv", index_label="Sample")
    y.to_csv(out / "labels.csv", index_label="Sample")
    folds = development_folds(X.loc[train], y.loc[train])
    write_json(out / "split_manifest.json", {
        "dataset": "Lab 6 simulation, 120 synthetic genes and 48 independent samples",
        "seed": SEED, "test_fraction": 0.25, "train_ids": train.tolist(), "test_ids": test.tolist(),
        "cv_folds": [{"fit_ids": train[a].tolist(), "validation_ids": train[b].tolist()} for a, b in folds],
        "label_mapping": {"control": "Control", "tratament": "Treated"},
        "source_sha256": {n: hashlib.sha256((SOURCE / n).read_bytes()).hexdigest() for n in ["expression_matrix.csv", "sample_metadata.csv"]},
        "versions": {"python": platform.python_version(), "numpy": np.__version__, "pandas": pd.__version__, "sklearn": sklearn.__version__},
    })
    return X, y, train, test, out


def benchmark(handle):
    X, y, train, test, out = prepare(handle)
    Xt, yt = X.loc[train], y.loc[train]
    folds = development_folds(Xt, yt)
    candidates = models()
    cv = {n: cross_val_score(m, Xt, yt, cv=folds, scoring="balanced_accuracy") for n, m in candidates.items()}
    chosen = max(["Logistic", "Forest"], key=lambda n: cv[n].mean())
    rows, reports, fitted = [], {}, {}
    for name, model in candidates.items():
        model.fit(Xt, yt)
        fitted[name] = model
        pred = model.predict(X.loc[test])
        rows.append({"model": name, "cv_balanced_accuracy": float(cv[name].mean()), "cv_sd": float(cv[name].std(ddof=1)), **scores(y.loc[test], pred)})
        reports[name] = classification_report(y.loc[test], pred, labels=LABELS, output_dict=True, zero_division=0)
        ConfusionMatrixDisplay.from_predictions(y.loc[test], pred, labels=LABELS, colorbar=False, cmap="Blues")
        plt.title(f"{name}: fixed test, n={len(test)}")
        plt.tight_layout()
        plt.savefig(out / f"confusion_{name.lower()}.png", dpi=180)
        plt.close()
    # Negative control uses development labels only, with fresh stratified folds.
    permuted = pd.Series(np.random.default_rng(SEED).permutation(yt.to_numpy()), index=yt.index)
    control = cross_val_score(models()["Logistic"], Xt, permuted, cv=development_folds(Xt, permuted), scoring="balanced_accuracy")
    pd.DataFrame({"Gene": X.columns, "ImpurityImportance": fitted["Forest"].feature_importances_}).sort_values("ImpurityImportance", ascending=False).to_csv(out / "development_forest_importance.csv", index=False)
    coef = fitted["Logistic"].named_steps["logisticregression"].coef_[0]
    pd.DataFrame({"Gene": X.columns, "CoefficientForTreated": coef}).to_csv(out / "development_logistic_coefficients.csv", index=False)
    result = {
        "selected_on_cv": chosen, "selection_rule": "Maximum CV balanced accuracy; Logistic wins ties; Dummy is a reference",
        "metrics": rows, "cv_folds": {n: v.tolist() for n, v in cv.items()},
        "permuted_training_labels_cv": control.tolist(), "classification_reports": reports,
        "test_ids": test.tolist(), "test_target": y.loc[test].tolist(),
        "test_predictions": {n: m.predict(X.loc[test]).tolist() for n, m in fitted.items()},
    }
    write_json(out / "benchmark.json", result)
    pd.DataFrame(rows).to_csv(out / "metrics.csv", index=False)
    print(pd.DataFrame(rows).to_string(index=False))
    print(f"Selected before test: {chosen}. Artifacts: {out}")
    return result


def explore(handle):
    X, y, train, _, out = prepare(handle)
    # Separate development-only view. Any predictive PCA must live inside CV.
    Z = StandardScaler().fit_transform(X.loc[train])
    pca = PCA(n_components=2)
    points = pca.fit_transform(Z)
    clusters = KMeans(n_clusters=2, n_init=10, random_state=SEED).fit_predict(Z)
    ctab = pd.crosstab(y.loc[train], clusters, colnames=["Cluster"])
    ctab.to_csv(out / "development_cluster_crosstab.csv")
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    for label, color in zip(LABELS, ["#315D9B", "#AF7B51"]):
        take = y.loc[train].to_numpy() == label
        axes[0].scatter(points[take, 0], points[take, 1], label=label, color=color)
    axes[1].scatter(points[:, 0], points[:, 1], c=clusters, cmap="coolwarm")
    axes[0].legend()
    for ax, title in zip(axes, ["Measured condition", "KMeans on all 120 scaled genes"]):
        ax.set(xlabel=f"PC1 ({pca.explained_variance_ratio_[0]:.1%})", ylabel=f"PC2 ({pca.explained_variance_ratio_[1]:.1%})", title=title)
    fig.tight_layout()
    fig.savefig(out / "development_pca.png", dpi=180)
    plt.close(fig)
    print(ctab.to_string())
    return points, clusters, pca.explained_variance_ratio_


def selftrain(handle):
    X, y, train, test, out = prepare(handle)
    labeled, unlabeled = train_test_split(train, test_size=0.5, random_state=SEED, stratify=y.loc[train])
    masked = y.loc[train].map({"Control": 0, "Treated": 1}).to_numpy()
    masked[np.isin(train, unlabeled)] = -1
    base = models()["Forest"]
    supervised = clone(base).fit(X.loc[labeled], y.loc[labeled].map({"Control": 0, "Treated": 1}))
    # base_estimator is the constructor name in the course-pinned sklearn 1.5.
    semi = SelfTrainingClassifier(base_estimator=clone(base), threshold=0.9, max_iter=10).fit(X.loc[train], masked)
    test_y = y.loc[test].map({"Control": 0, "Treated": 1})
    result = {
        "labeled_ids": labeled.tolist(), "masked_ids": unlabeled.tolist(), "test_ids": test.tolist(),
        "initial_balanced_accuracy": float(balanced_accuracy_score(test_y, supervised.predict(X.loc[test]))),
        "selftrained_balanced_accuracy": float(balanced_accuracy_score(test_y, semi.predict(X.loc[test]))),
        "accepted_pseudolabels": int(np.sum(semi.labeled_iter_ > 0)), "threshold": 0.9,
    }
    write_json(out / "selftraining.json", result)
    print(json.dumps(result, indent=2))
    return result


def main(default_command=None):
    parser = argparse.ArgumentParser(description=__doc__)
    if default_command is None:
        parser.add_argument("command", choices=["prepare", "benchmark", "explore", "selftrain"])
    parser.add_argument("--handle", default="demo")
    args = parser.parse_args()
    command = default_command or args.command
    {"prepare": prepare, "benchmark": benchmark, "explore": explore, "selftrain": selftrain}[command](args.handle)


if __name__ == "__main__":
    main()
