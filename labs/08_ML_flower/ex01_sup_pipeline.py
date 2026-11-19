"""Student implementation: training-only pipeline and shared-fold comparison.

Copy into submissions/<handle>/ and preserve the original teaching skeleton.
--check-input validates the infrastructure before completing the assessed TODOs.
"""
from pathlib import Path
import argparse
import sys

ROOT = next(p for p in Path(__file__).resolve().parents if (p / "labs/08_ML_flower/ml_lab.py").is_file())
sys.path.insert(0, str(ROOT / "labs/08_ML_flower"))
from ml_lab import prepare, development_folds  # noqa: E402


def compare_candidates(X_train, y_train, folds):
    """Return (candidates, cv_scores, selected_name).

    TODO: create Dummy, a scaled Logistic Pipeline, and Forest using the
    assignment settings. Reuse folds, compare balanced accuracy and choose
    between Logistic and Forest before opening the final test scores.
    """
    raise NotImplementedError("Implement the shared-fold comparison in your copy.")


def evaluate_frozen_models(candidates, X_train, y_train, X_test, y_test, out):
    """TODO: fit on development, save test metrics and confusion matrices.

    Preserve sample IDs in local predictions. All generated outputs go in out.
    """
    raise NotImplementedError("Implement evaluation after freezing the model choice.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--handle", required=True)
    parser.add_argument("--check-input", action="store_true")
    args = parser.parse_args()
    X, y, train, test, out = prepare(args.handle)
    print(f"{X.shape[0]} samples, {X.shape[1]} genes; development={len(train)}, test={len(test)}")
    if not args.check_input:
        candidates, cv_scores, selected = compare_candidates(X.loc[train], y.loc[train], development_folds(X.loc[train], y.loc[train]))
        print(f"Chosen using development CV: {selected}")
        evaluate_frozen_models(candidates, X.loc[train], y.loc[train], X.loc[test], y.loc[test], out)
