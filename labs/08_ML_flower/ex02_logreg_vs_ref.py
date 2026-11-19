"""Student extension: diagnosis of an apparently successful classifier."""
from pathlib import Path
import argparse
import sys

ROOT = next(p for p in Path(__file__).resolve().parents if (p / "labs/08_ML_flower/ml_lab.py").is_file())
sys.path.insert(0, str(ROOT / "labs/08_ML_flower"))
from ml_lab import prepare  # noqa: E402


def diagnose(X_train, y_train, out):
    """TODO: training-label permutation, PCA/KMeans and development feature ranks.

    Keep scaling inside CV for the permutation control. Fit exploratory PCA
    and KMeans on development data only. Use all scaled genes for KMeans,
    two PCs for the plot. Explain what feature rankings support biologically.
    """
    raise NotImplementedError("Implement the diagnostic analyses in your copy.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--handle", required=True)
    args = parser.parse_args()
    X, y, train, _, out = prepare(args.handle)
    diagnose(X.loc[train], y.loc[train], out)
