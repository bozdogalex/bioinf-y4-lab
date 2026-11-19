"""Student skeleton: expression signatures, negative controls and interpretation.

The CLI, validated input and development/test boundary are supplied. Complete
the expression-analysis functions and MLflow logging in your copied skeleton.
"""
from pathlib import Path
import argparse
import sys

ROOT = next(p for p in Path(__file__).resolve().parents if (p / "labs/08_ML_flower/ml_lab.py").is_file())
sys.path.insert(0, str(ROOT / "labs/08_ML_flower"))
sys.path.insert(0, str(ROOT))
import mlflow  # noqa: E402, F401
from ml_lab import prepare, development_folds, configure_tracking, EXPERIMENT, SEED  # noqa: E402
from mlops.mlflow_utils import start_run, log_params, log_metrics  # noqa: E402, F401


def compare_gene_panels(X_train, y_train, folds):
    """TODO: compare all genes with a 20-gene expression signature using CV.

    Return a table with columns panel, fold, balanced_accuracy (six rows).
    Reuse the same three development folds for both logistic pipelines.
    The 20-gene pipeline must fit SelectKBest(f_classif, k=20) and scaling
    INSIDE each fit fold. The full-gene pipeline uses scaling and logistic
    regression with the same settings as Exercise 1. The driver saves the table.
    """
    raise NotImplementedError("Compare expression signatures using development folds.")


def negative_control(X_train, y_train):
    """TODO: permute labels with seed 42; return three logistic CV scores.

    Recreate stratified folds for the permuted labels. Scaling belongs inside
    the pipeline. Use all genes and preserve the untouched test samples.
    """
    raise NotImplementedError("Implement the training-label negative control.")


def explore_expression(X_train, y_train, out):
    """TODO: scale expression, fit PCA/KMeans and return a label/cluster table.

    Use two PCs for the plot and all scaled genes for KMeans. Save the plot
    and cross-tab under out. Labels are used after fitting for interpretation.
    """
    raise NotImplementedError("Investigate expression structure in development samples.")


def profile_predictive_genes(X_train, y_train, out):
    """TODO: return (top_two_gene_table, pearson_correlation_between_them).

    Fit the prescribed Forest on development samples. For its top two genes,
    calculate mean expression in each condition and Treated minus Control.
    Table columns: Gene, ImpurityImportance, MeanControl, MeanTreated,
    MeanDifference. Compute their Pearson correlation on development samples.
    Save the table locally. These prepared simulation values do not establish
    a log2 fold change, real gene function or a causal treatment mechanism.
    """
    raise NotImplementedError("Connect predictive gene ranks to expression profiles.")


def log_diagnostics(panel_scores, null_scores, clusters, gene_profiles, correlation, out):
    """TODO: log summaries and saved analysis files to the active MLflow run.

    Log each panel's mean CV balanced accuracy, the negative-control mean,
    and top-pair correlation. Log the panel table, PCA plot, cluster cross-tab
    and gene-profile table as artifacts. Avoid final-test metrics in this run:
    every input to these diagnostics comes from development data only.
    """
    raise NotImplementedError("Record the bioinformatics diagnostics with MLflow.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--handle", required=True)
    args = parser.parse_args()
    X, y, train, _, out = prepare(args.handle)
    configure_tracking()
    Xt, yt = X.loc[train], y.loc[train]
    with start_run(EXPERIMENT, run_name=f"expression-diagnostics-{args.handle}",
                   tags={"analysis": "development-only", "handle": args.handle}):
        log_params({"seed": SEED, "n_development": len(train),
                    "n_genes": Xt.shape[1], "signature_size": 20})
        panel_scores = compare_gene_panels(Xt, yt, development_folds(Xt, yt))
        panel_scores.to_csv(out / "development_gene_panels.csv", index=False)
        null_scores = negative_control(Xt, yt)
        clusters = explore_expression(Xt, yt, out)
        gene_profiles, correlation = profile_predictive_genes(Xt, yt, out)
        log_diagnostics(panel_scores, null_scores, clusters, gene_profiles, correlation, out)
        mlflow.log_artifact(str(out / "split_manifest.json"), "data")
