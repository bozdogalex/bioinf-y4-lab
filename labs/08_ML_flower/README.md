# Lab 8: Predicting condition from gene expression

Build and diagnose a supervised classifier using the same documented simulation as Labs 6–7. Compare regularized logistic regression, Random Forest and a dummy baseline. The lecture develops the biological question and evaluation logic. The assignment asks you to implement and explain the workflow.

[Lecture presentation](../../docs/presentations/08_ML_flower.pptx) · [Assignment and rubric](assignment.md)

## Data and preparation

`data/sample/lab06/expression_matrix.csv` contains 120 synthetic genes × 48 independent samples. `sample_metadata.csv` contains 24 controls and 24 treated samples (`tratament` in the source). [Generation and provenance](../../data/sample/lab06/README.md). These are simulated continuous expression values, already on an analysis-ready scale. G001–G120 are synthetic identifiers. No GTEx, TP53 or clinical biomarker claims are attached to them. The old six-row tissue demo is not used for this lab.

`ml_lab.py` aligns metadata by sample ID and transposes the matrix to 48 samples × 120 features. IDs and condition labels never enter X. Its fixed stratified split reserves 12 samples for final testing and 36 for development. Three-fold stratified CV uses identical folds across models.

## Run in Codespaces

Create/open your Codespace from your private course repository and wait for its environment setup to finish. Use the repository-root terminal and its configured Python environment. The course pins scikit-learn 1.5.x. If working outside the provided environment, install `requirements.txt` first.

Replace `yourhandle` with your own handle. These commands run reference analyses and create local artifacts:

```bash
python labs/08_ML_flower/ml_lab.py prepare --handle yourhandle
python labs/08_ML_flower/demo02_expr.py --handle yourhandle
python labs/08_ML_flower/demo01_sup_unsup.py --handle yourhandle
```

Preparation writes `features.csv`, `labels.csv` and `split_manifest.json` under `data/work/yourhandle/lab08/`. The manifest records sample IDs, CV folds, source hashes and library versions. All scripts also work when launched by absolute path from another directory. Existing artifact names are replaced on a rerun. Keep separate handles for separate investigations.

The benchmark saves `metrics.csv`, `benchmark.json`, three confusion-matrix PNGs and development-only feature ranks. The exploratory demo saves a development-only PCA plot and cluster cross-tab. `demo03_logreg.py` is a compatibility entry point to the same model comparison. It intentionally reuses the same split rather than creating a second evaluation protocol.

## Track and inspect the comparison

The benchmark uses the existing course MLflow helpers and the local file-backed setup also documented in the master's BDHB lab. Default store: repository-root `mlruns/`, independent of the directory from which the Python script is launched. `MLFLOW_TRACKING_URI` takes precedence when explicitly configured. Otherwise `MLFLOW_TRACKING_DIR` selects a local directory, with relative paths resolved from the repository root.

Each benchmark execution creates a comparison run and three nested candidate runs in experiment **BIOINF-Y4 Lab08**. Parameters record the model settings and split design. Metrics distinguish development CV from final test evaluation, with fold histories retained. Artifacts include the split manifest, data hashes and versions, confusion matrices, predictions, feature ranks and the selected fitted model. `benchmark.json` includes the comparison and candidate run IDs. Rerunning creates new runs rather than replacing old tracking records.

After running the benchmark, open a second repository-root Codespaces terminal:

```bash
mlflow ui --backend-store-uri ./mlruns --host 0.0.0.0 --port 5000
```

Open port 5000 from the Ports panel and keep its visibility private. In **BIOINF-Y4 Lab08**, expand the comparison's candidate runs, compare `cv_balanced_accuracy` with the separately named `test_*` metrics, then inspect `data/split_manifest.json` on the parent and an `evaluation` artifact on a candidate. The UI reads previously logged records; it can be started after training. If using a configured store, point the UI at that same store. See the [existing MLflow guide](../../docs/mlflow.md) for background.

The winning model is stored as an MLflow model artifact. The fitted logistic candidate, if selected, includes its scaler. Run IDs identify records, while data hashes identify the exact source files. Keep `mlruns/` and generated artifacts outside your assessment commit.

Optional self-training reference, after completing the core comparison:

```bash
python labs/08_ML_flower/ml_lab.py selftrain --handle yourhandle
```

This hides half of the development labels, retains 18 labeled samples, accepts predictions above a fixed 0.9 confidence threshold, and compares against a model using the same 18 observed labels. Both models use the same untouched test set. A confidence score is model-dependent and can be overconfident. Report improvements, ties and declines as observed.

## Your implementation

Copy the two exercise skeletons into `labs/08_ML_flower/submissions/yourhandle/`. Preserve the originals. They locate the shared input helper from their copied location. Check setup before filling in the assessed TODOs:

```bash
mkdir -p labs/08_ML_flower/submissions/yourhandle
cp labs/08_ML_flower/ex0*.py labs/08_ML_flower/submissions/yourhandle/
python labs/08_ML_flower/submissions/yourhandle/ex01_sup_pipeline.py --handle yourhandle --check-input
```

A normal skeleton run raises `NotImplementedError` until you complete the modeling functions. That is intentional. Read [assignment.md](assignment.md) for the required analyses and rubric.

The runnable demos provide a baseline, as in Labs 6–7. They are distinct from the assessed skeletons. Exercise 1 leaves model comparison, final evaluation and candidate-run logging to you. Exercise 2 leaves a training-only 20-gene signature comparison, label-permutation control, expression exploration, condition-specific gene profiles and diagnostic MLflow logging to you. The supplied setup handles paths, validated input, deterministic splits and the diagnostic run context. Complete the functions in your own copies and explain the biological meaning and limits of the outputs.

## Evaluation contract

Scale inside the logistic Pipeline, so CV learns the scaler using only its fit fold. Fit any learned feature selection, PCA or module representation inside CV if used predictively. The exploratory PCA/KMeans demo uses development samples only and is not fed into the benchmark. Choose the classifier using development CV, then evaluate the frozen comparison once on the held-out test. Repeated test-driven changes turn that set into development data.

The simulation has a strong designed response. Perfect test performance is possible with just six test examples per class. Rankings identify model dependence in this simulation. For real cohorts, define the endpoint, normalization, donor grouping and independent validation before applying this workflow. Repeated samples from one donor belong in the same split.

## Submission

Use your private repository shared with `bozdogalex`. Submit your PR URL and deadline commit SHA through the LMS, following the [course guide](../../docs/git-workflow.md). Commit completed code, a concise Markdown report and at most ten rows of summary metadata. Keep generated matrices, full result CSVs and PNGs in ignored `data/work/yourhandle/lab08/`. Include exact reproduction commands and AI-assistance attribution.

## Teaching-material provenance

The 2026 Lab 8 redesign, reference workflow and assignment were prepared with OpenAI Codex assistance. Protocols and calculations were checked against the course simulation and scikit-learn 1.5 documentation, and the runnable workflow was tested locally.
