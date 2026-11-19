# Lab 8 assignment: a defensible expression classifier

Work in your private course repository shared with `bozdogalex`. Submit the assessment PR URL and deadline commit SHA via the LMS. Preserve the teaching skeletons and put completed copies in `labs/08_ML_flower/submissions/<handle>/`.

## 1. Establish the input and split (2 points)

Use `ml_lab.prepare(handle)` to obtain X, y, development IDs, test IDs and the local output folder. Explain the 48 × 120 data layout, sample-ID alignment and synthetic provenance. Verify the 36/12 split is disjoint, contains both classes, and has 18/18 development and 6/6 test labels. Save the manifest locally. A handle selects an output folder, never a predictor.

## 2. Compare three prespecified candidates (3 points)

Implement `compare_candidates` in your copy of `ex01_sup_pipeline.py`:

- `DummyClassifier(strategy="most_frequent")`.
- `Pipeline(StandardScaler(), LogisticRegression(C=1.0, max_iter=2000, random_state=42))`. Use `make_pipeline` or an equivalent named Pipeline.
- `RandomForestClassifier(n_estimators=200, min_samples_leaf=2, max_features="sqrt", random_state=42, n_jobs=1)`.

Use `development_folds` to reuse the same three stratified folds. Compare balanced accuracy, retaining the three fold values and their mean. Choose between logistic regression and Forest by mean CV balanced accuracy, with logistic regression winning a tie. Dummy is a reference. Return the candidates, scores and chosen name. Explain why a small difference across these folds would be weak evidence for choosing a more complex model. Additional tuning is outside the required comparison.

## 3. Freeze and evaluate (2 points)

Record the chosen classifier before evaluating the final test set. Implement `evaluate_frozen_models`: fit all three prespecified candidates on the 36 development samples, obtain test predictions on the same 12 IDs, save confusion matrices and report accuracy, balanced accuracy and macro-F1. Give class-specific recall and the number of errors. Explain what one extra mistake does to accuracy with n=12. Do not revise the pipeline using the test results.

## 4. Diagnose the result (3 points)

Implement `diagnose` in your copy of `ex02_logreg_vs_ref.py`, using development data only:

1. Permute training labels once with `np.random.default_rng(42).permutation(...)`. Rerun logistic CV with scaling inside its Pipeline and three-fold stratification of the permuted labels. Report the three scores. This is a negative control, not a permutation p-value and not a reason to select a favorable random seed.
2. Standardize development expression, fit PCA with two components for a plot, and fit `KMeans(n_clusters=2, n_init=10, random_state=42)` on all 120 scaled genes. Compare clusters with measured labels using a cross-tab. Explain why a cluster number cannot be read directly as a condition and why visible separation does not replace predictive evaluation.
3. Fit the prespecified Forest on development data and rank `feature_importances_`. Discuss two top-ranked synthetic features and the implications of correlated genes. Explain what additional evidence a real causal or biomarker claim would require. Feature rankings from test-driven selection receive no credit.

Full credit depends on a correct protocol and a supported interpretation, including ties or surprising results. A higher test score does not earn extra points. Required total: **10 points**.

## Optional extension (+1, subject to the course grade cap)

Mask 50% of development labels with a stratified seed-42 split. Train the initial Forest on the 18 observed labels. Apply self-training to the development pool only using threshold 0.9 and at most ten iterations. Keep the same 12 test IDs outside both fits. Record the labeled/masked IDs, number of accepted pseudo-labels and both test scores. Discuss how an incorrect confident prediction could be reinforced. An improvement is not required. The supplied `selftrain` command demonstrates the protocol.

## Deliverables

Commit completed Python copies, a Markdown report of about two pages, and a small summary table (at most ten rows). The report should identify the data, split, fold scores, frozen model choice, final errors, diagnostic results and limitations. Include commands, source hashes/library versions from the local manifest, and AI-assistance attribution when applicable.

Keep all generated PNGs, matrices, full prediction and importance tables, and manifests under ignored `data/work/<handle>/lab08/`. Refer to their reproduction commands and filenames in the report. Submit no clinical/personal data. Follow the [submission guide](../../docs/git-workflow.md) and [data policy](../../docs/GDPR_and_DataPolicy.md).
