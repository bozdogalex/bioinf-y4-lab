"""Regressions for Lab 8's data alignment and evaluation boundary."""
import importlib.util
import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd
from sklearn.model_selection import cross_validate

PATH = Path(__file__).resolve().parents[1] / "labs/08_ML_flower/ml_lab.py"
SPEC = importlib.util.spec_from_file_location("lab08", PATH)
lab = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(lab)


class Lab08ProtocolTests(unittest.TestCase):
    def test_comparison_is_reproducible_from_mlflow_runs(self):
        import mlflow
        from mlflow.tracking import MlflowClient

        previous_uri = mlflow.get_tracking_uri()
        try:
            with tempfile.TemporaryDirectory() as folder:
                base = Path(folder)
                outputs = base / "outputs"
                outputs.mkdir()
                uri = (base / "tracking").as_uri()
                with patch.dict(os.environ, {"MLFLOW_TRACKING_URI": uri}):
                    with patch.object(lab, "output_dir", return_value=outputs), patch("builtins.print"):
                        result = lab.benchmark("qa08")
                    client = MlflowClient(tracking_uri=uri)
                    parent = client.get_run(result["tracking"]["comparison_run_id"])
                    children = client.search_runs(
                        [parent.info.experiment_id],
                        filter_string=f"tags.mlflow.parentRunId = '{parent.info.run_id}'",
                    )
                    self.assertEqual(len(children), 3)
                    self.assertEqual(parent.data.tags["selected_on_cv"], "Forest")
                    self.assertEqual(parent.data.params["n_development"], "36")
                    downloaded = Path(client.download_artifacts(parent.info.run_id, "data/split_manifest.json"))
                    manifest = json.loads(downloaded.read_text(encoding="utf-8"))
                    self.assertIn("source_sha256", manifest)
                    self.assertEqual(manifest["test_ids"], result["test_ids"])
                    for child in children:
                        self.assertEqual(child.info.status, "FINISHED")
                        self.assertIn("cv_balanced_accuracy", child.data.metrics)
                        self.assertIn("test_balanced_accuracy", child.data.metrics)
                        history = client.get_metric_history(child.info.run_id, "cv_fold_balanced_accuracy")
                        self.assertEqual([m.step for m in history], [1, 2, 3])
                        self.assertTrue(client.list_artifacts(child.info.run_id, "evaluation"))
                    winning_id = result["tracking"]["model_run_ids"]["Forest"]
                    restored = mlflow.sklearn.load_model(f"runs:/{winning_id}/model")
                    X, _ = lab.load_dataset()
                    np.testing.assert_array_equal(
                        restored.predict(X.loc[result["test_ids"]]),
                        result["test_predictions"]["Forest"],
                    )
        finally:
            mlflow.set_tracking_uri(previous_uri)

    def test_selftraining_keeps_test_outside_both_training_pools(self):
        with tempfile.TemporaryDirectory() as folder:
            with patch.object(lab, "output_dir", return_value=Path(folder)), patch("builtins.print"):
                result = lab.selftrain("qa08")
        labeled = set(result["labeled_ids"])
        masked = set(result["masked_ids"])
        test = set(result["test_ids"])
        self.assertEqual((len(labeled), len(masked), len(test)), (18, 18, 12))
        self.assertFalse(labeled & masked)
        self.assertFalse((labeled | masked) & test)
        self.assertEqual(len(labeled | masked | test), 48)
        self.assertLessEqual(result["accepted_pseudolabels"], len(masked))

    def test_metadata_order_and_mismatched_ids(self):
        X, y = lab.load_dataset()
        with tempfile.TemporaryDirectory() as folder:
            source = Path(folder)
            X.T.to_csv(source / "expression_matrix.csv", index_label="Gene")
            meta = pd.read_csv(lab.SOURCE / "sample_metadata.csv").iloc[::-1]
            meta.to_csv(source / "sample_metadata.csv", index=False)
            X2, y2 = lab.load_dataset(source)
            pd.testing.assert_frame_equal(X, X2)
            pd.testing.assert_series_equal(y, y2)
            meta.iloc[0, 0] = "not-a-sample"
            meta.to_csv(source / "sample_metadata.csv", index=False)
            with self.assertRaisesRegex(ValueError, "exactly the same sample"):
                lab.load_dataset(source)

    def test_cv_scaler_sees_only_its_fit_fold(self):
        X, y = lab.load_dataset()
        train, test = lab.split_ids(X, y)
        self.assertFalse(set(train) & set(test))
        self.assertEqual((len(train), len(test)), (36, 12))
        self.assertEqual(y.loc[test].value_counts().to_dict(), {"Control": 6, "Treated": 6})
        Xt, yt = X.loc[train], y.loc[train]
        folds = lab.development_folds(Xt, yt)
        result = cross_validate(lab.models()["Logistic"], Xt, yt, cv=folds, return_estimator=True)
        validation_ids = []
        for (fit, valid), model in zip(folds, result["estimator"]):
            scaler = model.named_steps["standardscaler"]
            self.assertEqual(scaler.n_samples_seen_, 24)
            np.testing.assert_allclose(scaler.mean_, Xt.iloc[fit].mean().to_numpy())
            self.assertFalse(set(train[fit]) & set(train[valid]))
            self.assertFalse(set(test) & set(train[valid]))
            validation_ids.extend(train[valid])
        self.assertCountEqual(validation_ids, train)


if __name__ == "__main__":
    unittest.main()
