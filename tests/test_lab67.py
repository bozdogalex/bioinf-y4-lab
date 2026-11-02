"""Regression checks for the Lab 6 network and Lab 7 handoff."""
import importlib.util
from pathlib import Path
import unittest
import tempfile
import json

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def load_script(name, relative):
    spec = importlib.util.spec_from_file_location(name, ROOT / relative)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class CorrelationNetworkTests(unittest.TestCase):
    def test_existing_handoff_is_protected_before_export(self):
        lab6 = load_script("lab06", "labs/06_wgcna/ex01_gce_networks.py")
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "results"
            lab6.run_analysis(lab6.DEFAULT_INPUT, output, metadata_csv=lab6.DEFAULT_METADATA)
            original = {path.name: path.read_bytes() for path in output.iterdir()}
            # An identical rerun is reproducible and allowed.
            lab6.run_analysis(lab6.DEFAULT_INPUT, output, metadata_csv=lab6.DEFAULT_METADATA)
            self.assertEqual(original, {path.name: path.read_bytes() for path in output.iterdir()})
            changed_metadata = Path(temporary) / "metadata.csv"
            meta = pd.read_csv(lab6.DEFAULT_METADATA)
            meta.loc[0, "Condition"] = "changed"
            meta.to_csv(changed_metadata, index=False)
            changed_input = Path(temporary) / "expression.csv"
            expr = pd.read_csv(lab6.DEFAULT_INPUT, index_col=0)
            expr.iloc[0, 0] += 0.1
            expr.to_csv(changed_input)
            changes = [
                {"threshold": 0.5, "metadata_csv": lab6.DEFAULT_METADATA},
                {"input_scale": "normalized_counts", "metadata_csv": lab6.DEFAULT_METADATA},
                {"metadata_csv": changed_metadata},
                {"metadata_csv": None},
            ]
            for change in changes:
                with self.subTest(change=change):
                    with self.assertRaisesRegex(ValueError, "overwrite"):
                        lab6.run_analysis(lab6.DEFAULT_INPUT, output, **change)
                    self.assertEqual(original, {path.name: path.read_bytes() for path in output.iterdir()})
            with self.assertRaisesRegex(ValueError, "overwrite"):
                lab6.run_analysis(changed_input, output, metadata_csv=lab6.DEFAULT_METADATA)
            self.assertEqual(original, {path.name: path.read_bytes() for path in output.iterdir()})
            # Explicit replacement without metadata clears the old cohort file.
            report = lab6.run_analysis(lab6.DEFAULT_INPUT, output, overwrite=True)
            self.assertFalse((output / "sample_metadata.csv").exists())
            self.assertFalse(report["sample_metadata_included"])

    def test_lab7_rejects_graph_counts_from_a_different_run(self):
        lab6 = load_script("lab06", "labs/06_wgcna/ex01_gce_networks.py")
        lab7 = load_script("lab07", "labs/07_network_viz/ex01_network_viz.py")
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary)
            lab6.run_analysis(lab6.DEFAULT_INPUT, output)
            path = output / "parameters.json"
            params = json.loads(path.read_text(encoding="utf-8"))
            params["edges"] += 1
            path.write_text(json.dumps(params), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "parameters.json"):
                lab7.load_lab06_results(output)

    def test_prepared_scale_and_input_validation(self):
        lab = load_script("lab06", "labs/06_wgcna/ex01_gce_networks.py")
        expr = pd.DataFrame([[6, 7, 8], [2, 2, 2], [3, 4, 2]],
                            index=["variable", "constant", "other"])
        prepared = lab.log_and_filter(expr, variance_threshold=0)
        self.assertEqual(list(prepared.index), ["variable", "other"])
        self.assertEqual(prepared.loc["variable", 2], 8)
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "duplicate_samples.csv"
            path.write_text("Gene,s1,s1,s3\na,1,2,3\nb,4,5,6\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "probe"):
                lab.read_expression_matrix(path)
            path.write_text("Gene,s1,s2,s3\n   ,1,2,3\nb,4,5,6\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "gen"):
                lab.read_expression_matrix(path)

    def test_gene_rows_and_opposite_profiles(self):
        lab = load_script("lab06", "labs/06_wgcna/ex01_gce_networks.py")
        expr = pd.DataFrame(
            [[1, 2, 3, 4], [4, 3, 2, 1], [1, 1, 2, 2]],
            index=["up", "down", "step"],
            columns=["s1", "s2", "s3", "s4"],
        )
        corr = lab.correlation_matrix(expr, method="pearson")
        self.assertEqual(list(corr.index), ["up", "down", "step"])
        self.assertAlmostEqual(corr.loc["up", "down"], -1)
        positive = lab.adjacency_from_correlation(corr, 0.95, use_abs=False)
        unsigned = lab.adjacency_from_correlation(corr, 0.95, use_abs=True)
        self.assertEqual(positive.loc["up", "down"], 0)
        self.assertEqual(unsigned.loc["up", "down"], 1)
        np.testing.assert_array_equal(np.diag(unsigned), np.zeros(3))

    def test_saved_network_reaches_lab7_without_reconstruction(self):
        lab6 = load_script("lab06", "labs/06_wgcna/ex01_gce_networks.py")
        lab7 = load_script("lab07", "labs/07_network_viz/ex01_network_viz.py")
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary)
            report = lab6.run_analysis(lab6.DEFAULT_INPUT, output, metadata_csv=lab6.DEFAULT_METADATA)
            graph, modules, parameters = lab7.load_lab06_results(output)
            self.assertEqual(graph.number_of_nodes(), 120)
            self.assertEqual(graph.number_of_edges(), report["edges"])
            self.assertEqual(sum(value == -1 for value in modules.values()), report["isolates"])
            self.assertEqual(parameters["input_sha256"], report["input_sha256"])
            saved_edges = pd.read_csv(output / "edges.csv")
            edge = saved_edges.iloc[0]
            self.assertAlmostEqual(graph[edge.Source][edge.Target]["weight"], edge.Weight)
            metadata = pd.read_csv(output / "sample_metadata.csv")
            expr = pd.read_csv(output / "processed_expression.csv", index_col=0)
            self.assertEqual(list(metadata.Sample), list(expr.columns))
            # Missing module IDs must be an error, never a silent dropped gene.
            mapping = pd.read_csv(output / "modules.csv")
            mapping.iloc[1:].to_csv(output / "modules.csv", index=False)
            with self.assertRaisesRegex(ValueError, "gene"):
                lab7.load_lab06_results(output)

    def test_sparse_and_edgeless_networks_export_only_connected_hubs(self):
        lab6 = load_script("lab06", "labs/06_wgcna/ex01_gce_networks.py")
        lab7 = load_script("lab07", "labs/07_network_viz/ex01_network_viz.py")
        columns = ["Gene", "Module", "Degree", "Strength", "IntramodularStrength"]
        with tempfile.TemporaryDirectory() as temporary:
            base = Path(temporary)
            for threshold in [0.9, 1.0]:
                source = base / str(threshold)
                output = source / "visualization"
                report = lab6.run_analysis(lab6.DEFAULT_INPUT, source, threshold=threshold)
                lab7.visualize(source, output, topk=10)
                hubs = pd.read_csv(output / "hubs.csv")
                self.assertEqual(list(hubs.columns), columns)
                self.assertTrue((hubs.Degree > 0).all())
                self.assertTrue((hubs.Module >= 0).all())
                self.assertEqual(len(hubs), min(10, report["retained_genes"] - report["isolates"]))
                self.assertTrue((output / "network.png").is_file())
                if threshold == 1.0:
                    self.assertTrue(hubs.empty)


if __name__ == "__main__":
    unittest.main()
