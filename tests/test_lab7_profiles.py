"""Lecture-alignment checks through the optional Lab 7 command-line driver."""
import importlib.util
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DRIVER = ROOT / "labs/07_network_viz/ex02_module_profiles.py"


class ModuleProfileTests(unittest.TestCase):
    def test_saved_modules_produce_the_lecture_profiles_despite_metadata_row_order(self):
        spec = importlib.util.spec_from_file_location(
            "lab6_profiles_test", ROOT / "labs/06_wgcna/ex01_gce_networks.py",
        )
        lab6 = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(lab6)
        with tempfile.TemporaryDirectory() as temporary:
            folder = Path(temporary)
            source, output = folder / "lab06", folder / "profiles"
            lab6.run_analysis(lab6.DEFAULT_INPUT, source, metadata_csv=lab6.DEFAULT_METADATA)
            metadata = pd.read_csv(source / "sample_metadata.csv")
            metadata.iloc[::-1].to_csv(source / "sample_metadata.csv", index=False)
            result = subprocess.run(
                [sys.executable, "-X", "utf8", str(DRIVER), "--handle", "test", "--lab06-dir", str(source),
                 "--output-dir", str(output)],
                capture_output=True, text=True, env={**os.environ, "MPLBACKEND": "Agg"},
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            summary = pd.read_csv(output / "module_summary.csv").set_index("Module")
            self.assertEqual(summary.Genes.tolist(), [30, 29, 30])
            for module, expected in enumerate([0.84747284, 0.20792757, -0.84490808]):
                self.assertAlmostEqual(summary.loc[module, "TreatmentCorrelation"], expected, places=6)
            self.assertAlmostEqual(summary.loc[0, "PC1VarianceExplained"], 0.78302414, places=6)
            profiles = pd.read_csv(output / "module_profiles.csv", index_col=0)
            self.assertEqual(profiles.shape, (48, 3))
            self.assertEqual(profiles.index[0], "S01")
            self.assertTrue((output / "module_profiles.png").is_file())

    def test_incomplete_duplicate_or_unknown_sample_metadata_is_rejected_before_export(self):
        spec = importlib.util.spec_from_file_location(
            "lab6_bad_metadata", ROOT / "labs/06_wgcna/ex01_gce_networks.py",
        )
        lab6 = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(lab6)
        with tempfile.TemporaryDirectory() as temporary:
            folder = Path(temporary)
            source = folder / "lab06"
            lab6.run_analysis(lab6.DEFAULT_INPUT, source, metadata_csv=lab6.DEFAULT_METADATA)
            metadata = pd.read_csv(source / "sample_metadata.csv")
            duplicate = metadata.copy()
            duplicate.loc[1, "Sample"] = duplicate.loc[0, "Sample"]
            unknown = metadata.copy()
            unknown.loc[0, "Condition"] = "unknown"
            for label, invalid in [("missing", metadata.iloc[1:]), ("duplicate", duplicate),
                                   ("unknown", unknown)]:
                with self.subTest(case=label):
                    invalid.to_csv(source / "sample_metadata.csv", index=False)
                    output = folder / label
                    result = subprocess.run(
                        [sys.executable, "-X", "utf8", str(DRIVER), "--handle", "test",
                         "--lab06-dir", str(source), "--output-dir", str(output)],
                        capture_output=True, text=True, encoding="utf-8",
                        env={**os.environ, "MPLBACKEND": "Agg"},
                    )
                    self.assertNotEqual(result.returncode, 0)
                    self.assertIn("Folosiți exportul complet", result.stderr)
                    self.assertFalse(output.exists())


if __name__ == "__main__":
    unittest.main()
