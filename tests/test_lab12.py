"""Regression checks for teaching materials, not solutions to student TODOs."""
import importlib.util
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from Bio import SeqIO

ROOT = Path(__file__).resolve().parents[1]


def module(path):
    spec = importlib.util.spec_from_file_location("teaching_demo", ROOT / path)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


class Lab12Tests(unittest.TestCase):
    def test_correct_protein_and_combined_dataset(self):
        fish = SeqIO.read(ROOT / "data/sample/tp53_dr_protein.fasta", "fasta")
        self.assertEqual(fish.id, "sp|P79734|P53_DANRE")
        self.assertEqual(len(fish), 373)
        self.assertIn("OS=Danio rerio", fish.description)
        with (ROOT / "data/sample/tp53_protein_multi.fasta").open() as handle:
            records = list(SeqIO.parse(handle, "fasta"))
        self.assertEqual([r.id for r in records],
                         ["sp|P04637|P53_HUMAN", "sp|P02340|P53_MOUSE", fish.id])
        self.assertEqual(records[-1].seq, fish.seq)

    def test_offline_brca1_is_defined_transcript(self):
        demo = module("labs/01_intro&databases/demo01_entrez_brca1.py")
        record = demo.validate_record(SeqIO.read(demo.CACHE, "genbank"))
        self.assertEqual(len(record), 7088)
        self.assertTrue(any(f.type == "CDS" for f in record.features))
        run = subprocess.run([sys.executable, str(ROOT / "labs/01_intro&databases/demo01_entrez_brca1.py")],
                             capture_output=True, text=True, timeout=20)
        self.assertEqual(run.returncode, 0, run.stderr)
        self.assertIn("NM_007294.4", run.stdout)

    def test_distances_require_corresponding_sites(self):
        demo = module("labs/02_alignment/demo02_distance_matrix.py")
        self.assertEqual(demo.aligned_counts("AC-GTN", "ATCGTA"), (1, 4))
        self.assertEqual(demo.aligned_counts("NN--", "--NN"), (0, 0))
        with self.assertRaises(ValueError):
            demo.aligned_counts("ACG", "AC")
        # An inserted base must not shift every following comparison.
        self.assertEqual(demo.pairwise_counts("ACGT", "ATCGT"), (0, 4))

    def test_no_positive_local_match_and_zero_sites(self):
        with tempfile.TemporaryDirectory() as tmp:
            fasta = Path(tmp) / "sequences.fasta"
            fasta.write_text(">a\nAAAA\n>b\nCCCC\n")
            run = subprocess.run([sys.executable, str(ROOT / "labs/02_alignment/demo01_pairwise.py"),
                                  "--fasta", str(fasta)], capture_output=True, text=True, timeout=20)
            self.assertEqual(run.returncode, 0, run.stderr)
            self.assertIn("No positive-scoring local alignment", run.stdout)
            fasta.write_text(">a\nNN--\n>b\n--NN\n")
            run = subprocess.run([sys.executable, str(ROOT / "labs/02_alignment/demo02_distance_matrix.py"),
                                  "--fasta", str(fasta), "--aligned"],
                                 capture_output=True, text=True, timeout=20)
            self.assertEqual(run.returncode, 0, run.stderr)
            self.assertIn("a,b,0,NA,0", run.stdout)

    def test_update_manifest_does_not_touch_student_work(self):
        path = ROOT / "docs/updates/lab12-2026-10-05.paths"
        entries = path.read_text().splitlines()
        self.assertEqual(len(entries), len(set(entries)))
        for entry in entries:
            self.assertTrue((ROOT / entry).is_file(), entry)
            self.assertFalse(Path(entry).is_absolute(), entry)
            self.assertFalse(set(Path(entry).parts) & {"..", ".git", "submissions", "roster"}, entry)
            self.assertFalse(entry.startswith("data/work/"), entry)


if __name__ == "__main__":
    unittest.main()
