"""Small invariants for preserving cohort identity and recorded boundaries."""

import json
import tempfile
import unittest
from pathlib import Path

from export_exact_graphbpe import export_records
from summarize_token_lengths import read_records, same_molecule, statistics


class HelperTests(unittest.TestCase):
    def test_duplicate_records_rejected(self):
        with tempfile.TemporaryDirectory() as root:
            source = Path(root) / "records.jsonl"
            source.write_text('{"record_index":1}\n{"record_index":1}\n', encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "Duplicate"):
                read_records(source)

    def test_matching_id_does_not_hide_smiles_change(self):
        left = dict(record_index=1, chembl_id="X", canonical_smiles="CCO")
        right = dict(record_index=1, chembl_id="X", smiles="CCN")
        with self.assertRaisesRegex(ValueError, "Molecule mismatch"):
            same_molecule(left, right)

    def test_linear_percentile(self):
        self.assertAlmostEqual(statistics([1, 2, 3, 4])["p95"], 3.85)

    def test_partition_check_and_original_order(self):
        row = dict(
            bpe_boundary_exact=True,
            bpe_node_ids=[9, 2],
            bpe_groups=[[0, 1], [2]],
            safe_groups=[[2], [0, 1]],
            safe="C1.CC1",
            num_bpe_fragments=2,
            num_safe_fragments=2,
            roundtrip_ok=True,
        )
        with tempfile.TemporaryDirectory() as root:
            source = Path(root) / "records.jsonl"
            source.write_text(
                json.dumps(row) + "\n" + json.dumps(dict(bpe_boundary_exact=False)) + "\n",
                encoding="utf-8",
            )
            summary = export_records(source, Path(root) / "export")
            self.assertEqual(summary["kept"], 1)
            self.assertEqual(summary["excluded_boundary_changed"], 1)
            saved = json.loads((Path(root) / "export/records.jsonl").read_text())
            self.assertEqual(saved["node_ids"], [9, 2])
            self.assertFalse(saved["ids_are_lossless"])
            row["safe_groups"] = [[0, 2], [1]]
            source.write_text(json.dumps(row) + "\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "partitions differ"):
                export_records(source, Path(root) / "bad")


if __name__ == "__main__":
    unittest.main()
