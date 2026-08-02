import json
import tempfile
import unittest
from pathlib import Path

from campaign import grade_accepted, load_frontiers
from certificate_result import CertificateGrade


class CampaignGradeTests(unittest.TestCase):
    def test_legacy_success_is_not_rigorous(self):
        record = {"ok": True}
        self.assertFalse(grade_accepted(record, CertificateGrade.RIGOROUS))
        self.assertTrue(
            grade_accepted(record, CertificateGrade.EXPERIMENTAL)
        )

    def test_sampled_success_is_rejected_by_rigorous_policy(self):
        record = {"ok": True, "grade": "SAMPLED_BOUND"}
        self.assertFalse(grade_accepted(record, CertificateGrade.RIGOROUS))
        self.assertTrue(
            grade_accepted(record, CertificateGrade.SAMPLED_BOUND)
        )

    def test_resume_frontier_uses_only_accepted_grades(self):
        records = [
            {"ok": True, "grade": "SAMPLED_BOUND",
             "beta": [0.2, 1.0, 2.0], "hv": [0.1, 0.1, 0.1]},
            {"ok": True, "grade": "RIGOROUS",
             "beta": [0.1, 1.0, 2.0], "hv": [0.05, 0.1, 0.1]},
        ]
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "ledger.jsonl"
            path.write_text("".join(json.dumps(r) + "\n" for r in records))
            frontiers = load_frontiers(path, CertificateGrade.RIGOROUS)

        self.assertAlmostEqual(frontiers[(1.0, 2.0)], 0.15)


if __name__ == "__main__":
    unittest.main()
