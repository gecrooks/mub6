import json
import tempfile
import unittest
from pathlib import Path

from campaign import grade_accepted, load_frontiers
from campaign_coverage import analyze_ledger_records
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
             "beta": [0.05, 1.0, 2.0], "hv": [0.05, 0.1, 0.1]},
        ]
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "ledger.jsonl"
            path.write_text("".join(json.dumps(r) + "\n" for r in records))
            frontiers = load_frontiers(path, CertificateGrade.RIGOROUS)

        self.assertAlmostEqual(frontiers[(1.0, 2.0)], 0.1)

    def test_disconnected_success_cannot_jump_frontier(self):
        records = [
            {"ok": True, "grade": "RIGOROUS",
             "beta": [0.05, 1.0, 2.0], "hv": [0.05, 0.1, 0.1]},
            {"ok": True, "grade": "RIGOROUS",
             "beta": [0.35, 1.0, 2.0], "hv": [0.05, 0.1, 0.1]},
        ]
        report = analyze_ledger_records(records, 0.0, 0.5)[(1.0, 2.0)]

        self.assertAlmostEqual(report.frontier, 0.1)
        self.assertEqual(len(report.gaps), 2)
        self.assertEqual(len(report.islands), 1)
        for actual, expected in zip(report.gaps,
                                    ((0.1, 0.3), (0.4, 0.5))):
            self.assertAlmostEqual(actual[0], expected[0])
            self.assertAlmostEqual(actual[1], expected[1])
        self.assertAlmostEqual(report.islands[0][0], 0.3)
        self.assertAlmostEqual(report.islands[0][1], 0.4)

    def test_out_of_order_overlaps_form_one_component(self):
        records = [
            {"ok": True, "grade": "RIGOROUS",
             "beta": [0.25, 1.0, 2.0], "hv": [0.1, 0.1, 0.1]},
            {"ok": True, "grade": "RIGOROUS",
             "beta": [0.075, 1.0, 2.0], "hv": [0.075, 0.1, 0.1]},
        ]
        report = analyze_ledger_records(records, 0.0, 0.5)[(1.0, 2.0)]

        self.assertAlmostEqual(report.frontier, 0.35)
        self.assertEqual(report.islands, ())

    def test_duplicates_are_counted_but_do_not_change_union(self):
        record = {"ok": True, "grade": "RIGOROUS",
                  "beta": [0.05, 1.0, 2.0], "hv": [0.05, 0.1, 0.1]}
        report = analyze_ledger_records([record, record], 0.0, 0.2)[(1.0, 2.0)]

        self.assertEqual(report.accepted_records, 1)
        self.assertEqual(report.duplicate_records, 1)

    def test_boundary_overhang_connects_to_domain_start(self):
        record = {"ok": True, "grade": "RIGOROUS",
                  "beta": [0.01, 1.0, 2.0], "hv": [0.02, 0.1, 0.1]}
        report = analyze_ledger_records([record], 0.0, 0.1)[(1.0, 2.0)]

        self.assertAlmostEqual(report.frontier, 0.03)
        self.assertEqual(len(report.gaps), 1)
        self.assertAlmostEqual(report.gaps[0][0], 0.03)

    def test_bad_grade_is_rejected_fail_closed(self):
        record = {"ok": True, "grade": "CERTIFIED",
                  "beta": [0.05, 1.0, 2.0], "hv": [0.05, 0.1, 0.1]}
        report = analyze_ledger_records([record], 0.0, 0.1)[(1.0, 2.0)]

        self.assertEqual(report.frontier, 0.0)
        self.assertEqual(report.rejected_records, 1)


if __name__ == "__main__":
    unittest.main()
