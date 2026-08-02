import json
import unittest

from certificate_result import (
    CertificateGrade,
    CertificateResult,
    Evidence,
    combine_results,
)


class CertificateResultTests(unittest.TestCase):
    def test_sampled_dependency_downgrades_result(self):
        result = CertificateResult.from_evidence(
            True,
            [
                Evidence("interval sweep", CertificateGrade.RIGOROUS),
                Evidence("cache curvature", CertificateGrade.SAMPLED_BOUND),
            ],
        )

        self.assertEqual(result.grade, CertificateGrade.SAMPLED_BOUND)
        self.assertFalse(result.accepted_at(CertificateGrade.RIGOROUS))
        self.assertTrue(result.accepted_at(CertificateGrade.SAMPLED_BOUND))

    def test_constructor_rejects_grade_above_weakest_dependency(self):
        with self.assertRaisesRegex(ValueError, "exceeds weakest dependency"):
            CertificateResult(
                True,
                CertificateGrade.RIGOROUS,
                (Evidence("sample", CertificateGrade.EXPERIMENTAL),),
            )

    def test_failure_is_never_accepted(self):
        result = CertificateResult.from_evidence(
            False,
            [Evidence("interval sweep", CertificateGrade.RIGOROUS)],
            reason="unresolved boxes",
        )

        self.assertFalse(result.accepted_at(CertificateGrade.EXPERIMENTAL))
        self.assertFalse(result)

    def test_combination_inherits_weakest_stage_and_failure(self):
        rigorous = CertificateResult.from_evidence(
            True,
            [Evidence("coverage", CertificateGrade.RIGOROUS)],
        )
        sampled_failure = CertificateResult.from_evidence(
            False,
            [Evidence("curvature", CertificateGrade.SAMPLED_BOUND)],
            reason="cache remainder unavailable",
        )

        combined = combine_results("tile", [rigorous, sampled_failure])

        self.assertFalse(combined.ok)
        self.assertEqual(combined.grade, CertificateGrade.SAMPLED_BOUND)
        self.assertEqual(combined.reason, "cache remainder unavailable")
        self.assertEqual(combined.metadata["stage_count"], 2)

    def test_ceiling_can_mark_an_otherwise_rigorous_experiment(self):
        result = CertificateResult.from_evidence(
            True,
            [Evidence("sweep", CertificateGrade.RIGOROUS)],
            ceiling=CertificateGrade.EXPERIMENTAL,
        )

        self.assertEqual(result.grade, CertificateGrade.EXPERIMENTAL)

    def test_dict_form_is_json_serializable(self):
        result = CertificateResult.from_evidence(
            True,
            [Evidence("coverage", CertificateGrade.RIGOROUS, "all boxes")],
            metadata={"boxes": 12},
        )

        restored = json.loads(json.dumps(result.as_dict()))
        self.assertEqual(restored["grade"], "RIGOROUS")
        self.assertEqual(restored["metadata"], {"boxes": 12})
        self.assertTrue(result)

    def test_empty_evidence_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "declare its evidence"):
            CertificateResult.from_evidence(True, [])


if __name__ == "__main__":
    unittest.main()
