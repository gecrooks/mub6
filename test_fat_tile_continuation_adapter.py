import unittest

from certificate_result import CertificateGrade, CertificateResult, Evidence
from continuation_artifact import RegularContinuationPayload
from coverage_contract import RootCoverageZone, SweepCoverageWitness
from fat_tile_continuation_adapter import (
    ContinuationAdapterError,
    adapt_fat_tile_result,
    continuation_readiness,
)


def source_result(n_roots=1):
    return CertificateResult.from_evidence(
        True,
        [Evidence("coverage", CertificateGrade.SAMPLED_BOUND,
                  "current fat-tile output")],
        metadata={"n_roots": n_roots},
    )


def components(grade=CertificateGrade.RIGOROUS):
    zone = RootCoverageZone(
        "tube", (0.0,) * 5, (0.2,) * 5, (0.1,) * 5,
        (0.01,) * 5, True, CertificateGrade.RIGOROUS,
    )
    parent = SweepCoverageWitness(
        (zone,), True, CertificateGrade.RIGOROUS, 10,
        parameter_center=(1.0, 2.0, 3.0),
        parameter_half_widths=(0.5, 0.5, 0.5),
    ).artifact()
    linear = tuple((0.0, 0.0, 0.0) for _ in range(5))
    quadratic = tuple(tuple((0.0, 0.0, 0.0) for _ in range(3))
                      for _ in range(5))
    payload = RegularContinuationPayload(
        0, zone.center, linear, quadratic, 1e-8, (1e-3,) * 5,
        (1e-2,) * 5, grade,
    )
    return parent, (payload,)


class FatTileContinuationAdapterTests(unittest.TestCase):
    def test_current_fat_tile_shape_reports_every_missing_proof_class(self):
        report = continuation_readiness(source_result())
        self.assertFalse(report.ready)
        self.assertEqual(report.source_grade, CertificateGrade.SAMPLED_BOUND)
        self.assertEqual(len(report.missing), 5)
        self.assertIn("ParentCoverageArtifact", report.missing[0])
        self.assertIn("certified S", report.missing[2])
        with self.assertRaises(ContinuationAdapterError) as caught:
            adapt_fat_tile_result(source_result())
        self.assertEqual(caught.exception.report, report)

    def test_complete_components_adapt_despite_separate_sampled_pair_result(self):
        parent, payloads = components()
        artifact = adapt_fat_tile_result(
            source_result(), parent=parent, continuations=payloads
        )
        self.assertEqual(artifact.grade, CertificateGrade.RIGOROUS)

    def test_sampled_continuation_is_rejected_by_rigorous_policy(self):
        parent, payloads = components(CertificateGrade.SAMPLED_BOUND)
        report = continuation_readiness(
            source_result(), parent=parent, continuations=payloads
        )
        self.assertFalse(report.ready)
        self.assertIn("below RIGOROUS", report.missing[0])

    def test_root_count_mismatch_is_rejected(self):
        parent, payloads = components()
        report = continuation_readiness(
            source_result(n_roots=2), parent=parent,
            continuations=payloads,
        )
        self.assertFalse(report.ready)
        self.assertIn("root count", report.missing[0])


if __name__ == "__main__":
    unittest.main()
