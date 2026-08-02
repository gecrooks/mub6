import json
import unittest

import numpy as np

from certificate_result import CertificateGrade
from continuation_artifact import (
    ContinuationCoverageArtifact,
    FoldContinuationPayload,
    FoldOraclePayload,
    RegularContinuationPayload,
)
from coverage_contract import RootCoverageZone, SweepCoverageWitness


def regular_payload(index=0, seed=(0.0,) * 5):
    linear = tuple(tuple(1.0 if i == j else 0.0 for j in range(3))
                   for i in range(5))
    quadratic = tuple(tuple(tuple(0.0 for _ in range(3))
                            for _ in range(3)) for _ in range(5))
    return RegularContinuationPayload(
        index, seed, linear, quadratic, 1e-8, (1e-3,) * 5,
        (1e-2,) * 5, CertificateGrade.RIGOROUS,
    )


def parent(zones):
    return SweepCoverageWitness(
        zones=tuple(zones), global_sweep_complete=True,
        arithmetic_grade=CertificateGrade.RIGOROUS, boxes_processed=10,
        parameter_center=(1.0, 2.0, 3.0),
        parameter_half_widths=(0.5, 0.5, 0.5),
    ).artifact()


class ContinuationArtifactTests(unittest.TestCase):
    def test_regular_round_trip_and_predicted_child_seed(self):
        zone = RootCoverageZone(
            "tube", (0.0,) * 5, (0.2,) * 5, (0.1,) * 5,
            (0.01,) * 5, True, CertificateGrade.RIGOROUS,
        )
        artifact = ContinuationCoverageArtifact(
            parent((zone,)), (regular_payload(),)
        )
        restored = ContinuationCoverageArtifact.from_dict(
            json.loads(json.dumps(artifact.as_dict()))
        )
        self.assertEqual(restored, artifact)
        child = restored.restrict((1.1, 2.2, 3.3), (0.1, 0.1, 0.1))
        self.assertTrue(child.ready_for_child_verification)
        self.assertTrue(child.roots[0].requires_krawczyk)
        np.testing.assert_allclose(child.roots[0].proposed_seed[:3],
                                   (0.1, 0.2, 0.3), rtol=0.0, atol=1e-15)
        outside = restored.restrict((1.6, 2.0, 3.0), (0.1, 0.1, 0.1))
        self.assertFalse(outside.ready_for_child_verification)

    def test_sampled_payload_downgrades_v2_artifact(self):
        zone = RootCoverageZone(
            "tube", (0.0,) * 5, (0.2,) * 5, (0.1,) * 5,
            (0.01,) * 5, True, CertificateGrade.RIGOROUS,
        )
        payload = regular_payload()
        sampled = RegularContinuationPayload(
            payload.root_index, payload.seed, payload.linear,
            payload.quadratic, payload.residual_bound,
            payload.jacobian_remainder, payload.tube_radii,
            CertificateGrade.SAMPLED_BOUND,
        )
        artifact = ContinuationCoverageArtifact(
            parent((zone,)), (sampled,)
        )
        self.assertEqual(artifact.grade, CertificateGrade.SAMPLED_BOUND)

    def test_missing_or_mismatched_payload_fails_closed(self):
        zone = RootCoverageZone(
            "tube", (0.0,) * 5, (0.2,) * 5, (0.1,) * 5,
            (0.01,) * 5, True, CertificateGrade.RIGOROUS,
        )
        with self.assertRaisesRegex(ValueError, "one continuation"):
            ContinuationCoverageArtifact(parent((zone,)), ())
        with self.assertRaisesRegex(ValueError, "seed must equal"):
            ContinuationCoverageArtifact(
                parent((zone,)), (regular_payload(seed=(0.1,) * 5),)
            )

    def test_digest_tampering_is_rejected(self):
        zone = RootCoverageZone(
            "tube", (0.0,) * 5, (0.2,) * 5, (0.1,) * 5,
            (0.01,) * 5, True, CertificateGrade.RIGOROUS,
        )
        value = ContinuationCoverageArtifact(
            parent((zone,)), (regular_payload(),)
        ).as_dict()
        value["continuations"][0]["residual_bound"] = 2e-8
        with self.assertRaisesRegex(ValueError, "digest mismatch"):
            ContinuationCoverageArtifact.from_dict(value)

    def test_split_fold_requires_exactly_eight_octants(self):
        oracle = FoldOraclePayload(
            (-1, -1, -1), (0.0,) * 5, (1.0, 0.0, 0.0, 0.0, 0.0),
            tuple(tuple(1.0 if i == j else 0.0 for j in range(4))
                  for i in range(5)),
            (-1.0, 1.0), ((0.0,) * 4, (0.0,) * 4), (0.1,) * 4,
            1.0, CertificateGrade.RIGOROUS,
        )
        with self.assertRaisesRegex(ValueError, "all 8 octants"):
            FoldContinuationPayload(0, "fold-split", (oracle,))


if __name__ == "__main__":
    unittest.main()
