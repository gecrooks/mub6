import unittest
import json

import numpy as np

from certificate_result import CertificateGrade
from coverage_contract import (
    BoxDisposition,
    RootCoverageZone,
    SweepCoverageWitness,
    close_survivor_boxes,
)


class CoverageContractTests(unittest.TestCase):
    def test_mixed_blob_fails_if_one_member_is_unresolved(self):
        report = close_survivor_boxes(
            centers=[[0.01, 0.01], [0.30, 0.30]],
            half_widths=[[0.01, 0.01], [0.01, 0.01]],
            root_centers=[[0.0, 0.0]],
            zone_radii=[[0.05, 0.05]],
        )

        self.assertEqual(
            report.dispositions,
            (BoxDisposition.COLLECTED, BoxDisposition.UNRESOLVED),
        )
        self.assertFalse(report.complete)
        self.assertEqual(report.unresolved_indices, (1,))

    def test_candidate_cannot_be_satisfied_by_a_different_root(self):
        report = close_survivor_boxes(
            centers=[[1.0, 1.0]],
            half_widths=[[0.01, 0.01]],
            root_centers=[[0.0, 0.0], [1.0, 1.0]],
            zone_radii=[[0.1, 0.1], [0.1, 0.1]],
            candidate_roots=[0],
        )

        self.assertEqual(
            report.dispositions, (BoxDisposition.UNRESOLVED,)
        )
        self.assertEqual(report.root_indices, (None,))

    def test_componentwise_width_must_fit(self):
        report = close_survivor_boxes(
            centers=[[0.0, 0.0]],
            half_widths=[[0.04, 0.06]],
            root_centers=[[0.0, 0.0]],
            zone_radii=[[0.05, 0.05]],
        )

        self.assertFalse(report.complete)

    def test_independent_exclusion_closes_an_uncollected_box(self):
        report = close_survivor_boxes(
            centers=[[0.5, 0.5]],
            half_widths=[[0.01, 0.01]],
            root_centers=[[0.0, 0.0]],
            zone_radii=[[0.05, 0.05]],
            excluded=[True],
        )

        self.assertEqual(report.dispositions, (BoxDisposition.EXCLUDED,))
        self.assertTrue(report.complete)

    def test_torus_wraparound_is_respected(self):
        report = close_survivor_boxes(
            centers=[[np.pi - 0.01]],
            half_widths=[[0.01]],
            root_centers=[[-np.pi + 0.01]],
            zone_radii=[[0.04]],
        )

        self.assertEqual(report.dispositions, (BoxDisposition.COLLECTED,))

    def test_empty_roots_leave_survivors_unresolved(self):
        report = close_survivor_boxes(
            centers=[[0.0, 0.0]],
            half_widths=[[0.01, 0.01]],
            root_centers=np.empty((0, 2)),
            zone_radii=np.empty((0, 2)),
        )

        self.assertFalse(report.complete)
        self.assertEqual(report.unresolved_indices, (0,))

    def test_no_survivors_is_vacuously_complete(self):
        report = close_survivor_boxes(
            centers=np.empty((0, 2)),
            half_widths=np.empty((0, 2)),
            root_centers=np.empty((0, 2)),
            zone_radii=np.empty((0, 2)),
        )

        self.assertTrue(report.complete)

    def test_sweep_witness_exports_guard_handoff_and_tube(self):
        zone = RootCoverageZone(
            kind="tube",
            center=(0.0, 0.0),
            guard_radii=(0.2, 0.3),
            collected_reach=(0.1, 0.25),
            enclosure_radii=(0.01, 0.02),
            handoff_complete=True,
            grade=CertificateGrade.RIGOROUS,
        )
        witness = SweepCoverageWitness(
            zones=(zone,),
            global_sweep_complete=True,
            arithmetic_grade=CertificateGrade.RIGOROUS,
            boxes_processed=123,
            parameter_center=(1.0,),
            parameter_half_widths=(0.1,),
        )

        self.assertTrue(witness.complete)
        self.assertEqual(witness.evidence().grade, CertificateGrade.RIGOROUS)
        self.assertEqual(
            witness.as_dict()["zones"][0]["enclosure_radii"], [0.01, 0.02]
        )
        restored = SweepCoverageWitness.from_dict(
            json.loads(json.dumps(witness.as_dict()))
        )
        self.assertEqual(restored, witness)
        self.assertTrue(restored.matches((1.0,), (0.05,), [[0.0, 0.0]]))
        self.assertFalse(restored.matches((1.1,), (0.05,), [[0.0, 0.0]]))

    def test_open_handoff_keeps_sweep_witness_incomplete(self):
        zone = RootCoverageZone(
            kind="fold",
            center=(0.0,),
            guard_radii=(0.2,),
            collected_reach=(0.1,),
            enclosure_radii=None,
            handoff_complete=False,
            grade=CertificateGrade.SAMPLED_BOUND,
        )
        witness = SweepCoverageWitness(
            zones=(zone,),
            global_sweep_complete=True,
            arithmetic_grade=CertificateGrade.SAMPLED_BOUND,
            boxes_processed=5,
            parameter_center=(1.0,),
            parameter_half_widths=(0.1,),
        )

        self.assertFalse(witness.complete)

    def test_tube_zone_requires_final_enclosure(self):
        with self.assertRaisesRegex(ValueError, "require final enclosure"):
            RootCoverageZone(
                kind="tube",
                center=(0.0,),
                guard_radii=(0.2,),
                collected_reach=(0.1,),
                enclosure_radii=None,
                handoff_complete=True,
                grade=CertificateGrade.RIGOROUS,
            )


if __name__ == "__main__":
    unittest.main()
