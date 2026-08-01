import unittest

import numpy as np

from coverage_contract import BoxDisposition, close_survivor_boxes


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


if __name__ == "__main__":
    unittest.main()
