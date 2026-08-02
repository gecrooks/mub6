import unittest

from rigorous_census_report import summarize


class CensusReportTests(unittest.TestCase):
    def test_summary_tracks_largest_rigorous_rung_and_failures(self):
        records = [
            {"beta": [1, 2, 3], "hv": [0.1, 0.1, 0.1],
             "accepted_rigorous": True, "seconds": 2,
             "result": {"metadata": {"coverage_boxes": 100}}},
            {"beta": [1, 2, 3], "hv": [0.2, 0.2, 0.2],
             "accepted_rigorous": False, "seconds": 3,
             "result": {"reason": "stuck boxes", "metadata": {}}},
            {"beta": [4, 5, 6], "hv": [0.15, 0.15, 0.15],
             "accepted_rigorous": True, "seconds": 4,
             "result": {"metadata": {"coverage_boxes": 300}}},
        ]

        report = summarize(records)

        self.assertEqual(report["rigorous_passes"], 2)
        self.assertEqual(report["seconds_median"], 3)
        self.assertEqual(report["coverage_boxes_median"], 200)
        self.assertEqual(report["failures"], {"stuck boxes": 1})
        self.assertEqual(report["point_results"][0]["max_rigorous_h"], 0.1)

    def test_empty_summary_is_defined(self):
        report = summarize([])
        self.assertEqual(report["pass_fraction"], 0.0)
        self.assertIsNone(report["seconds_median"])


if __name__ == "__main__":
    unittest.main()
