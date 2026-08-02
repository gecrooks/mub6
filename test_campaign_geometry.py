import unittest

from campaign_geometry import (
    axis_tile_count,
    box_tile_count,
    parent_reuse_comparison,
    wall_band_comparison,
)


class CampaignGeometryTests(unittest.TestCase):
    def test_closed_cover_rounds_up(self):
        self.assertEqual(axis_tile_count(1.0, 0.3), 2)
        self.assertEqual(box_tile_count((1.0, 2.0), (0.3, 0.6)), 4)

    def test_anisotropic_wall_saves_two_tangential_ratios(self):
        result = wall_band_comparison(
            tangent_lengths=(1.0, 1.0),
            band_width=2e-6,
            normal_half_width=1e-6,
            tangent_half_width=25e-6,
        )

        self.assertEqual(result["anisotropic_tiles"], 400_000_000)
        self.assertEqual(result["isotropic_tiles"], 250_000_000_000)
        self.assertEqual(result["tile_reduction"], 625.0)

    def test_symmetry_reduces_count_fail_closed_by_ceiling(self):
        self.assertEqual(box_tile_count((1.0,), (0.1,), symmetry_factor=3), 2)

    def test_parent_reuse_amortizes_only_coverage_leg(self):
        result = parent_reuse_comparison(
            parent_tiles=10, children_per_parent=100,
            coverage_seconds=5.0, pair_seconds=1.0,
        )

        self.assertAlmostEqual(result["repeated_a100_hours"], 6000 / 3600)
        self.assertAlmostEqual(result["inherited_a100_hours"], 1050 / 3600)
        self.assertAlmostEqual(result["hour_reduction"], 6000 / 1050)


if __name__ == "__main__":
    unittest.main()
