import unittest

from campaign_geometry import (
    axis_tile_count,
    box_tile_count,
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


if __name__ == "__main__":
    unittest.main()
