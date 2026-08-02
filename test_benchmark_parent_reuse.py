import unittest

from benchmark_parent_reuse import representative_indices, summarize


class ParentReuseBenchmarkTests(unittest.TestCase):
    def test_representative_indices_include_endpoints(self):
        self.assertEqual(representative_indices(8, 3), (0, 4, 7))
        self.assertEqual(representative_indices(2, 3), (0, 1))

    def test_summary_includes_parent_cost_and_break_even(self):
        result = summarize(20.0, 8, (1.0, 1.0), (6.0, 6.0))
        self.assertEqual(result["break_even_children"], 4)
        self.assertEqual(result["projected_all_children_reuse_seconds"], 28)
        self.assertEqual(
            result["projected_all_children_independent_seconds"], 48
        )
        self.assertAlmostEqual(result["projected_all_children_speedup"],
                               48 / 28)


if __name__ == "__main__":
    unittest.main()
