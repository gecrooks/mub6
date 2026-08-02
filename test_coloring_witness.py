import unittest

import numpy as np

from coloring_witness import (
    ColoringStatus,
    critical_directions,
    find_coloring,
    verify_coloring,
)


def cycle(n):
    adj = np.zeros((n, n), dtype=bool)
    for i in range(n):
        adj[i, (i + 1) % n] = adj[(i + 1) % n, i] = True
    return adj


class ColoringWitnessTests(unittest.TestCase):
    def test_bipartite_graph_gets_verified_witness(self):
        result = find_coloring(cycle(20), k=2)
        self.assertEqual(result.status, ColoringStatus.COLORABLE)
        self.assertTrue(verify_coloring(cycle(20), result.colors, 2))

    def test_odd_cycle_is_not_two_colorable(self):
        result = find_coloring(cycle(5), k=2)
        self.assertEqual(result.status, ColoringStatus.UNCOLORABLE)
        self.assertIsNone(result.colors)

    def test_k6_is_not_five_colorable(self):
        adj = ~np.eye(6, dtype=bool)
        result = find_coloring(adj, k=5)
        self.assertEqual(result.status, ColoringStatus.UNCOLORABLE)

    def test_budget_exhaustion_is_unknown(self):
        adj = ~np.eye(8, dtype=bool)
        result = find_coloring(adj, k=4, node_budget=1)
        self.assertEqual(result.status, ColoringStatus.UNKNOWN)

    def test_critical_direction_scores_only_conflict_edges(self):
        adj = np.zeros((3, 3), dtype=bool)
        adj[0, 1] = adj[1, 0] = True
        uncertainty = np.zeros((3, 3, 3))
        uncertainty[0, 1] = uncertainty[1, 0] = (1.0, 4.0, 2.0)
        uncertainty[0, 2] = uncertainty[2, 0] = (100.0, 0.0, 0.0)

        ranked = critical_directions(adj, uncertainty)

        self.assertEqual(ranked, ((1, 4.0), (2, 2.0), (0, 1.0)))


if __name__ == "__main__":
    unittest.main()
