"""Deterministic coloring witnesses and refinement diagnostics.

The solver may prove that a conflict graph is k-colorable by exhibiting a
proper coloring. If its node budget expires it returns UNKNOWN, never a
false failure. Exhausting the finite search proves that this particular
graph is not k-colorable, which is a scheduling signal to refine the pair
bounds or parameter box—not a theorem failure.
"""

from dataclasses import dataclass
from enum import Enum

import numpy as np


class ColoringStatus(Enum):
    COLORABLE = "colorable"
    UNCOLORABLE = "uncolorable"
    UNKNOWN = "unknown"


@dataclass(frozen=True)
class ColoringWitness:
    status: ColoringStatus
    colors: tuple[int, ...] | None
    nodes: int
    critical_vertices: tuple[int, ...]

    @property
    def ok(self):
        return self.status is ColoringStatus.COLORABLE


def _adjacency(value):
    adj = np.asarray(value, dtype=bool)
    if adj.ndim != 2 or adj.shape[0] != adj.shape[1]:
        raise ValueError("adjacency must be square")
    if np.diag(adj).any():
        raise ValueError("adjacency diagonal must be false")
    if not np.array_equal(adj, adj.T):
        raise ValueError("adjacency must be symmetric")
    return adj


def verify_coloring(adjacency, colors, k=None):
    adj = _adjacency(adjacency)
    colors = np.asarray(colors, dtype=int)
    if colors.shape != (len(adj),) or (colors < 0).any():
        return False
    if k is not None and (colors >= k).any():
        return False
    rows, cols = np.where(np.triu(adj, 1))
    return not bool(np.any(colors[rows] == colors[cols]))


def find_coloring(adjacency, k=5, node_budget=1_000_000):
    """Exact DSATUR search with a fail-closed node budget."""
    adj = _adjacency(adjacency)
    n = len(adj)
    if k <= 0 or node_budget <= 0:
        raise ValueError("k and node_budget must be positive")
    if n == 0:
        return ColoringWitness(ColoringStatus.COLORABLE, (), 0, ())

    colors = -np.ones(n, dtype=int)
    degree = adj.sum(axis=1).astype(int)
    branch_count = np.zeros(n, dtype=int)
    nodes = 0
    budget_hit = False

    def choose_vertex():
        uncolored = np.where(colors < 0)[0]
        best = None
        best_key = None
        for vertex in uncolored:
            used = set(colors[adj[vertex]]) - {-1}
            key = (len(used), degree[vertex], -int(vertex))
            if best_key is None or key > best_key:
                best_key = key
                best = int(vertex)
        return best

    def search(colored_count, max_used):
        nonlocal nodes, budget_hit
        if colored_count == n:
            return True
        if nodes >= node_budget:
            budget_hit = True
            return False
        nodes += 1
        vertex = choose_vertex()
        branch_count[vertex] += 1
        forbidden = set(colors[adj[vertex]]) - {-1}
        existing = [c for c in range(max_used + 1) if c not in forbidden]
        # Least-constraining existing colors first.
        existing.sort(key=lambda c: int(np.sum(
            (colors < 0) & adj[vertex]
            & np.array([not np.any(colors[adj[u]] == c)
                        for u in range(n)], dtype=bool)
        )), reverse=True)
        choices = existing
        if max_used + 1 < k:
            choices.append(max_used + 1)  # one canonical new color
        for color in choices:
            colors[vertex] = color
            if search(colored_count + 1, max(max_used, color)):
                return True
            colors[vertex] = -1
            if budget_hit:
                return False
        return False

    found = search(0, -1)
    order = np.lexsort((np.arange(n), -degree, -branch_count))
    critical = tuple(int(v) for v in order[:min(10, n)])
    if found:
        witness = tuple(int(c) for c in colors)
        if not verify_coloring(adj, witness, k):
            raise AssertionError("internal coloring witness is invalid")
        return ColoringWitness(ColoringStatus.COLORABLE, witness, nodes,
                               critical)
    status = ColoringStatus.UNKNOWN if budget_hit else ColoringStatus.UNCOLORABLE
    return ColoringWitness(status, None, nodes, critical)


def critical_directions(adjacency, pair_uncertainty):
    """Rank parameter axes by uncertainty on current conflict edges.

    ``pair_uncertainty[a,b,j]`` is a nonnegative bound charged to pair
    ``(a,b)`` along parameter direction ``j``. The result is advisory
    scheduling data; subdivision still needs independently certified child
    bounds.
    """
    adj = _adjacency(adjacency)
    uncertainty = np.asarray(pair_uncertainty, dtype=float)
    if uncertainty.ndim != 3 or uncertainty.shape[:2] != adj.shape:
        raise ValueError("pair_uncertainty must have shape (n,n,d)")
    if not np.isfinite(uncertainty).all() or (uncertainty < 0).any():
        raise ValueError("pair uncertainty must be finite and nonnegative")
    edge_mask = np.triu(adj, 1)
    scores = uncertainty[edge_mask].sum(axis=0)
    order = np.argsort(-scores, kind="stable")
    return tuple((int(axis), float(scores[axis])) for axis in order)
