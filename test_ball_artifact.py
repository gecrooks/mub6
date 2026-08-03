"""Tests for ball_coverage_artifact + fat_tile emission wiring."""

import copy
import math
import unittest

import numpy as np

from ball_coverage_artifact import (BallCoverageArtifact, DriftBoundClaim,
                                    SweepReplaySpec, ball_readiness,
                                    float_hex, hex_float)
from certificate_result import CertificateGrade, CertificateResult, Evidence

# tj - ti stepping by pi/3 makes the six-term phasor sum vanish
# exactly: a possibly-orthogonal pair at any budget.
ORTH_STEP = tuple((k + 1) * math.pi / 3.0 for k in range(5))


def _drift(beta=(0.5, 1.0, 2.0), h=0.01, bu=1.2):
    return DriftBoundClaim(
        derivation="rates.certified_rates(beta, (h,h,h), mv=True)"
                   "['beta_unit_vec'] sup",
        beta=beta, half_widths=(h, h, h), bu_max=bu, drift=bu * h * 3.0)


def _sweep(uncovered=(), frontier=False):
    return SweepReplaySpec(
        sweep="starve.fat_sweep_hulls", wmin=0.025, cell=0.025,
        tax_derivation="rates.certified_rates(mv=True) tax constants",
        boxes_swept=1000, hull_cells=3, uncovered=uncovered,
        frontier_complete=frontier)


def _artifact(**kw):
    balls = kw.pop("balls", ((0.0,) * 5, (0.01, -0.02, 0.03, 0.0, 0.01),
                             ORTH_STEP))
    drift = kw.pop("drift", _drift())
    r_loc = kw.pop("r_loc", 0.06)
    defaults = dict(
        parameter_center=(0.5, 1.0, 2.0),
        parameter_half_widths=(0.01, 0.01, 0.01),
        r_loc=r_loc, balls=balls, drift=drift, sweep=_sweep(),
        budget_w=r_loc + drift.drift,
        # balls 0,1 are near-identical phases: certified
        # non-orthogonal; ball 2 is the exact-orthogonality pattern,
        # adjacent to both.
        nonorth_pairs=((0, 1, 0.9),),
        coloring=(0, 0, 1), chi_bound=2)
    defaults.update(kw)
    return BallCoverageArtifact(**defaults)


class TestFloatBits(unittest.TestCase):
    def test_round_trip_exact(self):
        for x in (0.1, 1.0 / 3.0, 1e-300, -math.pi, 5e-324, 0.0, -0.0):
            self.assertEqual(hex_float(float_hex(x)), x)
            self.assertEqual(
                np.float64(hex_float(float_hex(x))).tobytes(),
                np.float64(x).tobytes())

    def test_rejects_nonfinite_and_malformed(self):
        for bad in (float("nan"), float("inf")):
            with self.assertRaises(ValueError):
                float_hex(bad)
        with self.assertRaises(ValueError):
            hex_float("7ff0000000000000")  # +inf bit pattern
        with self.assertRaises(ValueError):
            hex_float("0.1")


class TestArtifact(unittest.TestCase):
    def test_round_trip_and_digest(self):
        a = _artifact()
        d = a.as_dict()
        b = BallCoverageArtifact.from_dict(d)
        self.assertEqual(a, b)
        self.assertEqual(a.artifact_id, b.artifact_id)
        d2 = copy.deepcopy(d)
        d2["artifact_id"] = a.artifact_id
        BallCoverageArtifact.from_dict(d2)  # digest verifies
        d2["balls"][0][0] = float_hex(1e-9)
        with self.assertRaises(ValueError):
            BallCoverageArtifact.from_dict(d2)

    def test_kernel_check_passes(self):
        self.assertEqual(_artifact().check(), [])

    def test_kernel_catches_improper_coloring(self):
        a = _artifact(coloring=(0, 0, 0), chi_bound=1)
        self.assertTrue(any("not proper" in f for f in a.check()))

    def test_kernel_catches_inflated_pair_claim(self):
        # claim balls 0 and 2 (exact orthogonality pattern) are
        # certified non-orthogonal: the replay must refute it
        a = _artifact(nonorth_pairs=((0, 1, 0.9), (0, 2, 0.5)),
                      coloring=(0, 0, 1))
        self.assertTrue(any("refutes" in f for f in a.check()))

    def test_fail_closed_validation(self):
        with self.assertRaises(ValueError):
            _artifact(r_loc=0.9)  # beyond ball-vertex lemma radius
        with self.assertRaises(ValueError):
            _artifact(budget_w=0.01)  # undercuts r_loc + drift
        with self.assertRaises(ValueError):
            DriftBoundClaim(derivation="x", beta=(0.5, 1.0, 2.0),
                            half_widths=(0.01,) * 3, bu_max=1.2,
                            drift=0.01)  # undercuts bu*h*3
        with self.assertRaises(ValueError):
            _artifact(coloring=(0, 0))  # wrong length

    def test_completeness_and_grade(self):
        a = _artifact()
        self.assertFalse(a.complete)  # frontier open
        self.assertEqual(a.grade, CertificateGrade.SAMPLED_BOUND)
        cell = ((0.0,) * 5, (0.001,) * 5)
        a2 = _artifact(sweep=_sweep(uncovered=(cell,), frontier=True))
        self.assertFalse(a2.complete)
        a3 = _artifact(sweep=_sweep(frontier=True))
        self.assertTrue(a3.complete)
        self.assertEqual(a3.grade, CertificateGrade.SAMPLED_BOUND)
        # rigor arrives only when the sweep arithmetic is rigorous
        rig = SweepReplaySpec(
            sweep="starve.fat_sweep_hulls", wmin=0.025, cell=0.025,
            tax_derivation="rates.certified_rates(mv=True) tax constants",
            boxes_swept=1000, hull_cells=3, frontier_complete=True,
            arithmetic_grade=CertificateGrade.RIGOROUS)
        a4 = _artifact(sweep=rig)
        self.assertEqual(a4.grade, CertificateGrade.RIGOROUS)
        self.assertTrue(a4.ok)


class TestReadiness(unittest.TestCase):
    def _result(self, ok=True):
        return CertificateResult(
            ok, CertificateGrade.SAMPLED_BOUND,
            (Evidence("test", CertificateGrade.SAMPLED_BOUND, ""),),
            reason="test")

    def test_open_frontier_not_ready(self):
        rep = ball_readiness(self._result(), _artifact())
        self.assertFalse(rep.ready)
        self.assertTrue(any("frontier" in m for m in rep.missing))
        self.assertTrue(any("arithmetic" in m for m in rep.missing))
        self.assertEqual(rep.ball_count, 3)
        self.assertIsNotNone(rep.artifact_id)

    def test_source_grade_never_overrides(self):
        # a failed source leg is reported but does not change the
        # artifact's own grade (FAT_TILE_V2_GAP.md)
        rep = ball_readiness(self._result(ok=False), _artifact())
        self.assertFalse(rep.source_ok)
        self.assertEqual(rep.artifact_grade,
                         CertificateGrade.SAMPLED_BOUND)

    def test_rigorous_artifact_ready(self):
        rig = SweepReplaySpec(
            sweep="starve.fat_sweep_hulls", wmin=0.025, cell=0.025,
            tax_derivation="rates.certified_rates(mv=True) tax constants",
            boxes_swept=1000, hull_cells=3, frontier_complete=True,
            arithmetic_grade=CertificateGrade.RIGOROUS)
        rep = ball_readiness(self._result(), _artifact(sweep=rig))
        self.assertTrue(rep.ready)
        self.assertEqual(rep.missing, ())


class TestFatTileEmission(unittest.TestCase):
    """End-to-end through fat_tile with the sweep monkeypatched."""

    def test_emission(self):
        import fat_tile as ft
        import rates
        import starve

        roots = np.array([[0.0] * 5,
                          [0.01, -0.02, 0.03, 0.0, 0.01],
                          list(ORTH_STEP)])
        hulls_c = roots + 0.001
        hulls_r = np.full((3, 5), 1e-4)

        orig_pool = ft._pool_phases
        orig_sweep = starve.fat_sweep_hulls
        orig_rates = rates.certified_rates
        try:
            ft._pool_phases = lambda beta, n_starts=4000: [
                tuple(r) for r in roots]
            starve.fat_sweep_hulls = lambda beta, h, wmin, cell: (
                hulls_c, hulls_r, 3, 1000)
            rates.certified_rates = lambda beta, hv, mv=None: {
                "beta_unit_vec": np.array([1.2, 0.8, 0.5])}
            out = {}
            res = ft.fat_tile((0.5, 1.0, 2.0), 0.01, artifact_out=out)
        finally:
            ft._pool_phases = orig_pool
            starve.fat_sweep_hulls = orig_sweep
            rates.certified_rates = orig_rates

        self.assertIn("artifact", out)
        art = out["artifact"]
        self.assertEqual(len(art.balls), 3)
        self.assertEqual(art.sweep.uncovered, ())
        self.assertEqual(art.check(), [])
        self.assertEqual(art.chi_bound, res.metadata["chi"])
        self.assertAlmostEqual(art.budget_w, res.metadata["w"])
        # round-trips through the canonical bit-pattern form
        again = BallCoverageArtifact.from_dict(art.as_dict())
        self.assertEqual(again.artifact_id, art.artifact_id)
        rep = out["report"]
        self.assertFalse(rep.ready)  # frontier honestly open
        self.assertTrue(any("frontier" in m for m in rep.missing))


if __name__ == "__main__":
    unittest.main()
