"""Smoke tests — fast invariants of the certificate substrate.

Run: .venv/bin/python test_smoke.py   (~2 min, no network)

These are NOT the numerical campaigns; they pin the load-bearing
invariants a refactor could silently break: kernel error bounds,
interval containment, map ports agreeing, fail-closed paths, and
one representative tile of each kind.
"""

import unittest

import numpy as np


class TestKernels(unittest.TestCase):
    def test_trig_kernel_bound(self):
        from trig_kernel import kexp_i, E_TRIG
        x = np.linspace(-np.pi, np.pi, 2001)
        err = np.abs(kexp_i(x[None, :]) - np.exp(1j * x)).max()
        self.assertLess(err, E_TRIG)

    def test_stable_map_matches_naive_generic(self):
        from karlsson import karlsson_map
        from karlsson_stable import stable_karlsson
        H1 = karlsson_map(0.3, 1.0, 2.0)
        H2 = stable_karlsson(0.3, 1.0, 2.0)
        self.assertLess(np.abs(H1 - H2).max(), 1e-12)

    def test_stable_map_is_hadamard_at_corner(self):
        from karlsson_stable import stable_karlsson
        H = stable_karlsson(1e-5, np.pi / 3, np.pi / 3)
        self.assertLess(np.abs(np.abs(H) - 1 / np.sqrt(6)).max(),
                        1e-10)
        self.assertLess(np.abs(H.conj().T @ H - np.eye(6)).max(),
                        1e-9)

    def test_interval_map_contains_float(self):
        from interval import IV
        from ivkarlsson import iv_karlsson
        from karlsson import karlsson_map
        H, _diag = iv_karlsson(IV(0.3), IV(1.0), IV(2.0))
        Hf = karlsson_map(0.3, 1.0, 2.0)
        for a in range(6):
            for b in range(6):
                civ = H[a][b]
                self.assertTrue(civ.re.contains(Hf[a, b].real))
                self.assertTrue(civ.im.contains(Hf[a, b].imag))

    def test_iv_stable_z3sq_tight_and_contains(self):
        # reference = 40-digit truth (the enclosure is now TIGHTER
        # than the float kernel's own error, so float values fall
        # outside it — correctly)
        from interval import IV
        from ivstable import iv_stable_z3sq
        from mpmath import mp, mpc, sqrt, exp, cos, sin
        z, dm = iv_stable_z3sq(IV(1e-4), IV(np.pi / 3),
                               IV(np.pi / 3))
        with mp.workdps(40):
            I = mpc(0, 1)
            th, ph, la = (mp.mpf(1e-4), mp.mpf(np.pi / 3),
                          mp.mpf(np.pi / 3))
            c, s = cos(th), sin(th)
            A11 = mp.mpf(-1) / 2 + I * (sqrt(3) / 2) * (
                c + exp(-I * ph) * s)
            A12 = mp.mpf(-1) / 2 + I * (sqrt(3) / 2) * (
                exp(I * ph) * s - c)
            w = exp(2 * I * la)
            cj = lambda x: x.conjugate()
            ref = ((A11 ** 2 - w * A12 ** 2)
                   / (cj(A12) ** 2 - w * cj(A11) ** 2))
        self.assertTrue(z.re.contains(float(ref.real)))
        self.assertTrue(z.im.contains(float(ref.imag)))
        self.assertLess(max(z.re.width, z.im.width), 1e-12)
        self.assertGreater(dm, 0)


class TestFailClosed(unittest.TestCase):
    def test_coarse_empty_enumeration_fails_even_without_survivors(self):
        from coarse_tile import _coarse_result
        result = _coarse_result(
            n_roots=0, bound=0, coverage_complete=True, uncovered=0,
        )
        self.assertFalse(result)
        self.assertIn("empty enumeration", result.reason)

    def test_coarse_result_grade_cannot_claim_rigorous(self):
        from certificate_result import (
            CertificateGrade, CertificateResult, Evidence,
        )
        result = CertificateResult.from_evidence(
            True,
            [Evidence("chain cache", CertificateGrade.EXPERIMENTAL)],
        )
        self.assertFalse(result.accepted_at(CertificateGrade.RIGOROUS))

    def test_collar_tile_empty_pool_fails(self):
        from collar_tile import collar_tile
        ok = collar_tile(0.005, 0.01, 1.0, 2.0, hf=5e-3,
                         adjacency="signed", pool=np.zeros((0, 5)))
        self.assertFalse(ok)

    def test_chain_step_rejects_failed_anchor(self):
        from coarse_chain import chain_step
        with self.assertRaises(RuntimeError):
            chain_step({"ok": False, "h": 3e-3}, (1.0, 1.0, 1.0))


class TestTiles(unittest.TestCase):
    def test_fully_rigorous_tile_with_global_coverage(self):
        from certificate_result import CertificateGrade
        from rigor_tile import fully_rigorous_signed_tile
        b0 = (5.978503016422594, 4.007534549834652,
              1.6327649325136653)
        h = 3e-4
        result = fully_rigorous_signed_tile(
            b0[0] - h, b0[0] + h, b0[1], b0[2],
            hf=h, hf3=h, verbose=False,
        )
        self.assertTrue(result)
        self.assertEqual(result.grade, CertificateGrade.RIGOROUS)
        coverage = next(
            item for item in result.dependencies
            if item.name == "enumeration-coverage"
        )
        self.assertIn("global sweep closed", coverage.detail)

    def test_signed_collar_tile_generic(self):
        from collar_tile import collar_tile
        ok = collar_tile(0.005, 0.01, 2.041236, float(np.pi),
                         hf=5e-3, adjacency="signed", hf3=2.5e-4)
        self.assertTrue(ok)

    def test_signed_bulk_tile(self):
        from collar_tile import collar_tile
        ok = collar_tile(0.5, 0.55, 1.0, 2.0, hf=5e-3,
                         adjacency="signed", hf3=5e-3)
        self.assertTrue(ok)

    def test_face_tile_generic(self):
        from facewalk import face_tile
        ok, nb, m = face_tile(1.0, 2.0, hf=0.05, th_tube=0.0025)
        self.assertTrue(ok)
        self.assertGreater(m, 0.05)


class TestCertifiedS(unittest.TestCase):
    def test_certified_S_encloses_central_diff(self):
        # reference = independent central differences of the
        # polished root (root_data2's internal continuation
        # carries an ~1e-4 systematic — measured 2026-08-01; the
        # certified path matches central differences to ~1e-9)
        from certpair import certified_S
        from collar_tile import _pool_phases, _polish
        from karlsson import karlsson_map
        beta = (0.05, 1.0, 2.0)
        ph = _pool_phases(beta, n_starts=1500)
        th0 = ph[0]
        S0, err, c, r = certified_S(beta, th0)
        d = 1e-6
        Sd = np.zeros((5, 3))
        for l in range(3):
            bp, bm = list(beta), list(beta)
            bp[l] += d
            bm[l] -= d
            tp = _polish(karlsson_map(*bp), th0)
            tm = _polish(karlsson_map(*bm), th0)
            Sd[:, l] = ((tp - tm + np.pi) % (2 * np.pi)
                        - np.pi) / (2 * d)
        self.assertLess(np.abs(S0 - Sd).max(), err + 1e-6)


if __name__ == "__main__":
    unittest.main(verbosity=2)
