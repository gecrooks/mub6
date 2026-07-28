"""Certified arc rates for the Beauchamp-Nicoara (B-arc) walks —
removing assumption 3 of the ledger (FD+PAD family rates).

dB/dtheta: dual-AD implementation of the B-N map on the interval
substrate (mpmath.iv transcendentals, CDual arithmetic; the csqrt
branch clearance and denominator nonvanishing are checked by the
arithmetic itself, not assumed).

dK/dtheta: K is implicit (partner basis of MU vectors), so its motion
over a theta-ball is bounded by parametric Krawczyk: krawczyk_verify
certifies a unique root in a box for EVERY Hadamard in the hslop-ball
around B; with hslop >= sup|dB/dtheta| * delta the box's tight
enclosure radius bounds each column's phase motion, hence its entry
drift (|e^{i a} - e^{i b}| <= |a-b|, scaled 1/sqrt6).

The triple sweep's ball slop is then rate_B * delta + K entry drift —
no FD, no PAD.

STATUS (2026-07-28): dual_bn is DONE and tight (sup|dB/dth| enclosure
0.567 over a 4e-3 ball at theta=1.6 vs FD 0.57 — no dependency
blowup). certified_k_drift works for well-conditioned columns but K
has ILL-CONDITIONED columns (||J^-1|| ~ 5e2 at theta=1.6): the
q-contraction cannot dominate there, and the certified treatment is
the 1-parameter frame-split continuation (strong 4x4 block contracts
with q ~ 0.02; the weak direction needs the signed interval dg with
cancellation via Y @ dg_iv, plus a 1-dim crossing argument — the
valley machinery specialized to one parameter). Remaining ledger
item 3 narrows to exactly that piece.
"""

import warnings

import numpy as np

from certify import krawczyk_verify
from dual import CDual, RDual, cdual_mag, expiD
from interval import CIV, IV

warnings.filterwarnings("ignore")

SQ6 = np.sqrt(6.0)


def dual_bn(theta_iv):
    """6x6 CDual matrix of B(theta)/1 over the theta interval (entries
    include the 1/sqrt6 normalization), with certified d/dtheta."""
    th = RDual.var(theta_iv, 0)
    y = expiD(th)
    one = CDual(CIV(1.0))
    two = CDual(CIV(2.0))
    y2 = y * y
    z = (one + two * y - y2) / (y * (-one + two * y + y2))
    disc = one + two * y + two * y * y2 + y2 * y2
    sq2 = CDual(CIV(IV.pad(np.sqrt(2.0), 8)))
    x = (one + two * y + y2 - sq2 * disc.csqrt()) / (one + two * y - y2)
    t = x * y * z
    xb, yb, zb, tb = x.conj(), y.conj(), z.conj(), t.conj()
    c = {1: one, -1: -one}
    rows = [
        [c[1], c[1], c[1], c[1], c[1], c[1]],
        [c[1], c[-1], -xb, -y, y, xb],
        [c[1], -x, c[1], y, zb, -tb],
        [c[1], -yb, yb, c[-1], -tb, tb],
        [c[1], yb, z, -t, c[1], -xb],
        [c[1], x, -t, t, -x, c[-1]],
    ]
    inv6 = CDual(CIV(IV.pad(1.0 / SQ6, 8)))
    return [[inv6 * e for e in row] for row in rows]


def certified_bn_rate(theta, delta):
    """Certified sup over [theta-delta, theta+delta] of the max entry
    magnitude of dB/dtheta (B including 1/sqrt6)."""
    Bd = dual_bn(IV(theta - delta, theta + delta))
    return max(cdual_mag(Bd[i][j].d[0])
               for i in range(6) for j in range(6))


def certified_k_drift(B, K, theta, delta, slop_extra=1e-11):
    """Certified entry-drift bound for K's columns over the theta-ball
    by certified continuation: tight Krawczyk at the anchor (existence
    + uniqueness), then |du/dtheta| <= ||J^{-1}|| ||dg/dtheta|| with
    sigma_min(J) lower-bounded over the swept tube (Hessian remainder
    r/3 + sum r/18 + family-drift L_H_J * hs). Self-consistent tube
    radius via two fixed-point passes. Raises on failure."""
    from certify import SLOP, _g_and_J
    Bd = dual_bn(IV(theta - delta, theta + delta))
    e1 = np.array([[cdual_mag(Bd[i][k].d[0]) for k in range(6)]
                   for i in range(6)])
    col1 = e1.sum(axis=0)                       # sum_i |dB_ik/dth|
    rate_B = float(e1.max())
    hs = rate_B * delta + slop_extra
    inv6 = 1.0 / SQ6
    worst_motion = 0.0
    for j in range(K.shape[1]):
        u = K[:, j] * np.exp(-1j * np.angle(K[0, j]))
        th = np.angle(u * SQ6)[1:]
        ok, c, r_encl, _R = krawczyk_verify(B, th, np.full(5, 2e-4),
                                            hslop=1e-11)
        if not ok:
            raise RuntimeError(f"K column {j}: anchor Krawczyk failed")
        _g, J = _g_and_J(B, c)
        Y = np.linalg.inv(J)
        hmagB = float(max(Bd[i][k].v.mag().hi
                          for i in range(6) for k in range(6)))
        # family drift of J per unit theta (rates.py dJ structure with
        # the B-N entry rates): |dJ_kl/dth| <= 2[(1/sqrt6) col1_k
        # hmagB/sqrt6 + e1[l,k]/sqrt6]
        dJ_fam = 2.0 * ((inv6 * col1[1:6, None] * hmagB * inv6)
                        + e1[1:6, 1:6].T / SQ6)
        r_tube = float(r_encl.max())
        for _ in range(3):
            s_enc = inv6 + inv6 * (5.0 * r_tube) + hs
            # |dg_k/dth| <= 2 |s_k|_enc (1/sqrt6) sum_i |dB_ik/dth|
            dg = 2.0 * min(s_enc, 1.0) * inv6 * col1[1:6]
            RJ_col = (r_tube / 3.0 + 5.0 * r_tube / 18.0
                      + dJ_fam.max() * delta + SLOP)
            q = float(np.max(np.abs(np.eye(5) - Y @ J).sum(axis=1))
                      + np.abs(Y).sum(axis=1).max() * RJ_col)
            if q >= 0.7:
                raise RuntimeError(f"K column {j}: contraction q = "
                                   f"{q:.2f} >= 0.7")
            rate_u = float(np.max(np.abs(Y) @ dg)) / (1.0 - q)
            r_tube = float(r_encl.max()) + 1.1 * rate_u * delta
        worst_motion = max(worst_motion, r_tube)
    return worst_motion / SQ6, worst_motion


def certified_pair_slop(theta, delta, B, K):
    """Ball slop for certified_triple_sweep over the theta-ball:
    sup entry drift of B plus of K — certified, no PAD."""
    rate_B = certified_bn_rate(theta, delta)
    dK, _ = certified_k_drift(B, K, theta, delta)
    return rate_B * delta + dK


def main():
    import time
    from layer3 import build_triple
    theta = 1.6
    delta = 4e-3
    t0 = time.time()
    B, K, _p, _b = build_triple(theta)
    rate_B = certified_bn_rate(theta, delta)
    dK, wph = certified_k_drift(B, K, theta, delta)
    # compare with the FD estimate the walks use
    from mub import beauchamp_nicoara
    eps = 1e-4
    B1 = beauchamp_nicoara(theta + eps)
    fd = np.max(np.abs(B1 - B)) / eps
    print(f"certified sup|dB/dth| = {rate_B:.4f}  (FD center {fd:.4f})")
    print(f"certified K drift over delta={delta:g}: {dK:.2e} "
          f"(phase radius {wph:.2e})")
    print(f"certified pair slop = {rate_B*delta + dK:.3e}  vs FD+PAD "
          f"~ {1.25*(fd+dK/delta)*delta:.3e}")
    print(f"[{time.time()-t0:.0f} s]")


if __name__ == "__main__":
    main()
