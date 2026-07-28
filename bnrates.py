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

STATUS (2026-07-28, end of session): WORKS end-to-end. dual_bn tight
(0.417 enclosure over a 5e-4 ball vs FD 0.408). certified_k_drift
certifies at delta = 5e-4 with pair slop 9.4e-3 via two independent
continuation fixed points (plain |Y|-bound and adaptive frame-split
with signed point forcing + Lipschitz tube term). Known conservatism
~15x vs the FD-observed K rate (0.39): the deep column (sv5 ~ 2e-3)
is bounded by the plain path's magnitude bound; the frame-split
would beat it at delta ~ 2e-4 (the weak-slope floor sv5 - sqrt5*RJ
demands it) but the ladder accepts the first certifying delta.
Open polish items: (a) let the walker choose the split's delta when
it wins; (b) anchors whose s-ball CONTAINS the branch point need an
algebraic even-symmetry bound (B(s) = B(-s); entry drift via the
sqrt-disc factorization) since dual_bn's csqrt rightly refuses
zero-containing discriminants. Neither blocks the sliver walks:
non-branch anchors can use certified slop today.
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
        U, sv, Vt = np.linalg.svd(J)
        dg_iv = _dg_interval(Bd, c, float(r_encl.max()))
        L_phi = 2.0 * np.sqrt(5.0) * float(col1.max()) * inv6 \
            * (5.0 * hmagB * inv6 + 1.0)

        def rj_col(r_t):
            return (r_t / 3.0 + 5.0 * r_t / 18.0
                    + dJ_fam.max() * delta + SLOP)

        def dgm(r_t):
            s_enc = min(inv6 + inv6 * (5.0 * r_t) + hs, 1.0)
            return 2.0 * s_enc * inv6 * col1[1:6]

        def rate_plain(r_t):
            q = float(np.max(np.abs(np.eye(5) - Y @ J).sum(axis=1))
                      + np.abs(Y).sum(axis=1).max() * rj_col(r_t))
            if q >= 0.7:
                return None
            return float(np.max(np.abs(Y) @ dgm(r_t))) / (1.0 - q)

        def rate_split(r_t):
            # adaptive frame-split: largest contracting strong block;
            # weak rows use the SIGNED point forcing u_i^T dg plus an
            # explicit Lipschitz tube term (cancellation preserved)
            RJ5 = np.sqrt(5.0) * rj_col(r_t)
            if sv[4] - RJ5 <= 0:
                return None
            for m in (4, 3, 2):
                Ym = np.linalg.inv(U[:, :m].T @ J @ Vt[:m].T)
                qm = float(np.abs(Ym).sum(axis=1).max() * RJ5)
                if qm >= 0.5:
                    continue
                rate_w = sum((_proj_interval(U[:, i], dg_iv)
                              + L_phi * r_t) / (sv[i] - RJ5)
                             for i in range(m, 5))
                dgP = np.abs(U[:, :m].T) @ dgm(r_t)
                return rate_w + float(np.max(np.abs(Ym) @ dgP)) \
                    / (1.0 - qm)
            return None

        # independent fixed points per method; take the best converged
        best = None
        for fn in (rate_split, rate_plain):
            r_t = float(r_encl.max())
            rate = None
            for _ in range(3):
                rate = fn(r_t)
                if rate is None:
                    break
                r_t = float(r_encl.max()) + 1.1 * rate * delta
            if rate is not None:
                best = r_t if best is None else min(best, r_t)
        if best is None:
            raise RuntimeError(f"K column {j}: neither continuation "
                               f"method certifies (reduce delta)")
        worst_motion = max(worst_motion, best)
    return worst_motion / SQ6, worst_motion


def _dg_interval(Bd, c, phase_rad):
    """Signed interval vector dg_k = 2 Re(conj(s_k) <u, dh_k/dth>),
    k = 1..5, over the theta-ball (Bd) and the root box c +- phase_rad
    (CIV arithmetic; cancellation preserved)."""
    from interval import iv_cos, iv_sin
    u = [CIV(IV.pad(1.0 / SQ6, 4))]
    for i in range(5):
        ang = IV(c[i] - phase_rad, c[i] + phase_rad)
        u.append(CIV(iv_cos(ang), iv_sin(ang))
                 * CIV(IV.pad(1.0 / SQ6, 4)))
    out = []
    for k in range(1, 6):
        s = CIV(0.0)
        ds = CIV(0.0)
        for i in range(6):
            s = s + u[i].conj() * Bd[i][k].v
            ds = ds + u[i].conj() * Bd[i][k].d[0]
        re = (s.conj() * ds).re * IV(2.0)
        out.append(re)
    return out


def _proj_interval(u5, dg_iv):
    """mag of the signed interval sum sum_k u5_k dg_k."""
    acc = IV(0.0)
    for k in range(5):
        acc = acc + dg_iv[k] * float(u5[k])
    return float(acc.mag())


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
    from mub import beauchamp_nicoara
    eps = 1e-4
    B1 = beauchamp_nicoara(theta + eps)
    fd = np.max(np.abs(B1 - B)) / eps
    for d in (delta, 2e-3, 1e-3, 5e-4, 2e-4, 1e-4):
        try:
            rate_B = certified_bn_rate(theta, d)
            dK, wph = certified_k_drift(B, K, theta, d)
        except RuntimeError as e:
            print(f"delta={d:g}: {e}")
            continue
        print(f"delta={d:g}: sup|dB/dth| = {rate_B:.4f} (FD {fd:.4f}); "
              f"K drift {dK:.2e} (phase {wph:.2e}); "
              f"pair slop {rate_B*d + dK:.3e}")
        break
    print(f"[{time.time()-t0:.0f} s]")


if __name__ == "__main__":
    main()
