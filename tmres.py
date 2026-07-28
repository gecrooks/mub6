"""Certified curve residuals via Taylor models: R(h) bounding
max |g(theta0 + S db + Q[db,db]/2, beta0 + db)| over |db|_inf <= h,
with no sampling and no PAD. Replaces `curve_residual` + PAD in the
Q-tube certificates (roadmap R4 tail item)."""

import warnings

import numpy as np

from tm import TM, tm_exp_i, tm_sincos

warnings.filterwarnings("ignore")

F2 = np.array([[1, 1], [1, -1]], dtype=complex)


def tm_karlsson(beta, h):
    """6x6 matrix of TMs enclosing H(beta + db) over |db|_inf <= h.
    Mirrors karlsson.karlsson_map / dual.dual_karlsson."""
    th0, ph0, la0 = beta
    st, ct = tm_sincos(h, th0, 0)
    sp, cp = tm_sincos(h, ph0, 1)
    sl, cl = tm_sincos(h, la0, 2)
    eip = cp + sp.scale(1j)
    eim = cp - sp.scale(1j)
    i_s3 = 1j * np.sqrt(3.0) / 2.0
    L11 = ct
    L12 = eip * st
    L21 = eim * st
    M11 = TM.const(h, -0.5) + L11.scale(i_s3)
    M12 = L12.scale(i_s3)
    M21 = L21.scale(i_s3)
    M22 = TM.const(h, -0.5) - L11.scale(i_s3)
    A11 = M11 + M21
    A12 = M12 + M22
    B11 = -TM.const(h, 1.0) - A11
    B12 = -TM.const(h, 1.0) - A12

    z1 = cl + sl.scale(1j)
    w = z1 * z1

    def moebius(P11, P12):
        num = P11 * P11 - w * (P12 * P12)
        den = (P12.conj() * P12.conj()) - w * (P11.conj() * P11.conj())
        return num / den

    z3 = moebius(A11, A12).csqrt_tm()
    z4 = moebius(B11, B12).csqrt_tm()
    z3sq = z3 * z3
    num2 = B11 * B11 - z3sq * (B12.conj() * B12.conj())
    den2 = B12 * B12 - z3sq * (B11.conj() * B11.conj())
    z2 = (num2 / den2).csqrt_tm()

    one = TM.const(h, 1.0)
    Z1 = [[one, one], [z1, -z1]]
    Z2 = [[one, one], [z2, -z2]]
    Z3 = [[one, z3], [one, -z3]]
    Z4 = [[one, z4], [one, -z4]]
    A = [[A11, A12], [A12.conj(), -A11.conj()]]
    B = [[B11, B12], [B12.conj(), -B11.conj()]]
    F2t = [[one, one], [one, -one]]

    def mm(X, Y, half=False):
        out = [[None, None], [None, None]]
        for i in range(2):
            for j in range(2):
                acc = X[i][0] * Y[0][j] + X[i][1] * Y[1][j]
                out[i][j] = acc.scale(0.5) if half else acc
        return out

    blocks = {(0, 0): F2t, (0, 1): Z1, (0, 2): Z2, (1, 0): Z3, (2, 0): Z4,
              (1, 1): mm(mm(Z3, A), Z1, half=True),
              (1, 2): mm(mm(Z3, B), Z2, half=True),
              (2, 1): mm(mm(Z4, B), Z1, half=True),
              (2, 2): mm(mm(Z4, A), Z2, half=True)}
    inv6 = 1.0 / np.sqrt(6.0)
    H = [[None] * 6 for _ in range(6)]
    for (bi, bj), blk in blocks.items():
        for a in range(2):
            for b in range(2):
                H[2 * bi + a][2 * bj + b] = blk[a][b].scale(inv6)
    return H


def certified_curve_residual(beta, h, th0, S, Q, Htm=None):
    """Certified bound on max_k |g_k| along the Q-curve over the tile.
    Htm may be shared across roots of the same tile."""
    if Htm is None:
        Htm = tm_karlsson(beta, h)
    inv6 = 1.0 / np.sqrt(6.0)
    # u_j(db) = e^{i theta_j(db)} / sqrt6, theta centered at th0
    us = [TM.const(h, inv6)]
    for j in range(5):
        t = TM(h)
        t.l = S[j].astype(complex)
        t.q = np.array([Q[j, i, k] * (1.0 if i == k else 1.0)
                        for (i, k) in [(0, 0), (0, 1), (0, 2),
                                       (1, 1), (1, 2), (2, 2)]],
                       complex)
        # off-diagonal q entries carry both (i,k) and (k,i) halves
        t.q[1] = Q[j, 0, 1]
        t.q[2] = Q[j, 0, 2]
        t.q[4] = Q[j, 1, 2]
        t.q[0] = 0.5 * Q[j, 0, 0]
        t.q[3] = 0.5 * Q[j, 1, 1]
        t.q[5] = 0.5 * Q[j, 2, 2]
        us.append(tm_exp_i(t).scale(np.exp(1j * th0[j]) * inv6))
    worst = 0.0
    for k in range(1, 6):
        s = TM.const(h, 0.0)
        for j in range(6):
            s = s + Htm[j][k].conj() * us[j]
        g = s * s.conj() - (1.0 / 6.0)
        # |g| over the box: coefficient magnitudes + remainder (g is real
        # up to roundoff; bound covers both parts)
        worst = max(worst, g.bound())
    return worst
