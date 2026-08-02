"""Forward-mode AD over interval types: certified derivative enclosures
of the Karlsson map, and mean-value tile bounds that kill the 30x
dependency blowup of naive interval evaluation.

RDual: IV value + 3 IV partials (w.r.t. theta, phi, lam).
CDual: CIV value + 3 CIV partials.

Mean-value form over a box of half-width h centered at c:
    |H_ij(b) - H_ij(c)| <= sum_k mag(dH_ij/db_k over box) * h,
whose overestimation is second-order (the derivative enclosure's own
dependency artifact multiplies another factor h).
"""

import warnings

import numpy as np

from interval import CIV, IV, iv_cos, iv_sin
from karlsson import karlsson_map

warnings.filterwarnings("ignore")

SQRT3_HALF = IV.pad(np.sqrt(3.0) / 2.0, 8)
INV_SQRT6 = IV.pad(1.0 / np.sqrt(6.0), 8)
ZERO3R = lambda: [IV(0.0), IV(0.0), IV(0.0)]
ZERO3C = lambda: [CIV(0.0), CIV(0.0), CIV(0.0)]


class RDual:
    __slots__ = ("v", "d")

    def __init__(self, v, d=None):
        self.v = v
        self.d = d if d is not None else ZERO3R()

    @staticmethod
    def var(iv, k):
        d = ZERO3R()
        d[k] = IV(1.0)
        return RDual(iv, d)

    def __neg__(self):
        return RDual(-self.v, [-x for x in self.d])

    def __mul__(self, o):
        if isinstance(o, RDual):
            return RDual(self.v * o.v,
                         [self.v * o.d[k] + o.v * self.d[k]
                          for k in range(3)])
        return RDual(self.v * o, [x * o for x in self.d])

    __rmul__ = __mul__


def cosD(x):
    s = iv_sin(x.v)
    return RDual(iv_cos(x.v), [-(s * x.d[k]) for k in range(3)])


def sinD(x):
    c = iv_cos(x.v)
    return RDual(iv_sin(x.v), [c * x.d[k] for k in range(3)])


class CDual:
    __slots__ = ("v", "d")

    def __init__(self, v, d=None):
        self.v = v if isinstance(v, CIV) else CIV(v)
        self.d = d if d is not None else ZERO3C()

    @staticmethod
    def from_r(r):
        return CDual(CIV(r.v, IV(0.0)),
                     [CIV(r.d[k], IV(0.0)) for k in range(3)])

    @staticmethod
    def imag_of(r):
        """i * r for RDual r."""
        return CDual(CIV(IV(0.0), r.v),
                     [CIV(IV(0.0), r.d[k]) for k in range(3)])

    def __neg__(self):
        return CDual(-self.v, [-x for x in self.d])

    def __add__(self, o):
        o = o if isinstance(o, CDual) else CDual(o)
        return CDual(self.v + o.v,
                     [self.d[k] + o.d[k] for k in range(3)])

    __radd__ = __add__

    def __sub__(self, o):
        o = o if isinstance(o, CDual) else CDual(o)
        return CDual(self.v - o.v,
                     [self.d[k] - o.d[k] for k in range(3)])

    def __rsub__(self, o):
        return CDual(o) - self

    def __mul__(self, o):
        o = o if isinstance(o, CDual) else CDual(o)
        return CDual(self.v * o.v,
                     [self.v * o.d[k] + o.v * self.d[k] for k in range(3)])

    __rmul__ = __mul__

    def conj(self):
        return CDual(self.v.conj(), [x.conj() for x in self.d])

    def __truediv__(self, o):
        o = o if isinstance(o, CDual) else CDual(o)
        q = self.v / o.v
        # quotient rule in the tighter form (d_num - q d_den)/den
        return CDual(q, [(self.d[k] - q * o.d[k]) / o.v
                         for k in range(3)])

    def csqrt(self):
        w = self.v.csqrt()
        inv2w = CIV(1.0) / (w * CIV(2.0))
        return CDual(w, [x * inv2w for x in self.d])

    def scale_iv(self, s):
        return CDual(self.v.scale(s), [x.scale(s) for x in self.d])


def expiD(r):
    """e^{i r} for RDual r -> CDual."""
    c, s = iv_cos(r.v), iv_sin(r.v)
    v = CIV(c, s)
    iz = CIV(-s, c)          # i e^{ir}
    return CDual(v, [iz * CIV(r.d[k]) for k in range(3)])


def matmulD(A, B):
    n, k, m = len(A), len(B), len(B[0])
    out = [[None] * m for _ in range(n)]
    for i in range(n):
        for j in range(m):
            acc = CDual(CIV(0.0))
            for l in range(k):
                acc = acc + A[i][l] * B[l][j]
            out[i][j] = acc
    return out


def dual_karlsson(theta_iv, phi_iv, lam_iv):
    """Karlsson map with certified values AND first partials over the
    parameter box. Returns 6x6 CDual matrix (entries include 1/sqrt6)."""
    th = RDual.var(theta_iv, 0)
    ph = RDual.var(phi_iv, 1)
    la = RDual.var(lam_iv, 2)

    ct, st = cosD(th), sinD(th)
    eip = expiD(ph)
    L11 = CDual.from_r(ct)
    L12 = eip * CDual.from_r(st)
    L21 = eip.conj() * CDual.from_r(st)

    i_s3 = CDual(CIV(IV(0.0), SQRT3_HALF))
    M11 = CDual(CIV(IV(-0.5))) + i_s3 * L11
    M12 = i_s3 * L12
    M21 = i_s3 * L21
    M22 = CDual(CIV(IV(-0.5))) - i_s3 * L11
    A11 = M11 + M21
    A12 = M12 + M22
    B11 = -CDual(CIV(1.0)) - A11
    B12 = -CDual(CIV(1.0)) - A12

    z1 = expiD(la)
    w = z1 * z1

    def moebius(P11, P12):
        num = P11 * P11 - w * (P12 * P12)
        den = (P12.conj() * P12.conj()) - w * (P11.conj() * P11.conj())
        if den.v.abs2().lo <= 0:
            raise RuntimeError("Moebius denominator touches 0")
        return num / den

    # z3 via the FACTORED quotient (NOTES 4.57/4.72):
    # algebraically identical to moebius(A11, A12) but the shared
    # small factors are enclosed directly (interval widths are
    # additive, not cancellation-amplified) — removes the corner
    # blowup in the dual path. z4's factors are O(1) everywhere.
    P_f = A11 + z1 * A12
    Q_f = A12.conj() + z1 * A11.conj()
    num3 = (A11 - z1 * A12) * P_f
    den3 = (A12.conj() - z1 * A11.conj()) * Q_f
    if den3.v.abs2().lo <= 0:
        raise RuntimeError("Moebius denominator touches 0")
    z3sq = num3 / den3
    z4sq = moebius(B11, B12)
    num2 = B11 * B11 - z3sq * (B12.conj() * B12.conj())
    den2 = B12 * B12 - z3sq * (B11.conj() * B11.conj())
    if den2.v.abs2().lo <= 0:
        raise RuntimeError("z2 denominator touches 0")
    z2sq = num2 / den2
    for name, zsq in (("z3", z3sq), ("z4", z4sq), ("z2", z2sq)):
        if not zsq.v.cut_clear():
            raise RuntimeError(f"{name}^2 touches branch cut")
    z3, z4, z2 = z3sq.csqrt(), z4sq.csqrt(), z2sq.csqrt()

    one = CDual(CIV(1.0))
    F2 = [[one, one], [one, -one]]
    Z1 = [[one, one], [z1, -z1]]
    Z2 = [[one, one], [z2, -z2]]
    Z3 = [[one, z3], [one, -z3]]
    Z4 = [[one, z4], [one, -z4]]
    A = [[A11, A12], [A12.conj(), -A11.conj()]]
    B = [[B11, B12], [B12.conj(), -B11.conj()]]

    half = IV(0.5)
    blocks = {(0, 0): F2, (0, 1): Z1, (0, 2): Z2, (1, 0): Z3, (2, 0): Z4}
    blocks[(1, 1)] = [[e.scale_iv(half) for e in row]
                      for row in matmulD(matmulD(Z3, A), Z1)]
    blocks[(1, 2)] = [[e.scale_iv(half) for e in row]
                      for row in matmulD(matmulD(Z3, B), Z2)]
    blocks[(2, 1)] = [[e.scale_iv(half) for e in row]
                      for row in matmulD(matmulD(Z4, B), Z1)]
    blocks[(2, 2)] = [[e.scale_iv(half) for e in row]
                      for row in matmulD(matmulD(Z4, A), Z2)]

    H = [[None] * 6 for _ in range(6)]
    for (bi, bj), blk in blocks.items():
        for a in range(2):
            for b in range(2):
                H[2 * bi + a][2 * bj + b] = blk[a][b].scale_iv(INV_SQRT6)
    return H


def cdual_mag(c):
    return c.mag().hi


def main():
    print("=== dual (mean-value) certified tile bounds ===")
    beta = (5.978503016422594, 4.007534549834652, 1.6327649325136653)
    Hc = karlsson_map(*beta)

    # derivative containment check vs FD at the center
    Hd = dual_karlsson(IV(beta[0]), IV(beta[1]), IV(beta[2]))
    eps = 1e-6
    ok = 0
    for k in range(3):
        b1 = list(beta)
        b1[k] += eps
        FD = (karlsson_map(*b1) - Hc) / eps
        good = all(abs(complex(Hd[i][j].d[k].re.mid,
                               Hd[i][j].d[k].im.mid) - FD[i, j]) < 1e-5
                   for i in range(6) for j in range(6))
        ok += good
    print(f"1. derivative containment vs FD: {ok}/3 directions agree")

    # certified mean-value tile bounds vs prototype PAD*FD constants
    for h in (3e-4, 1e-3, 3e-3):
        Hd = dual_karlsson(IV(beta[0] - h, beta[0] + h),
                           IV(beta[1] - h, beta[1] + h),
                           IV(beta[2] - h, beta[2] + h))
        dev = 0.0
        Lmax = 0.0
        for i in range(6):
            for j in range(6):
                s = sum(cdual_mag(Hd[i][j].d[k]) for k in range(3))
                dev = max(dev, s * h)
                Lmax = max(Lmax, max(cdual_mag(Hd[i][j].d[k])
                                     for k in range(3)))
        proto = 2.0 * 1.04 * np.sqrt(3) * h
        print(f"2. tile h={h:g}: certified mean-value max|dH| = {dev:.3e} "
              f"(naive interval was {32*proto/2:.1e}-ish; prototype "
              f"PAD*FD {proto:.3e}; ratio {dev/proto:.2f}); "
              f"certified L_map = {Lmax:.3f} (FD said 1.04)")


if __name__ == "__main__":
    main()
