"""Degree-2 Taylor models in 3 real variables (db, |db_i| <= h) with a
rigorous remainder — certified curve residuals for the Q-tubes.

TM = c0 + sum_i li db_i + sum_{i<=j} q_ij db_i db_j + [-rem, rem],
coefficients complex128, rem a magnitude bound covering (a) all degree>=3
content over the box, (b) elementary-series truncation, (c) accumulated
coefficient roundoff (Gamma-style, folded in per op; dwarfed by the h^3
terms it accompanies).

Soundness note: the tube curve theta(db) = theta0 + S db + Q[db,db]/2 is
a DEFINITION — S and Q need no certification; the certified object is
max |g(theta(db), beta0+db)| over the box, which this module bounds with
no sampling and no PAD.
"""

import numpy as np

EPS_OP = 3e-15          # per-operation coefficient roundoff allowance
IDX = [(0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2)]


class TM:
    __slots__ = ("c0", "l", "q", "rem", "h")

    def __init__(self, h, c0=0.0, l=None, q=None, rem=0.0):
        self.h = h
        self.c0 = complex(c0)
        self.l = np.zeros(3, complex) if l is None else np.asarray(l, complex)
        self.q = np.zeros(6, complex) if q is None else np.asarray(q, complex)
        self.rem = float(rem)

    @staticmethod
    def var(h, k):
        t = TM(h)
        t.l[k] = 1.0
        return t

    @staticmethod
    def const(h, z):
        return TM(h, c0=z)

    def mags(self):
        return (abs(self.c0), float(np.abs(self.l).sum()) * self.h,
                float(np.abs(self.q).sum()) * self.h * self.h)

    def bound(self):
        m0, m1, m2 = self.mags()
        return m0 + m1 + m2 + self.rem

    def bound_centered(self):
        """max |TM - c0| over the box."""
        _m0, m1, m2 = self.mags()
        return m1 + m2 + self.rem

    def __neg__(self):
        return TM(self.h, -self.c0, -self.l, -self.q, self.rem)

    def __add__(self, o):
        if not isinstance(o, TM):
            o = TM.const(self.h, o)
        return TM(self.h, self.c0 + o.c0, self.l + o.l, self.q + o.q,
                  self.rem + o.rem + EPS_OP * (self.bound() + o.bound()))

    __radd__ = __add__

    def __sub__(self, o):
        return self + (-o)

    def __rsub__(self, o):
        return (-self) + o

    def scale(self, z):
        z = complex(z)
        return TM(self.h, self.c0 * z, self.l * z, self.q * z,
                  self.rem * abs(z) + EPS_OP * self.bound() * abs(z))

    def conj(self):
        return TM(self.h, np.conj(self.c0), np.conj(self.l),
                  np.conj(self.q), self.rem)

    def __mul__(self, o):
        if not isinstance(o, TM):
            return self.scale(o)
        h = self.h
        r = TM(h)
        r.c0 = self.c0 * o.c0
        r.l = self.c0 * o.l + o.c0 * self.l
        # quadratic: c0*q + q*c0 + l x l
        r.q = self.c0 * o.q + o.c0 * self.q
        for n, (i, j) in enumerate(IDX):
            if i == j:
                r.q[n] += self.l[i] * o.l[i]
            else:
                r.q[n] += self.l[i] * o.l[j] + self.l[j] * o.l[i]
        # degree >= 3 spill + remainder cross terms
        _s0, s1, s2 = self.mags()
        _o0, o1, o2 = o.mags()
        spill = s1 * o2 + s2 * o1 + s2 * o2
        r.rem = (spill + self.rem * o.bound() + o.rem * self.bound()
                 + EPS_OP * self.bound() * o.bound() * 4.0)
        return r

    __rmul__ = __mul__

    def invert(self):
        """1/TM; requires the centered part small vs |c0|."""
        z0 = self.c0
        t = (self - z0).scale(1.0 / z0)          # s with |s| = u
        u = t.bound()
        if u > 0.35:
            raise RuntimeError(f"TM invert: |s| = {u:.3f} too large")
        # 1/(1+s) = 1 - s + s^2 + rem(|s|^3/(1-|s|))
        s2 = t * t
        out = (TM.const(self.h, 1.0) - t + s2)
        out.rem += u ** 3 / (1.0 - u)
        return out.scale(1.0 / z0)

    def __truediv__(self, o):
        if not isinstance(o, TM):
            return self.scale(1.0 / o)
        return self * o.invert()

    def csqrt_tm(self):
        """sqrt(TM) around the center (principal, center off the cut)."""
        z0 = self.c0
        w0 = np.sqrt(z0)
        t = (self - z0).scale(1.0 / z0)
        u = t.bound()
        if u > 0.3:
            raise RuntimeError(f"TM csqrt: |u| = {u:.3f} too large")
        # sqrt(1+u) = 1 + u/2 - u^2/8 + rem; |rem| <= 0.0625*|u|^3 for
        # |u| <= 0.3 (next Taylor coefficient 1/16 with tail factor < 1.0)
        s2 = t * t
        out = TM.const(self.h, 1.0) + t.scale(0.5) - s2.scale(0.125)
        out.rem += 0.0625 * u ** 3 / (1.0 - u)
        return out.scale(w0)


def tm_sincos(h, base, k):
    """(sin, cos) of (base + db_k) as TMs: exact deg-2 series in db_k."""
    x = TM.var(h, k)
    x2 = x * x
    sb, cb = np.sin(base), np.cos(base)
    sinx = x                                  # + rem h^3/6
    sinx = TM(h, 0.0, x.l.copy(), x.q.copy(), x.rem + h ** 3 / 6.0)
    cosx = TM.const(h, 1.0) - x2.scale(0.5)
    cosx.rem += h ** 4 / 24.0
    s = sinx.scale(cb) + cosx.scale(sb)
    c = cosx.scale(cb) - sinx.scale(sb)
    return s, c


def tm_exp_i(theta_tm):
    """e^{i t} for a CENTERED deg-2 TM t (c0 folded out by caller):
    1 + it - t^2/2 + rem(|t|^3/6 * e^{|t|})."""
    t = theta_tm
    u = t.bound()
    t2 = t * t
    out = TM.const(t.h, 1.0) + t.scale(1j) - t2.scale(0.5)
    out.rem += u ** 3 / 6.0 * float(np.exp(u))
    return out
