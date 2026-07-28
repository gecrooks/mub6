"""Directed-rounding interval arithmetic (scalar, prototype checker).

Real intervals [lo, hi] with outward rounding via np.nextafter; complex
intervals as axis-aligned rectangles (re, im).

R8: transcendentals (sin, cos, atan2) are backed by mpmath.iv —
guaranteed arbitrary-precision interval enclosures (the pure-Python
Arb equivalent) at 100-bit precision, converted outward to float64
with 1-ulp padding. NO libm assumption remains in the certified-
constants path. sqrt is IEEE-754 correctly rounded by specification
(padded outward anyway), so it carries no assumption either. The only
remaining libm dependence in the proof is the SWEEP's pointwise
float evaluations (covered by SLOP + the static rounding lemma, with
the faithful-libm assumption now scoped to that single use; the
C/GPU kernel with own argument reduction is the documented endgame).

Complex sqrt uses the principal branch and REFUSES rectangles touching
the branch cut (negative real axis) or containing 0 — the caller must
verify clearance, which is exactly the certificate's branch-tracking
obligation.
"""

import math

import mpmath
import numpy as np

ULPS = 4
_MPIV = mpmath.iv
_MPIV.prec = 100


def _mp_endpoints(r):
    """Outward float64 endpoints of an mpmath.iv result."""
    return (_dn(float(mpmath.mpf(r.a)), 1), _up(float(mpmath.mpf(r.b)), 1))


def _up(x, n=1):
    for _ in range(n):
        x = np.nextafter(x, np.inf)
    return float(x)


def _dn(x, n=1):
    for _ in range(n):
        x = np.nextafter(x, -np.inf)
    return float(x)


class IV:
    __slots__ = ("lo", "hi")

    def __init__(self, lo, hi=None):
        if hi is None:
            hi = lo
        if lo > hi:
            raise ValueError("empty interval")
        self.lo = float(lo)
        self.hi = float(hi)

    @staticmethod
    def pad(x, ulps=ULPS):
        return IV(_dn(x, ulps), _up(x, ulps))

    def __repr__(self):
        return f"[{self.lo:.17g}, {self.hi:.17g}]"

    @property
    def width(self):
        return self.hi - self.lo

    @property
    def mid(self):
        return 0.5 * (self.lo + self.hi)

    def contains(self, x):
        return self.lo <= x <= self.hi

    def __neg__(self):
        return IV(-self.hi, -self.lo)

    def __add__(self, o):
        o = o if isinstance(o, IV) else IV(o)
        return IV(_dn(self.lo + o.lo), _up(self.hi + o.hi))

    __radd__ = __add__

    def __sub__(self, o):
        o = o if isinstance(o, IV) else IV(o)
        return IV(_dn(self.lo - o.hi), _up(self.hi - o.lo))

    def __rsub__(self, o):
        return IV(o) - self

    def __mul__(self, o):
        o = o if isinstance(o, IV) else IV(o)
        ps = (self.lo * o.lo, self.lo * o.hi, self.hi * o.lo, self.hi * o.hi)
        return IV(_dn(min(ps)), _up(max(ps)))

    __rmul__ = __mul__

    def __truediv__(self, o):
        o = o if isinstance(o, IV) else IV(o)
        if o.lo <= 0.0 <= o.hi:
            raise ZeroDivisionError("interval divisor contains 0")
        ps = (self.lo / o.lo, self.lo / o.hi, self.hi / o.lo, self.hi / o.hi)
        return IV(_dn(min(ps)), _up(max(ps)))

    def __rtruediv__(self, o):
        return IV(o) / self

    def sq(self):
        m = self * self
        if self.lo <= 0.0 <= self.hi:
            return IV(0.0, m.hi)
        return IV(max(m.lo, 0.0), m.hi)

    def sqrt(self):
        if self.lo < 0:
            raise ValueError("sqrt of interval with negative part")
        return IV(max(_dn(math.sqrt(self.lo), ULPS), 0.0),
                  _up(math.sqrt(self.hi), ULPS))

    def mig(self):
        """min |x| over the interval."""
        if self.lo <= 0.0 <= self.hi:
            return 0.0
        return min(abs(self.lo), abs(self.hi))

    def mag(self):
        return max(abs(self.lo), abs(self.hi))


def iv_sin(x):
    """Certified interval sine via mpmath.iv (extremum logic and pi are
    internal to the guaranteed enclosure; no libm, no float-pi k-range)."""
    lo, hi = _mp_endpoints(_MPIV.sin(_MPIV.mpf([x.lo, x.hi])))
    return IV(max(lo, -1.0), min(hi, 1.0))


def iv_cos(x):
    """Certified interval cosine via mpmath.iv."""
    lo, hi = _mp_endpoints(_MPIV.cos(_MPIV.mpf([x.lo, x.hi])))
    return IV(max(lo, -1.0), min(hi, 1.0))


class CIV:
    """Complex interval as an axis-aligned rectangle."""
    __slots__ = ("re", "im")

    def __init__(self, re, im=None):
        if isinstance(re, complex):
            assert im is None
            self.re = IV(re.real)
            self.im = IV(re.imag)
            return
        self.re = re if isinstance(re, IV) else IV(re)
        self.im = (im if isinstance(im, IV) else IV(im)) \
            if im is not None else IV(0.0)

    @staticmethod
    def from_float(z, ulps=ULPS):
        return CIV(IV(_dn(z.real, ulps), _up(z.real, ulps)),
                   IV(_dn(z.imag, ulps), _up(z.imag, ulps)))

    def __repr__(self):
        return f"({self.re} + i{self.im})"

    @property
    def width(self):
        return max(self.re.width, self.im.width)

    def contains(self, z):
        return self.re.contains(z.real) and self.im.contains(z.imag)

    def __neg__(self):
        return CIV(-self.re, -self.im)

    def __add__(self, o):
        o = o if isinstance(o, CIV) else CIV(o)
        return CIV(self.re + o.re, self.im + o.im)

    __radd__ = __add__

    def __sub__(self, o):
        o = o if isinstance(o, CIV) else CIV(o)
        return CIV(self.re - o.re, self.im - o.im)

    def __rsub__(self, o):
        return CIV(o) - self

    def __mul__(self, o):
        o = o if isinstance(o, CIV) else CIV(o)
        return CIV(self.re * o.re - self.im * o.im,
                   self.re * o.im + self.im * o.re)

    __rmul__ = __mul__

    def conj(self):
        return CIV(self.re, -self.im)

    def abs2(self):
        return self.re.sq() + self.im.sq()

    def __truediv__(self, o):
        o = o if isinstance(o, CIV) else CIV(o)
        d = o.abs2()
        n = self * o.conj()
        return CIV(n.re / d, n.im / d)

    def scale(self, s):
        return CIV(self.re * s, self.im * s)

    def mag(self):
        return IV(math.sqrt(self.re.mig() ** 2 + self.im.mig() ** 2),
                  _up(math.sqrt(self.re.mag() ** 2 + self.im.mag() ** 2),
                      ULPS))

    def cut_clear(self):
        """True if the rectangle avoids the branch cut (-inf, 0]."""
        return self.re.lo > 0 or self.im.lo > 0 or self.im.hi < 0

    def arg(self):
        """Interval argument; requires cut_clear (arg then continuous and
        corner-extremal on the rectangle)."""
        if not self.cut_clear():
            raise ValueError("argument: rectangle touches branch cut")
        los, his = zip(*[_mp_endpoints(_MPIV.atan2(_MPIV.mpf(im),
                                                   _MPIV.mpf(re)))
                         for re in (self.re.lo, self.re.hi)
                         for im in (self.im.lo, self.im.hi)])
        return IV(min(los), max(his))

    def csqrt(self):
        """Principal square root; requires cut clearance."""
        r = self.mag().sqrt()
        half = self.arg() * IV(0.5)
        return CIV(r * iv_cos(half), r * iv_sin(half))


def civ_matmul(A, B):
    """(list-of-lists) complex-interval matrix product."""
    n, k, m = len(A), len(B), len(B[0])
    out = [[CIV(0.0) for _ in range(m)] for _ in range(n)]
    for i in range(n):
        for j in range(m):
            acc = CIV(0.0)
            for l in range(k):
                acc = acc + A[i][l] * B[l][j]
            out[i][j] = acc
    return out


def iv_exp_i(x):
    """e^{ix} for interval x."""
    return CIV(iv_cos(x), iv_sin(x))
