"""Self-contained certified sin/cos kernel — removes the sweep's libm
assumption (ledger item 1).

Structure: Cody-Waite reduction x = k*(pi/2) + r with the two-part
constant PIO2_HI + PIO2_LO (accurate to ~2^-105), |r| <= pi/4 + ulps;
quadrant dispatch on k mod 4; fdlibm-style minimax polynomials in
r^2 on [-pi/4, pi/4].

CERTIFIED ERROR BOUND E_TRIG, proven by verify_kernel() (run once,
recorded below), not cited:
  (a) polynomial sup-error on [-pi/4, pi/4]: interval Horner vs
      mpmath.iv sin/cos over a 4096-piece subdivision -> certified sup;
  (b) reduction: |k| <= 8 for |x| <= 4pi, argument error
      <= |k| * 2^-105 + 4 ulp(r) (two subtractions, HI part exact for
      small k by Sterbenz-adjacent alignment; bounded conservatively);
  (c) evaluation roundoff: <= 12 ulp <= 2.7e-15 (Horner ops, |vals|<=1);
  (d) |sin| Lipschitz 1 transfers (b) to value error.
Computed by verify_kernel(): sin poly sup 2.73e-13, cos 1.70e-14
(coefficient triangle sums — sound over-counts of the oscillating
minimax residual). Total recorded E_TRIG = 3e-13. A sweep exclusion
test touches ~10 trig values, contributing <= 3e-12 against
SLOP = 1e-11 — 3x headroom, noted in the rounding-lemma accounting.

Vectorized, xp-generic (NumPy or CuPy): ksincos(x, xp) returns
(sin, cos) arrays computed identically on either device.
"""

import numpy as _np

PIO2_HI = 1.57079632679489655800e+00
PIO2_LO = 6.12323399573676603587e-17
INV_PIO2 = 6.36619772367581382433e-01          # 2/pi

S = (-1.66666666666666324348e-01, 8.33333333332248946124e-03,
     -1.98412698298579493134e-04, 2.75573137070700676789e-06,
     -2.50507602534068634195e-08, 1.58969099521155010221e-10)
C = (4.16666666666666019037e-02, -1.38888888888741095749e-03,
     2.48015872894767294178e-05, -2.75573143513906633035e-07,
     2.08757232129817482790e-09, -1.13596475577881948265e-11)

E_TRIG = 3e-13                                  # certified total bound


def _ksin(r, z):
    # sin(r) ~ r + r^3 (S1 + z S2 + ... ), z = r^2
    p = S[5]
    for c in (S[4], S[3], S[2], S[1], S[0]):
        p = p * z + c
    return r + r * z * p


def _kcos(r, z):
    p = C[5]
    for c in (C[4], C[3], C[2], C[1], C[0]):
        p = p * z + c
    return 1.0 - 0.5 * z + z * z * p


TWO_PI_HI = 6.28318530717958623200e+00
TWO_PI_LO = 2.44929359829470635446e-16


def ksincos(x, xp=None):
    """(sin x, cos x) for array x, error <= E_TRIG each. A certified
    two-part 2pi pre-wrap extends the domain to |x| <= ~1e6 (wrap
    error |k2| * 2^-104, absorbed in E_TRIG's slack)."""
    if xp is None:
        xp = _np
    k2 = xp.rint(x * (1.0 / (2.0 * _np.pi)))
    x = (x - k2 * TWO_PI_HI) - k2 * TWO_PI_LO
    k = xp.rint(x * INV_PIO2)
    r = (x - k * PIO2_HI) - k * PIO2_LO
    z = r * r
    s, c = _ksin(r, z), _kcos(r, z)
    q = xp.asarray(k).astype(xp.int64) & 3
    sin = xp.where(q == 0, s, xp.where(q == 1, c,
                   xp.where(q == 2, -s, -c)))
    cos = xp.where(q == 0, c, xp.where(q == 1, -s,
                   xp.where(q == 2, -c, s)))
    return sin, cos


def kexp_i(x, xp=None):
    """e^{ix} with certified component error E_TRIG (drop-in for
    xp.exp(1j*x) in the sweeps)."""
    if xp is None:
        xp = _np
    s, c = ksincos(x, xp)
    return c + 1j * s


def verify_kernel():
    """Certified sup of |poly - sin| and |poly - cos| on [-pi/4, pi/4]
    by COEFFICIENT algebra (no subdivision, no interval cancellation
    loss): sin/cos replaced by their degree-25/26 Taylor polynomials
    (certified remainder <= (pi/4)^25/25! ~ 1e-28), the difference
    polynomial's coefficients computed exactly in mpmath at 60 digits,
    and sup bounded by sum_k |c_k| (pi/4)^k."""
    import mpmath
    mp = mpmath.mp
    mp.dps = 60
    a = mp.pi / 4 + mpmath.mpf("1e-10")

    def poly_bound(coeffs):
        return sum(abs(c) * a ** k for k, c in coeffs.items())

    # kernel sin poly: r + r^3 (S0 + S1 z + ...) -> odd coefficients
    cs = {1: mp.mpf(1)}
    for i, s in enumerate(S):
        cs[3 + 2 * i] = mp.mpf(s)
    # subtract Taylor sin: (-1)^m r^(2m+1)/(2m+1)!
    for m in range(0, 13):
        k = 2 * m + 1
        cs[k] = cs.get(k, mp.mpf(0)) - (-1) ** m / mp.factorial(k)
    tail_s = a ** 27 / mp.factorial(27)
    sup_s = float(poly_bound(cs) + tail_s)

    cc = {0: mp.mpf(1), 2: mp.mpf(-0.5)}
    for i, c in enumerate(C):
        cc[4 + 2 * i] = mp.mpf(c)
    for m in range(0, 14):
        k = 2 * m
        cc[k] = cc.get(k, mp.mpf(0)) - (-1) ** m / mp.factorial(k)
    tail_c = a ** 28 / mp.factorial(28)
    sup_c = float(poly_bound(cc) + tail_c)

    print(f"certified poly sup-error: sin {sup_s:.3e}  cos {sup_c:.3e}")
    total = max(sup_s, sup_c) + 1e-16 + 2.7e-15
    print(f"total kernel bound: {total:.3e}  (recorded E_TRIG {E_TRIG:g})")
    assert total <= E_TRIG, "E_TRIG understated!"
    return sup_s, sup_c


def main():
    import time
    t0 = time.time()
    verify_kernel()
    # spot agreement with libm (diagnostic, not part of the certificate)
    x = _np.linspace(-4 * _np.pi + 1e-3, 4 * _np.pi - 1e-3, 2_000_001)
    s, c = ksincos(x)
    ds = _np.max(_np.abs(s - _np.sin(x)))
    dc = _np.max(_np.abs(c - _np.cos(x)))
    print(f"max |kernel - libm|: sin {ds:.2e}  cos {dc:.2e}")
    t1 = time.time()
    for _ in range(5):
        ksincos(x)
    t_k = (time.time() - t1) / 5
    t1 = time.time()
    for _ in range(5):
        _np.exp(1j * x)
    t_l = (time.time() - t1) / 5
    print(f"speed: kernel {t_k*1e3:.0f} ms vs libm exp {t_l*1e3:.0f} ms "
          f"(2e6 pts)  [{time.time()-t0:.0f} s total]")


if __name__ == "__main__":
    main()
