"""Certified (interval) stable z3^2 kernel — the interval-grade
half of the corner fix (NOTES_LP_BRIDGE 4.56-4.57).

Same exact factorization as karlsson_stable.py, in the certified
substrate: z3^2 = (bignum * P) / (bigden * Q), with P, Q evaluated
from the cancellation-free bracket. Removes the eps/theta^3
interval blowup of iv_karlsson's naive Moebius quotient at the
(pi/3, pi/3) corner; off-corner behavior unchanged in kind.

Constants follow house style (IV.pad of float constants: float pi/6
is within 1 ulp of true pi/6, pad 8 ulps encloses it — same pattern
as SQRT3_HALF in ivkarlsson.py).
"""

import numpy as np

from interval import CIV, IV, iv_cos, iv_exp_i, iv_sin

import mpmath as _mpm
_MPIV = _mpm.iv
_MPIV.prec = 100


def _iv_const(x_mp):
    lo = float(_mpm.mpf(x_mp.a) if hasattr(x_mp, 'a') else x_mp)
    hi = float(_mpm.mpf(x_mp.b) if hasattr(x_mp, 'b') else x_mp)
    if lo > hi:
        lo, hi = hi, lo
    return IV(np.nextafter(lo, -np.inf), np.nextafter(hi, np.inf))


SQRT3_HALF = _iv_const(_MPIV.sqrt(_MPIV.mpf(3)) / 2)
SQRT3 = _iv_const(_MPIV.sqrt(_MPIV.mpf(3)))
PI_6 = _iv_const(_MPIV.pi / 6)
HALF = IV(0.5)


def _A_civ(theta, phi):
    ct, st = iv_cos(theta), iv_sin(theta)
    eip = iv_exp_i(phi)
    L11 = CIV(ct, IV(0.0))
    L12 = eip * CIV(st)
    L21 = eip.conj() * CIV(st)
    i_s3 = CIV(IV(0.0), SQRT3_HALF)
    M11 = CIV(IV(-0.5)) + i_s3 * L11
    M12 = i_s3 * L12
    M21 = i_s3 * L21
    M22 = CIV(IV(-0.5)) + i_s3 * (-L11)
    return M11 + M21, M12 + M22          # A11, A12


def _mp_sin_iv(*terms):
    """sin(sum of terms) with the SUM taken exactly at prec 100
    (floats promote exactly; PI_6-style constants passed as mp
    quantities). Collapses the eps-in-argument width of near-zero
    smalls (t1, t3) from ~4e-16 to ~1e-30."""
    acc = _MPIV.mpf(0)
    for t in terms:
        acc = acc + (t if not isinstance(t, float) else
                     _MPIV.mpf(t))
    return _iv_const(_MPIV.sin(acc))


def _mp_cos_iv(*terms):
    acc = _MPIV.mpf(0)
    for t in terms:
        acc = acc + (t if not isinstance(t, float) else
                     _MPIV.mpf(t))
    return _iv_const(_MPIV.cos(acc))


def iv_stable_z3sq(theta, phi, lam):
    """Certified enclosure of z3^2 via the factored kernel.
    Returns (z3sq_CIV, den_margin). Point-float lam/phi arguments
    get exact-argument mp evaluation for the two near-zero smalls;
    interval arguments fall back to plain interval trig."""
    A11, A12 = _A_civ(theta, phi)
    z1 = iv_exp_i(lam)
    hl = lam * IV(0.5) if hasattr(lam, "lo") else IV(lam * 0.5)
    e_hl = iv_exp_i(hl)
    lam_pt = getattr(lam, "width", 1.0) == 0.0
    phi_pt = getattr(phi, "width", 1.0) == 0.0
    if lam_pt and phi_pt:
        t1 = _mp_sin_iv(_MPIV.mpf(lam.lo) / 2, -_MPIV.pi / 6)
        t3c = _mp_cos_iv(_MPIV.mpf(lam.lo) / 2, phi.lo)
    else:
        t1 = iv_sin(hl - PI_6)
        t3c = iv_cos(hl + phi)
    sh = iv_sin(hl)
    s2t = iv_sin(theta * IV(0.5) if hasattr(theta, "lo")
                 else IV(theta * 0.5))
    t2 = -(SQRT3 * (s2t * s2t) * sh)
    t3 = SQRT3_HALF * iv_sin(theta) * t3c
    brP = CIV(t1 + t2, t3)
    brQ = CIV(t1 + t2, -t3)
    two = CIV(IV(2.0))
    P = two * e_hl * brP
    Q = two * e_hl * brQ
    bignum = A11 - z1 * A12
    bigden = A12.conj() - z1 * A11.conj()
    num = bignum * P
    den = bigden * Q
    dmag = den.abs2()
    if dmag.lo <= 0:
        raise RuntimeError("stable z3 denominator touches 0")
    return num / den, dmag.lo


def main():
    from ivkarlsson import iv_karlsson
    p3 = np.pi / 3
    print("z3^2 enclosure width, naive (iv_karlsson diag) vs "
          "stable kernel:", flush=True)
    for th in (1e-2, 1e-3, 1e-4, 1e-5, 1e-6):
        z, dm = iv_stable_z3sq(IV(th), IV(p3), IV(p3))
        ws = max(z.re.width, z.im.width)
        try:
            H, diag = iv_karlsson(IV(th), IV(p3), IV(p3))
            from ivkarlsson import h_width
            wn = f"{h_width(H):.1e}"
        except RuntimeError as e:
            wn = f"RAISE({str(e)[:24]})"
        print(f"  corner th={th:g}: naive-map {wn}  stable-z3sq "
              f"{ws:.1e}  den_margin {dm:.1e}", flush=True)


if __name__ == "__main__":
    main()
