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

SQRT3_HALF = IV.pad(np.sqrt(3.0) / 2.0, 8)
SQRT3 = IV.pad(np.sqrt(3.0), 8)
PI_6 = IV.pad(np.pi / 6.0, 8)
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


def iv_stable_z3sq(theta, phi, lam):
    """Certified enclosure of z3^2 via the factored kernel.
    Returns (z3sq_CIV, den_margin)."""
    A11, A12 = _A_civ(theta, phi)
    z1 = iv_exp_i(lam)
    hl = lam * IV(0.5) if hasattr(lam, "lo") else IV(lam * 0.5)
    e_hl = iv_exp_i(hl)
    t1 = iv_sin(hl - PI_6)
    sh = iv_sin(hl)
    s2t = iv_sin(theta * IV(0.5) if hasattr(theta, "lo")
                 else IV(theta * 0.5))
    t2 = -(SQRT3 * (s2t * s2t) * sh)
    t3 = SQRT3_HALF * iv_sin(theta) * iv_cos(hl + phi)
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
