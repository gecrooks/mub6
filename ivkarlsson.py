"""Certified interval enclosures of the Karlsson map H(theta, phi, lam)
over parameter boxes — the R4 substrate replacing the EMPIRICAL sampled
map constants of the prototype certificates.

Mirrors karlsson.karlsson_map with rectangle complex intervals; the
Moebius/sqrt construction's branch choices become explicit certified
obligations: every csqrt argument must clear the principal branch cut,
and every divisor's |.|^2 must exclude 0 — verified per box, else raise
(the campaign response is tile subdivision).

Validation harness:
  1. containment: float map inside point-width enclosures (~1e-15);
  2. tile enclosures vs the prototype's PAD*FD Lipschitz constants;
  3. branch/denominator headroom across the 24 scan points;
  4. certified g-evaluation widths vs the SLOP model of certify.py.
"""

import warnings

import numpy as np

from interval import CIV, IV, civ_matmul, iv_cos, iv_exp_i, iv_sin
from karlsson import karlsson_map

warnings.filterwarnings("ignore")

SQRT3_HALF = IV.pad(np.sqrt(3.0) / 2.0, 8)
INV_SQRT6 = IV.pad(1.0 / np.sqrt(6.0), 8)


def iv_karlsson(theta, phi, lam):
    """Interval Karlsson map. Parameters are IV; returns (H, diag) where
    H is a 6x6 list-of-lists of CIV (entries include the 1/sqrt6) and
    diag records branch/denominator margins."""
    ct, st = iv_cos(theta), iv_sin(theta)
    eip = iv_exp_i(phi)
    L11 = CIV(ct, IV(0.0))
    L12 = eip * CIV(st)
    L21 = eip.conj() * CIV(st)
    L22 = -L11
    # A = F2 (-e/2 + i sqrt3/2 L), F2 = [[1,1],[1,-1]]
    i_s3 = CIV(IV(0.0), SQRT3_HALF)
    M11 = CIV(IV(-0.5)) + i_s3 * L11
    M12 = i_s3 * L12
    M21 = i_s3 * L21
    M22 = CIV(IV(-0.5)) + i_s3 * L22
    A11 = M11 + M21
    A12 = M12 + M22
    B11 = -CIV(1.0) - A11
    B12 = -CIV(1.0) - A12

    z1 = iv_exp_i(lam)
    w = z1 * z1

    def moebius(P11, P12):
        num = P11 * P11 - w * (P12 * P12)
        den = (P12.conj() * P12.conj()) - w * (P11.conj() * P11.conj())
        dmag = den.abs2()
        if dmag.lo <= 0:
            raise RuntimeError("Moebius denominator interval touches 0")
        return num / den, dmag.lo

    # z3 via the FACTORED kernel (NOTES 4.57): algebraically
    # identical to the Moebius quotient but cancellation-free —
    # removes the eps/theta^3 interval blowup at the (pi/3, pi/3)
    # corner. z4's factors are O(1) there; it keeps the plain
    # quotient.
    from ivstable import iv_stable_z3sq
    z3sq, d3 = iv_stable_z3sq(theta, phi, lam)
    z4sq, d4 = moebius(B11, B12)
    num2 = B11 * B11 - z3sq * (B12.conj() * B12.conj())
    den2 = B12 * B12 - z3sq * (B11.conj() * B11.conj())
    d2m = den2.abs2()
    if d2m.lo <= 0:
        raise RuntimeError("z2 denominator interval touches 0")
    z2sq = num2 / den2

    cuts = {}
    for name, zsq in (("z3", z3sq), ("z4", z4sq), ("z2", z2sq)):
        if not zsq.cut_clear():
            raise RuntimeError(f"{name}^2 rectangle touches the branch cut")
        # cut margin: distance of the rectangle from (-inf, 0]
        m = zsq.re.lo if zsq.re.lo > 0 else min(abs(zsq.im.lo),
                                                abs(zsq.im.hi))
        cuts[name] = m
    z3 = z3sq.csqrt()
    z4 = z4sq.csqrt()
    z2 = z2sq.csqrt()

    one = CIV(1.0)
    F2 = [[one, one], [one, -one]]
    Z1 = [[one, one], [z1, -z1]]
    Z2 = [[one, one], [z2, -z2]]
    Z3 = [[one, z3], [one, -z3]]
    Z4 = [[one, z4], [one, -z4]]
    A = [[A11, A12], [A12.conj(), -A11.conj()]]
    B = [[B11, B12], [B12.conj(), -B11.conj()]]

    half = IV(0.5)
    blocks = {}
    blocks[(0, 0)] = F2
    blocks[(0, 1)] = Z1
    blocks[(0, 2)] = Z2
    blocks[(1, 0)] = Z3
    blocks[(2, 0)] = Z4
    blocks[(1, 1)] = [[e.scale(half) for e in row]
                      for row in civ_matmul(civ_matmul(Z3, A), Z1)]
    blocks[(1, 2)] = [[e.scale(half) for e in row]
                      for row in civ_matmul(civ_matmul(Z3, B), Z2)]
    blocks[(2, 1)] = [[e.scale(half) for e in row]
                      for row in civ_matmul(civ_matmul(Z4, B), Z1)]
    blocks[(2, 2)] = [[e.scale(half) for e in row]
                      for row in civ_matmul(civ_matmul(Z4, A), Z2)]

    H = [[None] * 6 for _ in range(6)]
    for (bi, bj), blk in blocks.items():
        for a in range(2):
            for b in range(2):
                H[2 * bi + a][2 * bj + b] = blk[a][b].scale(INV_SQRT6)
    diag = dict(den_margins=(d3, d4, d2m.lo), cut_margins=cuts)
    return H, diag


def h_width(H):
    return max(H[i][j].width for i in range(6) for j in range(6))


def main():
    print("=== certified Karlsson map: validation harness ===")

    # 1. containment at points
    rng = np.random.default_rng(4)
    worst_w, n_ok = 0.0, 0
    for _ in range(200):
        th, ph, la = rng.uniform(0.25, 2 * np.pi - 0.25, 3)
        Hf = karlsson_map(th, ph, la)
        try:
            Hi, _ = iv_karlsson(IV(th), IV(ph), IV(la))
        except RuntimeError:
            continue
        ok = all(Hi[i][j].contains(Hf[i, j])
                 for i in range(6) for j in range(6))
        n_ok += ok
        worst_w = max(worst_w, h_width(Hi))
        if not ok:
            print(f"  CONTAINMENT FAILURE at {th, ph, la}")
    print(f"1. point containment: {n_ok}/200 verified inside enclosures; "
          f"max enclosure width {worst_w:.2e}")

    # 2. tile enclosure vs prototype constants at the reference point
    beta = (5.978503016422594, 4.007534549834652, 1.6327649325136653)
    Hc = karlsson_map(*beta)
    for h in (3e-4, 1e-3, 3e-3):
        Hi, diag = iv_karlsson(IV(beta[0] - h, beta[0] + h),
                               IV(beta[1] - h, beta[1] + h),
                               IV(beta[2] - h, beta[2] + h))
        # certified bound on max |H(b) - H(center)| over the tile
        dev = max(max(abs(Hi[i][j].re.lo - Hc[i, j].real),
                      abs(Hi[i][j].re.hi - Hc[i, j].real),
                      abs(Hi[i][j].im.lo - Hc[i, j].imag),
                      abs(Hi[i][j].im.hi - Hc[i, j].imag))
                  for i in range(6) for j in range(6))
        proto = 2.0 * 1.04 * np.sqrt(3) * h        # PAD * L_map * sqrt3 h
        print(f"2. tile h={h:g}: certified max|dH| = {dev:.3e}  "
              f"(prototype PAD*FD bound {proto:.3e}; "
              f"ratio {dev/proto:.2f})")

    # 3. branch/denominator headroom across the scan points
    rng = np.random.default_rng(20260726)
    dmins, cmins, fails = [], [], 0
    for t in range(24):
        b = rng.uniform(0.25, 2 * np.pi - 0.25, 3)
        try:
            _, diag = iv_karlsson(IV(b[0] - 3e-4, b[0] + 3e-4),
                                  IV(b[1] - 3e-4, b[1] + 3e-4),
                                  IV(b[2] - 3e-4, b[2] + 3e-4))
            dmins.append(min(diag["den_margins"]))
            cmins.append(min(diag["cut_margins"].values()))
        except RuntimeError:
            fails += 1
    print(f"3. h=3e-4 boxes at 24 scan points: {fails} branch/denominator "
          f"failures; min |den|^2 {min(dmins):.3f}; "
          f"min cut margin {min(cmins):.3f}")

    # 4. certified g widths vs the SLOP model
    from certify import _uvec
    Hi, _ = iv_karlsson(IV(beta[0]), IV(beta[1]), IV(beta[2]))
    worst_g = 0.0
    for trial in range(50):
        th = rng.uniform(0, 2 * np.pi, 5)
        u = _uvec(th[None, :])[0]
        ui = [CIV.from_float(complex(x)) for x in u]
        for k in range(1, 6):
            s = CIV(0.0)
            for j in range(6):
                s = s + Hi[j][k].conj() * ui[j]
            gk = s.abs2() - IV.pad(1.0 / 6.0, 4)
            worst_g = max(worst_g, gk.width)
    print(f"4. certified g-evaluation width (point H, float theta): "
          f"max {worst_g:.2e}  (SLOP model assumed 1e-11)")


if __name__ == "__main__":
    main()
