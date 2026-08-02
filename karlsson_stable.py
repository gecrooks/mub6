"""Cancellation-free evaluation of the Karlsson map near the
(pi/3, pi/3) corner (NOTES_LP_BRIDGE 4.56-4.57).

The naive Moebius quotient for z3^2 loses eps/theta^3 at the
corner because numerator and denominator share vanishing factors.
Exact factorization (derived 2026-08-01):

  num = A11^2 - w A12^2 = (A11 - z1 A12)(A11 + z1 A12)
  den = cA12^2 - w cA11^2 = (cA12 - z1 cA11)(cA12 + z1 cA11)

with c* = conjugate, w = z1^2, and the SMALL factors given exactly
by the trig identity

  P = A11 + z1 A12
    = 2 e^{i lam/2} [ sin(lam/2 - pi/6)
                      - sqrt3 sin^2(theta/2) sin(lam/2)
                      + i (sqrt3/2) sin(theta) cos(lam/2 + phi) ]
  Q = cA12 + z1 cA11 = same bracket with the imaginary term
      negated.

Each bracket term is individually small near the special loci and
carries them algebraically: sin(lam/2 - pi/6) = 0 IS the wall
lam = pi/3; cos(lam/2 + phi) = 0 IS the 2 phi + lam = pi locus
through the corner; the theta^2 term is the map's even collapse.
The big factors never vanish there (checked: |A11 - z1 A12| = 2
at the corner).

Measured accuracy of stable_z3sq vs a 40-digit reference:
8.9e-12 relative at theta = 1e-4 (naive 1.5e-8), 8.9e-10 at 1e-6
(naive 1.3e-4), 1.1e-8 at 1e-8 (naive O(1)). z4^2 has no corner
degeneracy (its factors are O(1) there); z2^2 inherits stability
through z3^2 (residual relative error ~ 1e-14/|den2|, sufficient
to theta ~ 1e-6; the exact-delta form is the documented
completion).
"""

import numpy as np

SQRT3 = np.sqrt(3.0)


def _A(theta, phi):
    c, s = np.cos(theta), np.sin(theta)
    A11 = -0.5 + 1j * (SQRT3 / 2) * (c + np.exp(-1j * phi) * s)
    A12 = -0.5 + 1j * (SQRT3 / 2) * (np.exp(1j * phi) * s - c)
    return A11, A12


def stable_PQ(theta, phi, lam):
    hl = lam / 2
    e = np.exp(1j * hl)
    t1 = np.sin(hl - np.pi / 6)
    t2 = -SQRT3 * np.sin(theta / 2) ** 2 * np.sin(hl)
    t3 = (SQRT3 / 2) * np.sin(theta) * np.cos(hl + phi)
    return 2 * e * (t1 + t2 + 1j * t3), 2 * e * (t1 + t2 - 1j * t3)


def stable_z3sq(theta, phi, lam):
    A11, A12 = _A(theta, phi)
    z1 = np.exp(1j * lam)
    P, Q = stable_PQ(theta, phi, lam)
    return ((A11 - z1 * A12) * P) / ((np.conj(A12)
                                      - z1 * np.conj(A11)) * Q)


def _mp_z2(theta, phi, lam):
    from mpmath import mp, mpc, sqrt as msqrt, exp as mexp, \
        cos as mcos, sin as msin
    with mp.workdps(30):
        I = mpc(0, 1)
        c, s = mcos(theta), msin(theta)
        A11 = mp.mpf(-1) / 2 + I * (msqrt(3) / 2) * (
            c + mexp(-I * phi) * s)
        A12 = mp.mpf(-1) / 2 + I * (msqrt(3) / 2) * (
            mexp(I * phi) * s - c)
        B11, B12 = -1 - A11, -1 - A12
        w = mexp(2 * I * lam)
        cj = lambda z: z.conjugate()
        z3sq = (A11 ** 2 - w * A12 ** 2) / (cj(A12) ** 2
                                            - w * cj(A11) ** 2)
        z2 = msqrt((B11 ** 2 - z3sq * cj(B12) ** 2)
                   / (B12 ** 2 - z3sq * cj(B11) ** 2))
        return complex(z2)


def stable_karlsson(theta, phi, lam):
    """Karlsson map with the stable z3^2 kernel; z4, z2 as in the
    naive map (z4 has no corner degeneracy; z2 inherits stability
    through z3sq)."""
    F2 = np.array([[1, 1], [1, -1]], dtype=complex)
    c, s = np.cos(theta), np.sin(theta)
    L = np.array([[c, np.exp(1j * phi) * s],
                  [np.exp(-1j * phi) * s, -c]], dtype=complex)
    A = F2 @ (-0.5 * np.eye(2) + 1j * (SQRT3 / 2) * L)
    B = -F2 - A
    B11, B12 = B[0, 0], B[0, 1]
    z1 = np.exp(1j * lam)
    w = z1 ** 2
    z3sq = stable_z3sq(theta, phi, lam)
    z4sq = (B11 ** 2 - w * B12 ** 2) / (np.conj(B12) ** 2
                                        - w * np.conj(B11) ** 2)
    z3, z4 = np.sqrt(z3sq), np.sqrt(z4sq)
    den2 = B12 ** 2 - z3sq * np.conj(B11) ** 2
    if abs(den2) < 1e-3:
        # corner basin: the fully factored stable_z2sq (4.65/4.67,
        # exact odd-series small factors). Residual float floor:
        # eps in the near-zero trig ARGUMENTS amplified by 1/s^2 —
        # keep the mp path only in the deepest sliver.
        if abs(np.sin(theta)) < 3e-5:
            z2 = _mp_z2(theta, phi, lam)
        else:
            z2 = np.sqrt(stable_z2sq(theta, phi, lam))
    else:
        z2 = np.sqrt((B11 ** 2 - z3sq * np.conj(B12) ** 2) / den2)
    Z1 = np.array([[1, 1], [z1, -z1]], dtype=complex)
    Z2 = np.array([[1, 1], [z2, -z2]], dtype=complex)
    Z3 = np.array([[1, z3], [1, -z3]], dtype=complex)
    Z4 = np.array([[1, z4], [1, -z4]], dtype=complex)
    H = np.zeros((6, 6), dtype=complex)
    H[0:2, 0:2] = F2
    H[0:2, 2:4] = Z1
    H[0:2, 4:6] = Z2
    H[2:4, 0:2] = Z3
    H[4:6, 0:2] = Z4
    H[2:4, 2:4] = 0.5 * Z3 @ A @ Z1
    H[2:4, 4:6] = 0.5 * Z3 @ B @ Z2
    H[4:6, 2:4] = 0.5 * Z4 @ B @ Z1
    H[4:6, 4:6] = 0.5 * Z4 @ A @ Z2
    return H / np.sqrt(6)

W3 = np.exp(2j * np.pi / 3)


def _g1(phi, lam):
    z1 = np.exp(1j * lam)
    return 1j * (SQRT3 / 2) * (
        -4 * W3 * np.exp(1j * phi) * (W3 - z1 * np.conj(W3))
        - 2 * W3 ** 2 * (np.exp(-1j * phi) - z1 * np.exp(1j * phi)))


def _g21(phi, lam):
    z1 = np.exp(1j * lam)
    mu = 1j * (SQRT3 / 2) * np.exp(-1j * phi)
    nu = 1j * (SQRT3 / 2) * np.exp(1j * phi)
    n0 = W3 - z1 * np.conj(W3)
    n1 = mu - z1 * nu
    return -4 * np.conj(W3) * mu * n0 - 2 * np.conj(W3) ** 2 * n1


def _G_pair(s, phi, lam):
    """(G-, G+, G2-, G2+) at sin(theta) = s, naive evaluation
    (used at the safe extraction point s0 and for the O(1)
    plus-combinations)."""
    c = np.sqrt(max(1.0 - s * s, 0.0))
    A11 = -0.5 + 1j * (SQRT3 / 2) * (c + np.exp(-1j * phi) * s)
    A12 = -0.5 + 1j * (SQRT3 / 2) * (np.exp(1j * phi) * s - c)
    B11, B12 = -1 - A11, -1 - A12
    z1 = np.exp(1j * lam)
    bignum = A11 - z1 * A12
    bigden = np.conj(A12) - z1 * np.conj(A11)
    Gm = B12 ** 2 * bigden - np.conj(B11) ** 2 * bignum
    Gp = B12 ** 2 * bigden + np.conj(B11) ** 2 * bignum
    G2m = B11 ** 2 * bigden - np.conj(B12) ** 2 * bignum
    G2p = B11 ** 2 * bigden + np.conj(B12) ** 2 * bignum
    return Gm, Gp, G2m, G2p


def stable_z2sq(theta, phi, lam):
    """Fully factored z2^2 (NOTES 4.65/4.67). The common bigden*Q
    cancels, and the shared small factor s = sin(theta) is
    cancelled ANALYTICALLY (tau3 = t3/s is exact trig):
      z2^2 = (t ghat2 - i tau3 G2+) / (t ghat - i tau3 G+),
    ghat = g1 + s^2 g3 (G- exactly odd, g3 extracted at s0 = 0.02,
    exact to O(s^4) relative), t = t1 + t2. For s > 0.05 the naive
    quotient is well-conditioned and used directly."""
    s = np.sin(theta)
    if abs(s) > 0.05:
        Gm, Gp, G2m, G2p = _G_pair(s, phi, lam)
        hl = lam / 2
        t = np.sin(hl - np.pi / 6) \
            - SQRT3 * np.sin(theta / 2) ** 2 * np.sin(hl)
        t3 = (SQRT3 / 2) * s * np.cos(hl + phi)
        return (t * G2m - 1j * t3 * G2p) / (t * Gm - 1j * t3 * Gp)
    hl = lam / 2
    t = np.sin(hl - np.pi / 6) \
        - SQRT3 * np.sin(theta / 2) ** 2 * np.sin(hl)
    tau3 = (SQRT3 / 2) * np.cos(hl + phi)
    s0 = 0.02
    Gm0, _, G2m0, _ = _G_pair(s0, phi, lam)
    g1 = _g1(phi, lam)
    g21 = _g21(phi, lam)
    g3 = (Gm0 - s0 * g1) / s0 ** 3
    g23 = (G2m0 - s0 * g21) / s0 ** 3
    ghat = g1 + s * s * g3
    ghat2 = g21 + s * s * g23
    _, Gp, _, G2p = _G_pair(s, phi, lam)
    num = t * ghat2 - 1j * tau3 * G2p
    den = t * ghat - 1j * tau3 * Gp
    return num / den
