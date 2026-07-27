"""High-precision analysis of degenerate MU-vector solutions.

Q1: For {I, C_Bjorck}: are the two Jacobian-rank-3 solutions exact double
    roots (isolated, multiplicity >= 2) or pairs of nearby simple roots?
Q2: For {I, F(1/6,0)}: do solutions sit on positive-dimensional continua?

Method: residual growth along the numerical null direction. For an isolated
multiplicity-2 root, moving by t along the null direction grows the residual
~ c*t^2. For a genuine continuum, the residual stays ~ 0 to all orders
(numerically: growth like t^4 or below noise).
"""

import numpy as np

from mub import (bjorck_c, find_mu_vectors, fourier_family,
                 mu_vector_residuals)


def null_directions(v, mats, rank_tol=1e-4):
    th = np.angle(v * np.sqrt(6))[1:]
    eps = 1e-7
    f0 = mu_vector_residuals(th, mats)
    J = np.zeros((len(f0), 5))
    for k in range(5):
        t = th.copy()
        t[k] += eps
        J[:, k] = (mu_vector_residuals(t, mats) - f0) / eps
    U, s, Vt = np.linalg.svd(J)
    null = Vt[s.shape[0]:] if J.shape[0] < 5 else Vt[np.sum(s > rank_tol * s[0]):]
    return th, null, s


def residual_growth(th, direction, mats, ts=(1e-4, 1e-3, 1e-2)):
    out = []
    for t in ts:
        r = mu_vector_residuals(th + t * direction, mats)
        out.append(np.max(np.abs(r)))
    return out


def analyse(name, H, vecs):
    print(f"--- {name}: {len(vecs)} solutions ---")
    n_deg = 0
    for i, v in enumerate(vecs):
        th, null, s = null_directions(v, [H])
        if len(null) == 0:
            continue
        n_deg += 1
        growth = residual_growth(th, null[0], [H])
        # fit exponent between t=1e-3 and 1e-2
        expo = np.log10(growth[2] / growth[1])
        kind = ("continuum?" if growth[2] < 1e-10 else
                f"isolated deg. root (residual ~ t^{expo:.1f})")
        print(f"  sol {i:3d}: J-singvals {np.round(s, 6)}  "
              f"growth {['%.1e' % g for g in growth]} -> {kind}")
    if n_deg == 0:
        print("  no degenerate solutions (all Jacobians full rank)")


def main():
    C = bjorck_c()
    vecs_c = find_mu_vectors([C], n_starts=40000, seed=11)
    # dedup harder for analysis clarity
    print(f"Bjorck raw count {len(vecs_c)}")
    analyse("C Bjorck", C, vecs_c)

    F16 = fourier_family(1 / 6, 0)
    vecs_f = find_mu_vectors([F16], n_starts=30000, seed=17)
    print(f"\nF(1/6,0) raw count {len(vecs_f)}")
    analyse("F(1/6,0)", F16, vecs_f)


if __name__ == "__main__":
    main()
