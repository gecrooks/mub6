"""Probe pairs {I, H} across all known order-6 Hadamard families:
count MU vectors, orthonormal bases among them (=> MU triples), and MU
base-pairs (=> MU quadruples). A quadruple containing {I, H} exists iff
two of those bases are mutually unbiased.

Literature values (McNulty-Weigert review): F6(a,b): 48 vectors always
(16 bases at (0,0), 70 at (1/6,0), typically 8); Dita D(x): 120 -> 72 -> 48;
Bjorck C: 56; Tao S6: 90 (no bases at all).
The non-affine families (B, M, X, K) and generic G6^(4) points have NO
rigorous results -- this is the unexplored corner.
"""

import sys

import numpy as np

from mub import (bases_matrix, beauchamp_nicoara, bjorck_c, defect, dephase,
                 dita, find_bases, find_mu_vectors, fourier_family,
                 fourier_family_T, mu_base_pairs, random_hadamard, tao_s6,
                 unbiasedness_defect)


def probe(name, H, n_starts=3000, seed=0):
    vecs = find_mu_vectors([H], n_starts=n_starts, seed=seed)
    bases = find_bases(vecs)
    pairs = mu_base_pairs(vecs, bases)
    min_def = np.inf
    for a in range(len(bases)):
        A = bases_matrix(vecs, bases[a])
        for b in range(a + 1, len(bases)):
            B = bases_matrix(vecs, bases[b])
            min_def = min(min_def, unbiasedness_defect(A, B))
    dph = dephase(H) * np.sqrt(6)
    has_m1 = bool(np.min(np.abs(dph + 1)) < 1e-6)      # H2-reducibility marker
    print(f"{name:28s}  vecs={len(vecs):4d}  bases={len(bases):3d}  "
          f"MUpairs={len(pairs):2d}  minDefect={min_def if bases else float('nan'):.4f}  "
          f"defect(H)={defect(H)}  H2red~{has_m1}")
    sys.stdout.flush()
    return dict(name=name, n_vecs=len(vecs), n_bases=len(bases),
                n_pairs=len(pairs), min_defect=float(min_def))


def main():
    theta_min = np.arccos((np.sqrt(3) - 1) / 2)
    results = []

    print("=== affine families (rigorous results known) ===")
    results.append(probe("F(0,0) = F6", fourier_family(0, 0)))
    results.append(probe("F(1/6,0)", fourier_family(1 / 6, 0)))
    results.append(probe("F(0.1043,0.0567) generic", fourier_family(0.1043, 0.0567)))
    results.append(probe("F^T(0.1043,0.0567)", fourier_family_T(0.1043, 0.0567)))
    results.append(probe("D(0) Dita", dita(0)))
    results.append(probe("D(1/8) endpoint", dita(1 / 8)))
    results.append(probe("D(0.05)", dita(0.05)))
    results.append(probe("C Bjorck", bjorck_c()))
    results.append(probe("S6 Tao", tao_s6()))

    print("=== non-affine Hermitian family B(theta) (no rigorous results) ===")
    for th in [theta_min + 0.15, theta_min + 0.6, np.pi / 2 + 0.5, np.pi - 0.4]:
        results.append(probe(f"B({th:.3f})", beauchamp_nicoara(th)))

    print("=== generic random Hadamards (G6^(4) territory, unexplored) ===")
    for k in range(12):
        H = random_hadamard(6, seed=1000 + k)
        results.append(probe(f"random #{k}", H, seed=k))

    np.save("families_results.npy", results, allow_pickle=True)
    print("\nsaved -> families_results.npy")


if __name__ == "__main__":
    main()
